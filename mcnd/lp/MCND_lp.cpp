// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#include <cmath>
#include <algorithm>
#include <numeric>

#include "BCP_math.hpp"
#include "BCP_lp.hpp"
#include "BCP_lp_node.hpp"

#include "MCND_lp.hpp"

//#############################################################################

void
MCND_lp::unpack_module_data(BCP_buffer & buf){
    //std::cout<<"try unpack to lp "<<std::endl;
    data.unpack(buf); 
    bool has_init_sol;
    buf.unpack(has_init_sol);
    if(has_init_sol) best_sol.unpack(buf);

    y.resize(data.narcs,0);
    yfix = new int[data.narcs];
    std::fill(yfix,yfix+data.narcs, -1);
    ninsp.resize(data.narcs,Pair(0,0));
    psdcost.resize(data.narcs,PairF(0,0));
    
    ss_manager.initialize(&data);
    cover_manager.initialize(&data, 20);
    pump_heur.set_data(&data, 0.5);
    lpfeaschecker.initialize(&data);
    flwconnect.initialize(&data);
    localc_manager.initialize(&data, 50);
    globalc_manager.initialize(&data, 50);
    
    AppVolData.data = &data;
    AppVolData.ss_manager = &ss_manager;
    AppVolData.cover_manager = &cover_manager;
    AppVolData.localc_manager = &localc_manager;
    AppVolData.globalc_manager = &globalc_manager;

    srand(10);
     
}

//-------------------------------------------------------------------------------------------

OsiSolverInterface *
MCND_lp::initialize_solver_interface(){
    OsiVolSolverInterface* volsolver = new OsiVolSolverInterface("volmcnd.par");
    setOsiBabSolver(volsolver);
    //MaxIt = volsolver->volprob_.parm.maxsgriters;
    //std::cout<<"vol instantiated "<<MaxIt<<std::endl;
    return volsolver;
    
}

//#############################################################################
//#############################################################################
//#############################################################################
// Override the initializer so that we can choose between vol and simplex
// at runtime.
void
MCND_lp::initialize_new_search_tree_node(const BCP_vec<BCP_var*>& vars,
                                         const BCP_vec<BCP_cut*>& cuts,
                                         const BCP_vec<BCP_obj_status>& var_status,
                                         const BCP_vec<BCP_obj_status>& cut_status,
                                         BCP_vec<int>& var_changed_pos,
                                         BCP_vec<double>& var_new_bd,
                                         BCP_vec<int>& cut_changed_pos,
                                         BCP_vec<double>& cut_new_bd){
    
    //std::cout<<"initialize_new_search_tree_node "<<cuts.size()<<std::endl;
    ++num_nodes;
    lp_mode = LP_Normal;
    std::fill(yfix,yfix+data.narcs, -1);
    has_sol = getLpProblemPointer()->has_ub();
   
    bool testconn=false;
    MCND_node_branch_data* nodedata = dynamic_cast<MCND_node_branch_data*>( get_user_data());
    if(nodedata!=0){
        //std::cout<<"user_data hotstart: "<<nodedata->hs<<std::endl;
        LBi = nodedata->min_lb;
        if(nodedata->hs!=0){
            getLpProblemPointer()->lp_solver->setWarmStart(nodedata->hs);
        }
        //std::cout<<"node comes from branch on: "<<nodedata->branch_var;
        //std::cout<<" of "<<nodedata->pos_neg<<" side .. parent: "<<nodedata->parent<<std::endl;
        testconn = nodedata->test_conn ;
        if(nodedata->reduced_run)lp_mode |= LP_ReducedRun;
        //std::cout<<"reduced run? "<<nodedata->reduced_run<<std::endl;

    }else LBi=0;
    if(testconn ){
    	if(!flwconnect.check_connectivity(vars, var_changed_pos, var_new_bd, testconn, yfix)){
			var_changed_pos.clear();
			var_new_bd.clear();
			//std::cout<<"fathom node flwconnect.check_connectivity"<<std::endl;
			lp_mode |= LP_ForceNodeAbort;
			return;
		}	
    } 
    
    cut_varfix_and_updt( vars, cuts, var_changed_pos, var_new_bd);
    /*
    for (int a=data.narcs; a--;){
		if(vars[a]->lb()==1.0 || yfix[a]==1) std::cout<<"fixed "<<a<<" to 1"<<std::endl;
		else if(vars[a]->ub()==0.0 || yfix[a]==0) std::cout<<"fixed "<<a<<" to 0"<<std::endl;
	}*/
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
                      BCP_var_set& vars, BCP_cut_set& cuts){
	//if(lp_mode & LP_ForceNodeAbort) return;
    //std::cout<<"load problem "<<core->varnum()<<" "<<cuts.size()<<std::endl;
    cover_manager.clean_collection();
    localc_manager.clean_collection();
	globalc_manager.clean_collection();
	
    int cutnum = cuts.size();
  	int core_rownum = core->cutnum();
  	
  	if (cutnum > core_rownum){
  		for(int i=core_rownum; i<cutnum;++i){
  			MCND_Cut * mcnd_cut = dynamic_cast<MCND_Cut*>(cuts[i]);
  			switch(mcnd_cut->type){
  				case 1:{
					CoverCut * cut = dynamic_cast<CoverCut*>(mcnd_cut);
					//std::cout<<"collec: "<<i<<" "<<cut->get_cover()->serial_nmbr<<std::endl;
					cut->get_cover()->id_vi = i;
					cover_manager.covers.insert_end(cut->get_cover());
					break;
  				}case 2:{
					LocalCCut * lcut = dynamic_cast<LocalCCut*>(mcnd_cut);
					lcut->get_localc()->id_vi = i;
					localc_manager.locals.insert_end(lcut->get_localc());
					//std::cout<<"colleclc: "<<i<<" "<<lcut->get_localc()->serial_nmbr<<" sz: "<<localc_manager.locals.sizeOfCollection<<std::endl;
					break; 
				}
  			}
		}
  	}
  	globalc_manager.collect_globals(vars, cutnum);
  	
    osi.setApplicationData( &AppVolData);
    osi.loadProblem(core->varnum(), core->cutnum(),
                    0, 0, /*const int* start, const int* index*/
                    0, /*const double* value*/
                    core->matrix->ColLowerBound().begin(),core->matrix->ColUpperBound().begin(), /*const double* collb, const double* colub*/
                    0, /*const double* obj*/
                    core->matrix->RowLowerBound().begin(),core->matrix->RowUpperBound().begin() /*const double* rowlb, const double* rowub*/);
    
    //std::cout<<"ok "<<topo->sznz<<" "<<topo->szunfxd<<std::endl;
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::modify_lp_parameters(OsiSolverInterface* lp, const int changeType,
                              bool in_strong_branching){
    //std::cout<<"modify lp param "<<changeType<<" "<<in_strong_branching<<" iter: "<<current_level()<<std::endl;
    lp->setDblParam(OsiPrimalTolerance, 1e-4);
    OsiVolSolverInterface* vollp = getOsiVolBabSolver();
    VOL_parms& par = vollp->volprob_.parm;
    
    vollp->has_sol = has_sol;
    vollp->recheck_collct=false;
    vollp->upper_bound = LBi >0 ? fmin(upper_bound(), LBi*5) : upper_bound();
    par.dual_limit = vollp->upper_bound*5; 
    //std::cout<<"limit: "<<par.dual_limit<<" "<<vollp->upper_bound<<" strongbranching:"<<in_strong_branching<<std::endl;

    if(lp_mode & LP_ReducedRun){
    	 par.maxsgriters = 100;
    	 lp_mode &= ~LP_ReducedRun;
    }else if(!vollp->HotStartSet){
    	if(current_level()==0) par.maxsgriters = 1000;
    	else par.maxsgriters = 500;
    }else par.maxsgriters = 250;
     
    if(current_level()>0){ 
		par.lambdainit=1.0;
		par.alphainit=0.5;
		par.alphamin=1e-2;
		par.alphafactor=0.3;
		par.alphaint=10;

     	AppVolData.intvlVI = 200; /*vollp->mode=2;*/
    }
	else{ vollp->mode=1;}
	
	if(lp_mode & LP_LogicalFixed){
		vollp->recheck_collct=true;  
		lp_mode &= ~LP_LogicalFixed; /*vollp->mode=3; par.maxsgriters = 100;*/
	}
    //else
    if(in_strong_branching){
    	vollp->recheck_collct=true; 
    	vollp->min_lower_bound = LBi;  
    	vollp->mode=0;
	}else if(changeType==1 || lp_mode & LP_ForceNodeAbort ){ 
    	vollp->recheck_collct=true; 
    	vollp->mode=-1; 
    	//std::cout<<"changeType "<<changeType<<" lpmode: "<<(lp_mode & LP_ForceNodeAbort)<<std::endl; 
    }else if((lp_mode & LP_CutAddedFromHeuristic) ){
    	vollp->mode=-1; 
    	//std::cout<<"LP_CutAddedFromHeuristic"<<std::endl;
    }
    vollp->in_strong_branch = in_strong_branching;

   
}

//#############################################################################
//#############################################################################
//#############################################################################

void
MCND_lp::generate_cuts_in_lp(const BCP_lp_result& lpres,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts,
                             BCP_vec<BCP_cut*>& new_cuts,
                             BCP_vec<BCP_row*>& new_rows){
    
    //std::cout<<"generate_cuts_in_lp "<<cuts.size()<<std::endl;
	if((lp_mode & LP_StrongBranch)||
		lp_mode & LP_ForceNodeAbort) return;
	
    int sz = cover_manager.covers.sizeOfCollection;
 	Cover *vi = cover_manager.covers.begin;
 	LocalCut * lcvi = localc_manager.locals.begin;
    if(lp_mode & LP_ForceNodeAbort){
    	if(cover_manager.gend>0){
    		cover_manager.gend =0;
    		for(int i=sz;i--;){
				if(vi->toadd){ keep_track(vi); vi->toadd=false;}
				vi = vi->next;
			}
    	} 
    	if(localc_manager.gend>0){
    		localc_manager.gend =0;
			sz = localc_manager.locals.sizeOfCollection;
			for(int i=sz;i--;){
				if(lcvi->toadd){ keep_track(lcvi); lcvi->toadd=false;}
				lcvi = lcvi->next;
			}
		}
     	return;
    }
    int newly_added =0;
    
    if(	lp_mode & LP_tighterBounds){
        getOsiVolBabSolver()->add_external_localc(yfix);
		lp_mode &= ~LP_tighterBounds;
    }
    
    if(cover_manager.gend==0 && localc_manager.gend==0) return;
	
	if(cover_manager.gend>0){
		vi = cover_manager.covers.begin;
		for(int i=sz;i--;){
			if(vi->toadd){ 
				keep_track(vi);
				new_cuts.push_back(new CoverCut(vi)); 
				vi->toadd=false;
				++newly_added;
			}
			vi = vi->next;
		}
	}
	
	if(localc_manager.gend>0){
		lcvi = localc_manager.locals.begin;
		sz = localc_manager.locals.sizeOfCollection;
		for(int i=sz;i--;){
			if(lcvi->toadd){ 
				keep_track(lcvi);
				//std::cout<<"generate lc "<<lcvi->id_vi<<" "<<lcvi->serial_nmbr<<std::endl;
				new_cuts.push_back(new LocalCCut(lcvi)); 
				lcvi->toadd=false;
				++newly_added;
			}
			lcvi = lcvi->next;
		}
	}
    
    //std::cout<<"generate cuts: "<<newly_added<<" cutsvecsz: "<<new_cuts.size()<<" "<<sz<<std::endl;
 
    
    if(lp_mode & LP_CutAddedFromHeuristic) lp_mode &= ~LP_CutAddedFromHeuristic;
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::cuts_to_rows(const BCP_vec<BCP_var*>& vars,
                      BCP_vec<BCP_cut*>& cuts,
                      BCP_vec<BCP_row*>& rows,
                      const BCP_lp_result& lpres,
                      BCP_object_origin origin, bool allow_multiple){
    
    const int cutnum = cuts.size();
    //std::cout<<"cuts_to_rows: "<<cutnum<<std::endl;
    for (int i=cutnum; i--;) {
        BCP_row * row = new BCP_row();
        row->LowerBound(1e10);
        row->UpperBound(1e10);
        rows.push_back(row);
    }
}

//==========================================================================================

void
MCND_lp::select_vars_to_delete(const BCP_lp_result& lpres,
				   const BCP_vec<BCP_var*>& vars,
				   const BCP_vec<BCP_cut*>& cuts,
				   const bool before_fathom,
				   BCP_vec<int>& deletable){
	//std::cout<<"slect vars"<<std::endl;
	if(lp_mode & LP_ForceNodeAbort) return;
    BCP_lp_user::select_vars_to_delete(lpres,vars,cuts,before_fathom, deletable);
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::select_cuts_to_delete(const BCP_lp_result& lpres,
				   const BCP_vec<BCP_var*>& vars,
				   const BCP_vec<BCP_cut*>& cuts,
				   const bool before_fathom,
				   BCP_vec<int>& deletable){
	//std::cout<<"select_cuts_to_delete "<<before_fathom<<std::endl;
	if(lp_mode & LP_ForceNodeAbort) return;
	//if(!before_fathom){
	const int cutnum = cuts.size();
	deletable.reserve(cutnum-getLpProblemPointer()->core->cutnum()+1);

	for (int i = getLpProblemPointer()->core->cutnum(); i < cutnum; ++i) {
		MCND_Cut * cut = dynamic_cast<MCND_Cut*>(cuts[i]);
		//std::cout<<"cut: "<<cut->type<<" "<<cut->id_vi()<<" "<<cut->purgbl()<<std::endl;
 		if ( cut->purgbl() /*|| (cut->check_viol(lpres.x())==0 && lpres.pi()[cut->id_vi()]==0)*/) {
			//std::cout<<"out / dual: "<<" id: "<<cut->id_vi()<<" srnb: "<<cut->serial_nmbr()<<std::endl;
			switch(cut->type){
				case 1:{ cover_manager.purgbl.push_back((dynamic_cast<CoverCut *>(cut))->get_cover()); break;}
				case 2:{ localc_manager.purgbl.push_back((dynamic_cast<LocalCCut *>(cut))->get_localc()); break;}
			}
			deletable.unchecked_push_back(i);
			cut->mark_unpurgbl();
		}
	}
	if(getOsiVolBabSolver()->deleteGlobalRows())
		if(deletable.size()==0) getOsiVolBabSolver()->rebuild_collections();
	getOsiVolBabSolver()->num_purgbl=0;

}

//-------------------------------------------------------------------------------------------

void
MCND_lp::purge_slack_pool(const BCP_vec<BCP_cut*>& slack_pool,
		     BCP_vec<int>& to_be_purged){
	//std::cout<<"try to purge "<<slack_pool.size()<<std::endl;	 
	BCP_lp_user::purge_slack_pool(slack_pool, to_be_purged);	 
}

//#############################################################################
//#############################################################################
//#############################################################################

double
MCND_lp::compute_lower_bound(const double old_lower_bound,
                             const BCP_lp_result& lpres,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts){
    

    const int tc = lpres.termcode();
    //std::cout<<"compute lower bound: "<<std::endl;
    //if((lpres.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
    	//std::cout<<"compute_lower_bound aborted"<<std::endl;
		//abort();
	//}
	 
    //std::cout<<"compute lower bound: "<<(tc & BCP_ProvenOptimal)<<" "<<(tc & BCP_PrimalObjLimReached)<<std::endl;
    LBi = std::max(lpres.objval(),LBi);
    
    return LBi;
    
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::process_lp_result(const BCP_lp_result& lpres,
                           const BCP_vec<BCP_var*>& vars,
                           const BCP_vec<BCP_cut*>& cuts,
                           const double old_lower_bound,
                           double& true_lower_bound,
                           BCP_solution*& sol,
                           BCP_vec<BCP_cut*>& new_cuts,
                           BCP_vec<BCP_row*>& new_rows,
                           BCP_vec<BCP_var*>& new_vars,
                           BCP_vec<BCP_col*>& new_cols){
	
	getLpProblemPointer()->user_has_lp_result_processing = false;
 	if(lp_mode & LP_ForceNodeAbort) return;
	
	/*if(current_level()==0){
		for (int a=data.narcs; a--;)
			std::cout<<"sol y_"<<a<<" "<<lpres.x()[a]<<std::endl;
	}*/
	lp_mode |= LP_tighterBounds;
    const double *y_vol = lpres.x();
    y.assign(y_vol, y_vol+data.narcs);
    
    return;
}

//#############################################################################
//#############################################################################
//#############################################################################

BCP_solution*
MCND_lp::test_feasibility(const BCP_lp_result& lp_result,
                          const BCP_vec<BCP_var*>& vars,
                          const BCP_vec<BCP_cut*>& cuts){
     if((lp_mode & LP_StrongBranch || 
    	lp_mode & LP_ForceNodeAbort ||
    	lp_mode & LP_HeuristicRunned)) return 0;

	const double *sol = lp_result.x();
	double *fixd = new double [data.narcs];  
	double fathmval=0; 
	int cont=0;
	int closed=0;	
    bool integer=true;
	for (int a=data.narcs; a--;){
         if(vars[a]->lb()>0.5){fathmval += data.arcs[a].f; 
         }else if(vars[a]->ub()>0.5){ 
         	++cont;
         	if(sol[a]>1e-8 && sol[a]<(1-1e-8)) integer=false;
         	else if(sol[a]<1e-8){ fixd[closed++]=a;}
         }else{
          	fixd[closed++]=a;  
         }
    }
    //std::cout<<"violation: "<<getOsiVolBabSolver()->getViolation()<<" integ: "<<integer<<" unfx: "<<cont<<std::endl;

    int ret=-1;
    if(getOsiVolBabSolver()->getViolation()==0 && integer) {
    	ret = lpfeaschecker.solve_opt(vars, sol);
    	//std::cout<<"optimal::no_branch"<<std::endl;
    	lp_mode |= LP_ForceNodeAbort ;
    }else if(cont==0 || integer){ ret = lpfeaschecker.solve_opt(vars, 0);}
 
    if(ret<0){ 
     	if(ret<-1){
     		//std::cout<<"MCND_lp::test_feasibility::add_external_globalc type 0"<<std::endl;
     		getOsiVolBabSolver()->add_external_globalc( fixd, 0, closed);
     	}
     	delete []fixd;
    	return 0;
    } 
    delete []fixd;
    MCND_solution* mipsol = new MCND_solution(data.narcs+data.narcs*data.ndemands);
	ret = lpfeaschecker.getSolution(vars, mipsol->xy, mipsol->cost, fathmval);
	//std::cout<<"test feasibility unfx "<<cont<<" int "<<integer<<" : "<<mipsol->cost<<std::endl;
	if(ret<0 || fathmval >= upper_bound()){
		//std::cout<<"fathom value: "<<fathmval<<" ret: "<<ret<<" "<<cont<<std::endl;
		lp_mode |= LP_ForceNodeAbort ;
	}
	
	if(cont ==0){ lp_mode |= LP_ForceNodeAbort;}
 
	
	lp_mode |= LP_HeuristicRunned;
	    		
	if(upper_bound() > mipsol->cost){
		has_sol =true;
		//std::cout<<"MCND_lp::test_feasibility::add_external_globalc type 1"<<std::endl;
		getOsiVolBabSolver()->add_external_globalc(best_sol.xy, 1);
		best_sol.copy(*sol, data.narcs);
		lp_mode |= LP_tighterBounds;
		return mipsol;
	}else{
    	//std::cout<<"MCND_lp::test_feasibility::add_external_globalc type 1"<<std::endl;
		getOsiVolBabSolver()->add_external_globalc( mipsol->xy, 1);
	}
	delete mipsol;
	return 0;
	
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::restore_feasibility(const BCP_lp_result& lpres,
                             const std::vector<double*> dual_rays,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts,
                             BCP_vec<BCP_var*>& vars_to_add,
                             BCP_vec<BCP_col*>& cols_to_add){
    
    //std::cout<<"restore_feasibility"<<std::endl;
}

//==========================================================================================

BCP_solution*
MCND_lp::generate_heuristic_solution(const BCP_lp_result& lpres,
                                     const BCP_vec<BCP_var*>& vars,
                                     const BCP_vec<BCP_cut*>& cuts){
    //if(current_level() >100000){ lp_mode=LP_Normal;   return 0;}
    if(lp_mode & LP_ForceNodeAbort ) return 0;

    const double * x = lpres.x();
    double *fixd = new double [data.narcs];
    std::deque<Pair2> topo;
     
    int closed=0;
    int unfixed=0;
    int cont=0;
    for (int a=data.narcs; a--;){
        if(vars[a]->lb()==1.0){
            topo.push_back(Pair2(a, 1));
        }else if(vars[a]->ub()==1.0){
            if(x[a]>=0.9){ 
            	topo.push_back(Pair2(a, 1));
            }else if(x[a]>=0.3){
            	if(x[a]>0.5){
            		topo.push_front(Pair2(a, 1));
            	}else{
            		topo.push_front(Pair2(a, 0));
            	} 
            	++unfixed;
            }else{fixd[closed++]=a; }//std::cout<<"cand: "<<a<<" "<<x[a]<<std::endl;}
        	++cont;
        }else fixd[closed++]=a; 
    }
    
    if(cont>0){
    	if(lp_mode & LP_HeuristicRunned || num_nodes%10 || unfixed>=data.narcs*0.1){ 
    		topo.clear(); delete []fixd;  return 0;}
	}
    
	//std::cout<<"try heuristic "<<std::endl;

    MCND_solution* sol = new MCND_solution(data.narcs+data.narcs*data.ndemands);
    pump_heur.reset( unfixed, topo);
    int retval = pump_heur.solve(vars, sol->xy, sol->cost);
    topo.clear();
    lp_mode |= LP_HeuristicRunned;
    if(retval>=0){
    	delete []fixd;
    	if(upper_bound() > sol->cost){
    		has_sol =true;
    		//std::cout<<"MCND_lp::generate_heuristic_solution::add_external_globalc type 1"<<std::endl;
    		getOsiVolBabSolver()->add_external_globalc(best_sol.xy, 1);
    		best_sol.copy(*sol, data.narcs);
    		lp_mode |= LP_tighterBounds;
    		return sol;
    	}else{
    		//std::cout<<"MCND_lp::generate_heuristic_solution::add_external_globalc type 1"<<std::endl;
    		getOsiVolBabSolver()->add_external_globalc(sol->xy, 1);
    	}
    }else{
    	//std::cout<<"MCND_lp::generate_heuristic_solution::add_external_globalc type 0"<<std::endl;
		//getOsiVolBabSolver()->add_external_globalc(fixd, 0, closed );
		delete []fixd;
     }
    delete sol;
    return 0;
}


//#############################################################################
//#############################################################################
//#############################################################################
// helper functions

void
MCND_lp::pack_feasible_solution(BCP_buffer& buf, const BCP_solution* sol){
    //std::cout<<"pack feas sol: "<<std::endl;
    const MCND_solution* mcndsol = dynamic_cast<const MCND_solution*>(sol);
    
    mcndsol->pack(buf);
    
}

//=======================================================================================

void
MCND_lp::update_branch_data(MCND_node_branch_data * ndata){
    //if(ndata){
        //std::cout<<"update branch data"<<std::endl;
        //int index = ndata->pos_neg*data.narcs+ndata->branch_var;
        
        //std::cout<<index<<" new phi "<<branch_data.phi[index]<<" "<<ndata->branch_var;
        //delete ndata;
    //}
}

//=======================================================================================


void 
MCND_lp::keep_track(const MCND_CutUnit* vi){
	track.push_back(vi);
}


//#############################################################################


OsiVolSolverInterface*
MCND_lp::getOsiVolBabSolver(){
    return dynamic_cast<OsiVolSolverInterface*>(getLpProblemPointer()->lp_solver);
}

//#############################################################################


MCND_lp::~MCND_lp(){
	y.clear(); 
	ninsp.clear(); 
	psdcost.clear();
	mapd.clear();
    for(int i=track.size();i--;){
    	delete track[i];
    }
    track.clear();
     
    delete [] yfix;
}



