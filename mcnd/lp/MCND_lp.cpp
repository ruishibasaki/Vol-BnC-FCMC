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
    buf.unpack(no_heur);
    buf.unpack(has_sol);
    if(has_sol) best_sol.unpack(buf);
	
    y.resize(data.narcs,0);
    yfix = new int[data.narcs];
    std::fill(yfix,yfix+data.narcs, -1);
    ninsp.resize(data.narcs,Pair(0,0));
    psdcost.resize(data.narcs,PairF(0,0));
    tabu.resize(data.narcs,0);
    
    cover_manager.initialize(&data, 20);
    localc_manager.initialize(&data, 50);
    globalc_manager.initialize(&data, 50);
    ss_manager.initialize(&data);
    
    pump_heur.initialize(&data, &best_sol);
    lpfeaschecker.initialize(&data);
    flwconnect.initialize(&data);
    lpchecker.initialize(&data, &cover_manager.covers, &localc_manager.locals, &globalc_manager.globals);
    //topo_heur.initialize(&data, &best_sol, &cover_manager.covers, &localc_manager.locals, &globalc_manager.globals);
    
    AppVolData.data = &data;
    AppVolData.ss_manager = &ss_manager;
    AppVolData.cover_manager = &cover_manager;
    AppVolData.localc_manager = &localc_manager;
    AppVolData.globalc_manager = &globalc_manager;
    AppVolData.lpchecker = &lpchecker;
    
    nomgappres=0.001;
    max_cand = 5 ;
    maxszunfx=0.3*data.narcs;
	pump_heur.maxunfix = maxszunfx;
     
}

//-------------------------------------------------------------------------------------------

OsiSolverInterface *
MCND_lp::initialize_solver_interface(){
    OsiVolSolverInterface* volsolver = new OsiVolSolverInterface();
	volsolver->initialize( AppVolData);
    setOsiBabSolver(volsolver);
    //MaxIt = volsolver->volprob_.parm.maxsgriters;
    //std::cout<<"vol instantiated "<<volsolver<<std::endl;
    return volsolver;
    
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::process_message(BCP_buffer& buf){
	buf.unpack(lower_bound);
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
    
    if(has_sol){
    	double ub = upper_bound();
    	if(!(lp_mode & LP_Dive)){
			if((nomgap - (ub-lower_bound))/nomgap < nomgappres){
				no_gap_reduct += 1;
			}else no_gap_reduct=0;
			std::cout<<"no_gap_reduct: "<<no_gap_reduct<<" impv: "<<(nomgap - (ub-lower_bound))/nomgap<<std::endl;
			if(no_gap_reduct<=5){ 
				nomgap = (ub-lower_bound); 
				nomgappres = 0.001;
				reduced_run=false;
			}else{
				reduced_run=true;
				nomgappres = 0.005;
			}
			times_dived=0;
			//max_cand = 5 ;
		}else{ reduced_run=true;}//max_cand = 1;}
		double gap = (ub  - lower_bound)/ub;
		if(gap<0.005){
			maxszunfx=0;
			pump_heur.maxunfix = maxszunfx;
		} 
    	
    }
    lp_mode = LP_Normal;
    std::fill(yfix,yfix+data.narcs, -1);
    has_sol = getLpProblemPointer()->has_ub();
   
    bool testconn=false;
    MCND_node_branch_data* nodedata = dynamic_cast<MCND_node_branch_data*>( get_user_data());
    if(nodedata!=0){ 
        //std::cout<<"user_data hotstart: "<<nodedata->hs<<std::endl;
        LBi = nodedata->min_lb;
        lp_mode |= LP_TighterBounds;
        if(nodedata->hs!=0){
            getLpProblemPointer()->lp_solver->setWarmStart(nodedata->hs);
        }
    	std::cout<<"node comes from branch on: "<<nodedata->branch_var;
        std::cout<<" of "<<nodedata->pos_neg<<" side .. parent: "<<nodedata->parent<<std::endl;
        testconn = nodedata->test_conn ;
        //if(nodedata->reduced_run)lp_mode |= LP_ReducedRun;
        //std::cout<<"reduced run? "<<nodedata->reduced_run<<std::endl;

    }else LBi=0;
    int ret ;
    bool arcfx = false;
    if(testconn){
     	ret = flwconnect.check_connectivity(vars, var_changed_pos, var_new_bd, arcfx, yfix);
    	if(ret<0){
			if(ret==-1){ globalc_manager.globalc_generation_main2(0, yfix); }
			var_changed_pos.clear();
			var_new_bd.clear();
			std::cout<<"MCND_lp::initialize_new_search_tree_node::flwconnect.check_connectivity NOT"<<std::endl;
			lp_mode |= LP_ForceNodeAbort;
			return;
		}
		if(arcfx) lp_mode |= LP_LogicalFixed;	
    } 
    
    if(!cut_varfix_and_updt( vars, cuts, var_changed_pos, var_new_bd)) return;
    
    if(testconn || (lp_mode & LP_LogicalFixed)){
		flwconnect.reset(vars, yfix);
		ret = flwconnect.check_upperbounds(vars, var_changed_pos, var_new_bd, yfix);
		flwconnect.clear_memory();
		if(ret<0){
			var_changed_pos.clear();
			var_new_bd.clear();
			std::cout<<"MCND_lp::initialize_new_search_tree_node::flwconnect.check_upperbounds NOT"<<std::endl;
			lp_mode |= LP_ForceNodeAbort;
			return;
		}
	} 
    
    lp_mode |= LP_TestFeasibility;
    /*for (int a=var_changed_pos.size(); a--;){
    	if((var_new_bd[a*2]< vars[var_changed_pos[a]]->lb()) || (var_new_bd[a*2+1] > vars[var_changed_pos[a]]->ub()) ){
    		std::cout<<"problema: "<<var_changed_pos[a]<<" "<<var_new_bd[a*2]<<" "<<vars[var_changed_pos[a]]->lb();
    		std::cout<<" "<<var_new_bd[a*2+1]<<" "<<vars[var_changed_pos[a]]->ub()<<std::endl;

    	}
    }*/

     
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
                      BCP_var_set& vars, BCP_cut_set& cuts){
	//if(lp_mode & LP_ForceNodeAbort) return;
    //std::cout<<"load problem "<<core->varnum()<<" "<<cuts.size()<<" "<<getOsiVolBabSolver()<<std::endl;
    cover_manager.clean_collection();
    localc_manager.clean_collection();
	globalc_manager.clean_collection();
	
    int cutnum = cuts.size();
  	int core_rownum = core->cutnum();
  	
  	if (cutnum > core_rownum){
  		for(int i=core_rownum; i<cutnum;++i){
  			MCND_Cut * mcnd_cut = dynamic_cast<MCND_Cut*>(cuts[i]);
  			switch(mcnd_cut->cut_type){
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
    VOL_parms& par = AppVolData.volprob.parm;
    if(lp_mode & LP_ForceNodeAbort){ vollp->mode=-2;  return;}

    vollp->has_sol = has_sol;
    vollp->recheck_collct=false;
    vollp->upper_bound = LBi >0 ? fmin(upper_bound(), LBi*5) : upper_bound();
    par.dual_limit = vollp->upper_bound*5; 
    //if(no_gap_reduct >3 && !lpchecker.solved && !in_strong_branching) vollp->exact_solve=true;
    //else 
    vollp->exact_solve=false;
    //std::cout<<"limit: "<<par.dual_limit<<" "<<vollp->upper_bound<<" strongbranching:"<<in_strong_branching<<std::endl;     
    
    if(current_level()>0){ 
		par.lambdainit=1.0;
		par.alphainit=0.5;
		par.alphamin=1e-2;
		par.alphafactor=0.3;
		par.alphaint=10;
		AppVolData.minIterVI = 50;
     	AppVolData.intvlVI = 200; 
     	if(reduced_run){
     		//par.ascent_first_check = 50;
    		par.ascent_check_invl = 50;
     	}else{
     		//par.ascent_first_check = 100;
    		par.ascent_check_invl = 100;
     	}
     	vollp->mode=2;
    }else vollp->mode=1;
   
    if(in_strong_branching){
    	vollp->recheck_collct=true; 
    	vollp->min_lower_bound = LBi;  
    	vollp->mode=0;
    	return;
	}else if(changeType==1){ 
     	if(!(lp_mode & LP_LogicalFixed)) vollp->mode=-1;  
    }else if((lp_mode & LP_CutAddedFromHeuristic) && !(lp_mode & LP_LogicalFixed)){
    	vollp->mode=-1; 
    }else if(current_iteration()>1 && changeType==0) vollp->mode=-1;
    
    if(lp_mode & LP_LogicalFixed){
		vollp->recheck_collct=true;  
		if(lp_mode & LP_SecondIter)vollp->mode=3;
		lp_mode &= ~LP_LogicalFixed; /*vollp->mode=3; par.maxsgriters = 100;*/
	}
    vollp->in_strong_branch = in_strong_branching;
    lp_mode |= LP_SecondIter;
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
    LBi = std::max(lpres.objval(),LBi);
    if(lp_mode & LP_ForceNodeAbort) return LBi;

    double ub = upper_bound();
	double gap = (ub  - LBi)/ub;
    if(((gap < 0.01) && !(lp_mode & LP_OptSolved) && !reduced_run)|| times_dived>=3){
    	int ret  = lpchecker.solve(ub, vars, 0,0, LBi);
    	if(ret<0){
    		lp_mode |= LP_ForceNodeAbort;
    		std::cout<<"MCND_lp::compute_lower_bound:: compute opt lower bound: INFEASIBLE"<<std::endl;
    	}else{
     		lpchecker.solved=true;
     		lp_mode |= LP_OptSolved;
    	} 
    }else lpchecker.solved=false;
   
   /* if(force_dive){
    	std::cout<<"MCND_lp::compute_lower_bound break dive? "<<(nomgap - (upper_bound()-LBi))/nomgap <<std::endl;
    	if((nomgap - (upper_bound()-LBi))/nomgap < 0.01)
			break_diving=false;
		else{
			break_diving=true;
			//std::cout<<"MCND_lp::compute_lower_bound break dive"<<std::endl;
		}
	}*/

	 
    //std::cout<<"compute lower bound dev: "<<(LBi - lower_bound)/lower_bound<<std::endl;
    
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
	if((lpres.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
    	lp_mode |= LP_ForceNodeAbort;
    	std::cout<<"MCND_lp::process_lp_result:: BCP_ProvenPrimalInf"<<std::endl;
    	return;
	}else if((lpres.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached){
    	lp_mode |= LP_ForceNodeAbort;
    	std::cout<<"MCND_lp::process_lp_result:: BCP_PrimalObjLimReached"<<std::endl;
    	return;
	}else if(lp_mode & LP_ForceNodeAbort) return;
	
	std::cout<<"process_lp_result:: "<<lpres.objval()<<std::endl;
	if(lpres.objval() > LBi+1e-4) lp_mode |= LP_TighterBounds;
	/*if(current_level()==0){
		for (int a=data.narcs; a--;)
			std::cout<<"sol y_"<<a<<" "<<lpres.x()[a]<<std::endl;
	}*/
	//if(!(lp_mode & LP_Solved)) return;
	
	//std::cout<<"process_lp_result"<<std::endl;
	//lp_mode |= LP_TighterBounds;
    const double *y_vol = lpres.x();
    double val;
    for (int a=data.narcs; a--;) {
    	val = y_vol[a];
    	y[a] = val;
		//if(vars[a]->lb()==0 && vars[a]->ub()==1){
		//	if(val>1e-8) --tabu[a];
		//	else if(val<=1e-8) ++tabu[a];
		//}
	}
    return;
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
	if(lp_mode & LP_CutAddedFromHeuristic) lp_mode &= ~LP_CutAddedFromHeuristic;
	if((lp_mode & LP_StrongBranch) ) return;
    
    int sz;
 	Cover *vi = cover_manager.covers.begin;
 	LocalCut * lcvi = localc_manager.locals.begin;
    if(lp_mode & LP_ForceNodeAbort){
    	if(cover_manager.gend>0){
    		cover_manager.gend =0;
    		sz = cover_manager.covers.sizeOfCollection;
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
    
    if(	lp_mode & LP_TighterBounds){
    	//std::cout<<"generate_cuts_in_lp locals type 1"<<std::endl;
        getOsiVolBabSolver()->add_external_localc(yfix, 0, data.narcs, 1);
		lp_mode &= ~LP_TighterBounds;
    }
    //std::cout<<"generate_cuts_in_lp "<<localc_manager.gend<<std::endl;
    if(cover_manager.gend==0 && localc_manager.gend==0) return;
	
	if(cover_manager.gend>0){
		vi = cover_manager.covers.begin;
		sz = cover_manager.covers.sizeOfCollection;
		for(int i=sz;i--;){
			if(vi->toadd){ 
				keep_track(vi);
				new_cuts.push_back(new CoverCut(vi)); 
				vi->toadd=false;
				++newly_added;
			}
			vi = vi->next;
		}
		cover_manager.gend=0;
	}
	
	if(localc_manager.gend>0){
		lcvi = localc_manager.locals.begin;
		sz = localc_manager.locals.sizeOfCollection;
		for(int i=sz;i--;){
			if(lcvi->toadd){ 
				keep_track(lcvi);
				//std::cout<<"generate lc "<<lcvi->serial_nmbr<<" "<<lcvi->type<<std::endl;
				//lcvi->print();
				new_cuts.push_back(new LocalCCut(lcvi)); 
				lcvi->toadd=false;
				++newly_added;
			}
			lcvi = lcvi->next;
		}
		localc_manager.gend=0;
	}
    
    //std::cout<<"generate cuts: "<<newly_added<<" cutsvecsz: "<<new_cuts.size()<<" "<<sz<<std::endl;
 
    
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
			switch(cut->cut_type){
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

BCP_solution*
MCND_lp::test_feasibility(const BCP_lp_result& lp_result,
                          const BCP_vec<BCP_var*>& vars,
                          const BCP_vec<BCP_cut*>& cuts){
     if(lp_mode & LP_StrongBranch || 
    	lp_mode & LP_ForceNodeAbort) return 0;
    
    MCND_solution* mipsol;
    if(!(lp_mode & LP_TestFeasibility)) return 0;
 	lp_mode &= ~LP_TestFeasibility;

	const double *sol = lp_result.x();
	int *fixd = new int [data.narcs];  
	double fathmval=0; 
	int cont=0;
	int closed=0;	
    bool integer=true;
	for (int a=data.narcs; a--;){
		if(vars[a]->lb()>0.5){fathmval += data.arcs[a].f; 
		}else if(vars[a]->ub()>0.5){ 
			++cont;
			if(sol[a]>0.1 && sol[a]<(1-0.1)) integer=false;
		} 
    }
	std::cout<<"test_feasibility unfx: "<<cont<<std::endl;
    int ret=10;
    if(getOsiVolBabSolver()->getViolation()==0 && integer) {
    	ret = lpfeaschecker.solve_opt(cont, vars, sol, fixd, closed, mipsol, fathmval, upper_bound());
    	std::cout<<"MCND_lp::test_feasibility intger && no violation"<<std::endl;
    	lp_mode |= LP_ForceNodeAbort ;
    }else if(cont<= maxszunfx || integer ){ 
    	ret = lpfeaschecker.solve_opt(cont, vars, 0, fixd, closed, mipsol, fathmval, upper_bound());
    	if(cont ==0){ lp_mode |= LP_ForceNodeAbort; std::cout<<"MCND_lp::test_feasibility allfixed"<<std::endl;}
    }else {
    	delete []fixd;
    	return 0;
    }
 
 	if(ret<0){
		if(ret==-1){
			std::cout<<"MCND_lp::test_feasibility::add_external_globalc type 0"<<std::endl;
			getOsiVolBabSolver()->add_external_globalc( fixd, closed);
		}
		std::cout<<"MCND_lp::test_feasibility NOT"<<std::endl;
		lp_mode |= LP_ForceNodeAbort;
		delete []fixd;
    	return 0;
	}else if(ret==0 || ret==2){lp_mode |= LP_ForceNodeAbort;
	}else if(!(lp_mode & LP_ForceNodeAbort)){
		std::cout<<"test_feasibility integer? "<<integer<<std::endl;
		if(lpfeaschecker.test_high_lowerbounds(vars, 0, sol , upper_bound())<0){
			lp_mode |= LP_ForceNodeAbort;
			return 0;
		}
			
	}
	delete []fixd;
	
	 
    //std::cout<<"test_feasibility: "<<getOsiVolBabSolver()->getViolation()<<" integ: "<<integer<<" unfx: "<<cont<<std::endl;
    std::cout<<"test_feasibility: sol: "<<mipsol->cost<<" "<<cont<<" ub: "<<upper_bound()<<std::endl;
    lp_mode |= LP_HeuristicRunned;
	if(upper_bound() >= mipsol->cost){
		if(!(lp_mode & LP_ForceNodeAbort)){
			std::cout<<"MCND_lp::test_feasibility::add_external_localc1 type 2"<<std::endl;
			getOsiVolBabSolver()->add_external_localc( 0, best_sol.xy, data.narcs, 2);
		}
		best_sol.copy(*mipsol, data.narcs);
		lp_mode |= LP_TighterBounds;
		return mipsol;
	}else if(!(lp_mode & LP_ForceNodeAbort)){
    	std::cout<<"MCND_lp::test_feasibility::add_external_localc type 2"<<std::endl;
    	getOsiVolBabSolver()->add_external_localc( 0, mipsol->xy, data.narcs, 2);
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
    //return 0;
    if(lp_mode & LP_ForceNodeAbort  ) return 0;
    
	MCND_solution* sol=0;
	if(lpchecker.solved && current_iteration()==1){
    	sol = new MCND_solution(data.narcs+data.narcs*data.ndemands);
    	lpchecker.get_solution(sol);
    	if(upper_bound() > sol->cost){
			std::cout<<"MCND_lp::lpchecker.get_solution::add_external_localc11 type 2"<<std::endl;
			getOsiVolBabSolver()->add_external_localc( 0, best_sol.xy, data.narcs, 2);
			lp_mode |= LP_CutAddedFromHeuristic;
			best_sol.copy(*sol, data.narcs);
			lp_mode |= LP_TighterBounds; 
			return sol;
		}else{
			std::cout<<"MCND_lp::lpchecker.get_solution::add_external_localc12 type 2"<<std::endl;
			getOsiVolBabSolver()->add_external_localc( 0, sol->xy, data.narcs, 2);
			lp_mode |= LP_CutAddedFromHeuristic;
		}
		delete sol;
		//return 0;
    }
     
    const double * x = lpres.x();
    int * fixd0 = 0;
    int closed=0;
	int cont= pump_heur.make_topo( x, vars);
    
	//if(cont!=0) return 0; //
 	if(cont > 0){
 		//std::cout<<"Pump::solve trypump? "<<cont<<" "<<maxszunfx<<std::endl;
 		if(no_heur) return 0; 
    	else if((lp_mode & LP_HeuristicRunned) && (current_level()>0 || current_iteration()>15)) return 0;
    	else if( current_level()%50>0 || (current_iteration()>15)) return 0;
    	
	}
    //if(pump_heur.validate_topology()< 0) return 0;
    
	//std::cout<<"try heuristic "<<cont<<std::endl;
    int retval = pump_heur.solve( fixd0, closed, vars, sol); 
    lp_mode |= LP_HeuristicRunned;
    if(retval>0){
     	std::cout<<"heuristic sol: "<<sol->cost<<std::endl;
     	if(retval==2){
     	   std::cout<<"MCND_lp::generate_heuristic_solution::add_external_globalc00 type 0"<<std::endl;
     		getOsiVolBabSolver()->add_external_globalc( fixd0, closed);
			if(fixd0) delete [] fixd0;
		}
    	if(upper_bound() >= sol->cost){
    		has_sol =true;
    		std::cout<<"MCND_lp::generate_heuristic_solution::add_external_localc1 type 2"<<std::endl;
    		if(getOsiVolBabSolver()->add_external_localc( 0, best_sol.xy, data.narcs, 2))
    			lp_mode |= LP_CutAddedFromHeuristic;
    		best_sol.copy(*sol, data.narcs);
    		lp_mode |= LP_TighterBounds;
    		return sol;
    	}else{
    		std::cout<<"MCND_lp::generate_heuristic_solution::add_external_lobalc2 type 2"<<std::endl;
    		if(getOsiVolBabSolver()->add_external_localc( 0, sol->xy, data.narcs, 2))
    			lp_mode |= LP_CutAddedFromHeuristic;
    		delete sol;
    	}
    }else if(retval==0 || retval==-1){
     		std::cout<<"MCND_lp::generate_heuristic_solution::add_external_globalc1 type 0"<<std::endl;
			getOsiVolBabSolver()->add_external_globalc( fixd0, closed);
			if(fixd0) delete [] fixd0;
    }else if(retval>-10){
    	std::cout<<"MCND_lp::generate_heuristic_solution::add_external_localc3 type 0"<<std::endl;
		if(getOsiVolBabSolver()->add_external_localc(fixd0, 0, closed, 0))
			lp_mode |= LP_CutAddedFromHeuristic;
		if(fixd0) delete [] fixd0;
    }
   
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
MCND_lp::keep_track(MCND_CutUnit* vi){
	track.push_back(vi);
}


//#############################################################################


OsiVolSolverInterface*
MCND_lp::getOsiVolBabSolver(){
    return dynamic_cast<OsiVolSolverInterface*>(getLpProblemPointer()->lp_solver);
}

//#############################################################################


MCND_lp::~MCND_lp(){
	localc_manager.clean_collection();
	cover_manager.clean_collection();
	y.clear(); 
	ninsp.clear(); 
	psdcost.clear();
	mapd.clear();
   for(int i=track.size();i--;){
    	delete track[i];
    }
    track.clear();
    to_logical_fix.clear();
    tabu.clear();
    
    delete [] yfix;
}



