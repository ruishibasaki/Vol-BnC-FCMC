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
    ninsp.resize(data.narcs,Pair(0,0));
    psdcost.resize(data.narcs,PairF(0,0));
    
    ss_manager.initialize(&data);
    cover_manager.initialize(&data, 10);
    pump_heur.set_data(&data, 0.5);
    lpfeaschecker.initialize(&data);
    flwconnect.initialize(&data);
    
    AppVolData.data = &data;
    AppVolData.ss_manager = &ss_manager;
    AppVolData.cover_manager = &cover_manager;
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
    
    
    ++num_nodes;
    lp_mode = LP_Normal;
    std::cout<<"initialize_new_search_tree_node "<<num_nodes<<" "<<lp_mode <<std::endl;
    
    has_sol = getLpProblemPointer()->has_ub();
    MCND_node_branch_data* nodedata = dynamic_cast<MCND_node_branch_data*>( get_user_data());
    if(nodedata!=0){
        //std::cout<<"user_data hotstart: "<<nodedata->hs<<std::endl;
        LBi = nodedata->min_lb;
        if(nodedata->hs!=0){
            getLpProblemPointer()->lp_solver->setWarmStart(nodedata->hs);
        }
        std::cout<<"node comes from branch on: "<<nodedata->branch_var;
        std::cout<<" of "<<nodedata->pos_neg<<" side "<<std::endl;
        if(nodedata->test_conn){
        	if(!flwconnect.check_connectivity(vars, var_changed_pos, var_new_bd)){
        		var_changed_pos.clear();
        		var_new_bd.clear();
        		std::cout<<"fathom node flwconnect.check_connectivity"<<std::endl;
        		lp_mode |= LP_ForceNodeAbort;
        	}	
        }
        if(nodedata->reduced_run)lp_mode |= LP_ReducedRun;
        std::cout<<"reduced run? "<<nodedata->reduced_run<<std::endl;

    }else LBi=0;
     
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
                      BCP_var_set& vars, BCP_cut_set& cuts){
	//if(lp_mode & LP_ForceNodeAbort) return;
    std::cout<<"load problem "<<core->varnum()<<" "<<cuts.size()<<std::endl;
    cover_manager.clean_collection();
    int cutnum = cuts.size();
  	int core_rownum = core->cutnum();
  	
  	if (cutnum > core_rownum){
  		for(int i=core_rownum; i<cutnum;++i){
  			CoverCut * cut = dynamic_cast<CoverCut*>(cuts[i]);
  			if(cut){
  				//std::cout<<"collec: "<<i<<" "<<cut->get_cover()->serial_nmbr<<" "<<cuts[i]->bcpind()<<std::endl;
  				cut->get_cover()->id_vi = i;
  				cover_manager.covers.insert_end(cut->get_cover());
  			}
  			else std::cout<<"load_problem::problem"<<std::endl;
  		}
  	}
	
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
    vollp->upper_bound = LBi >0 ? fmin(upper_bound(), LBi*5) : upper_bound();
    par.dual_limit = vollp->upper_bound*5; 
    std::cout<<"limit: "<<par.dual_limit<<" "<<vollp->upper_bound<<" strongbranching:"<<in_strong_branching<<std::endl;

    if(lp_mode & LP_ReducedRun){
    	 par.maxsgriters = 100;
    	 lp_mode &= ~LP_ReducedRun;
    }else if(!vollp->HotStartSet){
    	if(current_level()==0) par.maxsgriters = 1000;
    	else par.maxsgriters = 500;
    }else par.maxsgriters = 250;
     
    if(current_level()>0){  AppVolData.intvlVI = 100;}
    
    if(in_strong_branching){ vollp->min_lower_bound = LBi;  vollp->mode=0;}
    else if(changeType==1 || lp_mode & LP_ForceNodeAbort ){ vollp->mode=-1; std::cout<<"changeType "<<changeType<<" lpmode: "<<(lp_mode & LP_ForceNodeAbort)<<std::endl; }
    else if((lp_mode & LP_CutAddedFromHeuristic) ){vollp->mode=-1; std::cout<<"LP_CutAddedFromHeuristic"<<std::endl;}
    else{ vollp->mode=1;}
    vollp->in_strong_branch = in_strong_branching;

    //if(in_strong_branching ){
    //  std::cout<<"setAuxiliaryInfo "<<lp->getAuxiliaryInfo()<<std::endl;
    //lp->setAuxiliaryInfo(new MCND_parent_branch_data());
    //}
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
    

    int sz = cover_manager.covers.sizeOfCollection;
    int newly_added = cover_manager.gend;
	Cover *vi = cover_manager.covers.begin;
    cover_manager.covers.advance(vi, cover_manager.num_actv-1);
    if(lp_mode & LP_ForceNodeAbort){
     	for(int i=newly_added;i--;){
			keep_track(vi);
			vi = vi->prev;
		}
     	return;
    }

    //std::cout<<"new added "<<newly_added<<" "<<sz<<" "<<cover_manager.num_actv<<std::endl;
    for(int i=newly_added;i--;){
    	keep_track(vi);
    	//std::cout<<"generate "<<vi->id_vi<<" "<<vi->serial_nmbr<<std::endl;
        new_cuts.push_back(new CoverCut(vi)); 
        vi = vi->prev;
    }

    if(newly_added){
        std::cout<<"generate cuts: "<<newly_added<<" sz: "<<sz<<" cutsvecsz: "<<new_cuts.size()<<std::endl;
	    cover_manager.gend =0;
    } 
    
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
        row->LowerBound(1);
        row->UpperBound(1e40);
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
		CoverCut * cut = dynamic_cast<CoverCut*>(cuts[i]);
		if (!cut->check_viol(vars) || cut->get_cover()->prgbl ||
			(cut->check_viol(lpres.x())==0 && lpres.pi()[cut->get_cover()->id_vi]==0)) {
			//std::cout<<"out / dual: "<<lpres.pi()[cut->get_cover()->id_vi]<<" id: "<<cut->get_cover()->id_vi<<" srnb: "<<cut->get_cover()->serial_nmbr<<std::endl;
			cover_manager.purgbl.push_back(cut->get_cover());
			deletable.unchecked_push_back(i);
			cut->get_cover()->prgbl=false;
		}
	}
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
    if((lpres.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
    	std::cout<<"compute_lower_bound aborted"<<std::endl;
		abort();
	}
	 
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

	if(lp_mode & LP_ForceNodeAbort) return;

    const double *y_vol = lpres.x();
    y.assign(y_vol, y_vol+data.narcs);
    getLpProblemPointer()->user_has_lp_result_processing = false;
    
    return;
}

//#############################################################################
//#############################################################################
//#############################################################################

BCP_solution*
MCND_lp::test_feasibility(const BCP_lp_result& lp_result,
                          const BCP_vec<BCP_var*>& vars,
                          const BCP_vec<BCP_cut*>& cuts){
    
    if(lp_mode & LP_StrongBranch || 
    	lp_mode & LP_ForceNodeAbort ||
    	lp_mode & LP_HeuristicRunned ) return 0;

    std::cout<<"test feasibility "<<std::endl;
	const double *sol = lp_result.x();
	double fathmval=0; 
	int cont=0;
	bool integer=true;
	for (int a=data.narcs; a--;){
         if(sol[a]>1e-8 && sol[a]<(1-1e-8)) integer=false;
         if(vars[a]->lb()>0.5)fathmval += data.arcs[a].f;
         else if(vars[a]->ub()>0.5) ++cont;
    }
    std::cout<<"violation: "<<getOsiVolBabSolver()->getViolation()<<" integ: "<<integer<<" unfx: "<<cont<<std::endl;

    int ret=-1;
    if(getOsiVolBabSolver()->getViolation()==0 && integer) {
    	ret = lpfeaschecker.solve_opt(vars, sol);
    	std::cout<<"optimal::no_branch"<<std::endl;
    	lp_mode |= LP_ForceNodeAbort ;
    }else if(cont==0 || integer){ ret = lpfeaschecker.solve_opt(vars, 0);}

    if(ret<0){ return 0;} 
    
    MCND_solution* mipsol = new MCND_solution(data.narcs+data.narcs*data.ndemands);
	ret = lpfeaschecker.getSolution(vars, mipsol->xy, mipsol->cost, fathmval);
	if(ret<0 || fathmval >= upper_bound()){
		std::cout<<"fathom value: "<<fathmval<<" ret: "<<ret<<" "<<cont<<std::endl;
		lp_mode |= LP_ForceNodeAbort ;
	}
	
	if(cont ==0) lp_mode |= LP_ForceNodeAbort;
	lp_mode |= LP_HeuristicRunned;
	    		
	if(upper_bound() > mipsol->cost){
		has_sol =true;
		best_sol.copy(*sol, data.narcs);
		return mipsol;
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
    double * yl = new double [data.narcs];
    std::deque<Pair2> topo;
    std::deque<Pair2> forced0;
    
    int unfixed=0;
    int cont=0;
    for (int a=data.narcs; a--;){
    	yl[a] = 0;
        if(vars[a]->lb()==1.0){
            topo.push_back(Pair2(a, 1));
        }else if(vars[a]->ub()==1.0){
            if(x[a]>=0.9){ 
            	topo.push_back(Pair2(a, 1));
            }else if(x[a]>=0.1){
            	if(x[a]>0.5){
            		topo.push_front(Pair2(a, 1));
            	}else topo.push_front(Pair2(a, 0));
            	++unfixed;
            }else{ forced0.push_back(Pair2(a,1.0)); }//std::cout<<"cand: "<<a<<" "<<x[a]<<std::endl;}
        	++cont;
        }
    }
    
    if(cont>0){
    if(lp_mode & LP_HeuristicRunned){ topo.clear(); forced0.clear(); return 0;}
	else if(num_nodes%100){ topo.clear(); forced0.clear();  return 0;}
    else if(unfixed>=data.narcs*0.1){ topo.clear(); forced0.clear(); return 0; }}
    
	std::cout<<"try heuristic "<<std::endl;

    MCND_solution* sol = new MCND_solution(data.narcs+data.narcs*data.ndemands);
    pump_heur.reset( unfixed, topo);
    int retval = pump_heur.solve(vars, yl, sol->xy, sol->cost);
    delete [] yl;
    topo.clear();
    lp_mode |= LP_HeuristicRunned;
    if(retval>=0){
    	if(upper_bound() > sol->cost){
    		has_sol =true;
    		best_sol.copy(*sol, data.narcs);
    		return sol;
    	}
    		
    }else{
    	OsiVolSolverInterface* vollp = getOsiVolBabSolver();
    	vollp->add_external_vi(forced0);
    	lp_mode |= LP_CutAddedFromHeuristic;
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
MCND_lp::keep_track(const Cover* vi){
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
    std::deque<const Cover *> track;
    mapd.clear();
    for(int i=track.size();i--;){
    	delete track[i];
    }
    track.clear();
}



