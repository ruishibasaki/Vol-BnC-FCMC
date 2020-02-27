// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#include <cmath>
#include <algorithm>
#include <numeric>

#include "BCP_math.hpp"
#include "BCP_lp.hpp"
#include "BCP_lp_node.hpp"

#include "MCND_lp.hpp"
#include "MCND_mincostflow.hpp"

//#############################################################################

void
MCND_lp::unpack_module_data(BCP_buffer & buf){
    //std::cout<<"try unpack to lp "<<std::endl;
    data.unpack(buf);
    y.resize(data.narcs,0);
    
    ss_manager.initialize(&data);
    cover_manager.initialize(&data, 10);
    pump_heur.set_data(&data, 0.5);
    
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
    MaxIt = volsolver->volprob_.parm.maxsgriters;
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
    MCND_node_branch_data* nodedata = dynamic_cast<MCND_node_branch_data*>( get_user_data());
    if(nodedata!=0){
        //std::cout<<"user_data hotstart: "<<nodedata->hs<<std::endl;
        if(nodedata->hs!=0){
            getLpProblemPointer()->lp_solver->setWarmStart(nodedata->hs);
        }
    }
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
                      BCP_var_set& vars, BCP_cut_set& cuts){
    
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
    std::cout<<"modify lp param "<<changeType<<" "<<in_strong_branching<<" iter: "<<current_level()<<std::endl;
    lp->setDblParam(OsiPrimalTolerance, 1e-4);
    OsiVolSolverInterface* vollp = getOsiVolBabSolver();
    VOL_parms& par = vollp->volprob_.parm;
    
    vollp->has_sol = has_sol;
    vollp->upper_bound = LBi >0 ? std::min(upper_bound(), LBi*5) : upper_bound();
    par.dual_limit = vollp->upper_bound*5; std::cout<<"limit: "<<par.dual_limit<<" "<<vollp->upper_bound<<std::endl;

    
    if(!vollp->HotStartSet){
    	if(current_level()==0) par.maxsgriters = 1000;
    	else par.maxsgriters = 500;
    }else par.maxsgriters = 250;
     
    if(current_level()>0){  AppVolData.intvlVI = 100;}
    if(in_strong_branching){ vollp->mode=0;}
    else if(changeType==1 && !in_strong_branching){ vollp->mode=-1;}
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

double
MCND_lp::compute_lower_bound(const double old_lower_bound,
                             const BCP_lp_result& lpres,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts){
    
    const int tc = lpres.termcode();
    std::cout<<"compute lower bound: "<<(tc & BCP_ProvenOptimal)<<" "<<(tc & BCP_PrimalObjLimReached)<<std::endl;

    if (tc & BCP_ProvenOptimal)
        return lpres.objval();
    
    
    if (tc & BCP_PrimalObjLimReached)
        return best_soln.cost + 1;
    
    
    return old_lower_bound;
    
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

    const double *y_vol = lpres.x();
    y.assign(y_vol, y_vol+data.narcs);
    getLpProblemPointer()->user_has_lp_result_processing = false;
    return;
    /*std::cout<<"process_lp_result and generate"<<std::endl;
    getLpProblemPointer()->user_has_lp_result_processing = false;
    const double *y_vol = lpres.x();
    y.assign(y_vol, y_vol+data.narcs);
    sol = test_feasibility(lpres, vars, cuts);
    
    
    
    if(!((lpres.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf))
        true_lower_bound = lpres.objval();*/
    
    
    
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
    std::cout<<"generate cuts: "<<newly_added<<" sz: "<<sz<<" cutsvecsz: "<<new_cuts.size()<<std::endl;
    for(int i=newly_added;i--;){
    	keep_track(vi);
        new_cuts.push_back(new CoverCut(vi)); 
        vi = vi->prev;
    }
    if(newly_added){
    	lp_mode = LP_CutAdded;
    	cover_manager.gend =0;
    } 
    
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
    BCP_lp_user::select_vars_to_delete(lpres,vars,cuts,before_fathom, deletable);
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::select_cuts_to_delete(const BCP_lp_result& lpres,
				   const BCP_vec<BCP_var*>& vars,
				   const BCP_vec<BCP_cut*>& cuts,
				   const bool before_fathom,
				   BCP_vec<int>& deletable){
	std::cout<<"select_cuts_to_delete "<<before_fathom<<std::endl;
	
	//if(!before_fathom){
	const int cutnum = cuts.size();
	deletable.reserve(cutnum-getLpProblemPointer()->core->cutnum()+1);

	for (int i = getLpProblemPointer()->core->cutnum(); i < cutnum; ++i) {
		CoverCut * cut = dynamic_cast<CoverCut*>(cuts[i]);
		if (!cut->check_viol(vars)) {
			std::cout<<"out: "<<i<<" id: "<<cut->get_cover()->id_vi<<" srnb: "<<cut->get_cover()->serial_nmbr<<std::endl;
			cover_manager.purgbl.push_back(cut->get_cover());
			deletable.unchecked_push_back(i);
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

BCP_solution*
MCND_lp::test_feasibility(const BCP_lp_result& lp_result,
                          const BCP_vec<BCP_var*>& vars,
                          const BCP_vec<BCP_cut*>& cuts){
    
    std::cout<<"test feasibility "<<std::endl;
    MCND_node_branch_data * ndata  = dynamic_cast<MCND_node_branch_data *>(get_user_data());
    if(ndata){
        std::cout<<"node comes from branch on: "<<ndata->branch_var;
        std::cout<<" of "<<ndata->pos_neg<<" side "<<std::endl;
    }
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
    std::cout<<"try heuristic "<<candidates.size()<<std::endl;
    if(current_level() >100000){ lp_mode=LP_Normal;   return 0;}
	if(lp_mode == LP_CutAdded){ lp_mode = LP_Normal; return 0;}
	
    const double * x = lpres.x();
    double * yl = new double [data.narcs];
    std::deque<Pair2> topo;
    int unfixed=0;
    for (int a=data.narcs; a--;){
    	yl[a] = 0;
        if(vars[a]->lb()==1.0){
            topo.push_back(Pair2(a, 1));
        }else if(vars[a]->ub()==1.0){
            if(x[a]>=0.7) topo.push_back(Pair2(a, 1));
            else if(x[a]>=0.3){
            	if(x[a]>0.5){
            		topo.push_front(Pair2(a, 1));
            	}else topo.push_front(Pair2(a, 0));
            	++unfixed;
            }else{ candidates.push_back(Pair2(a,x[a])); }//std::cout<<"cand: "<<a<<" "<<x[a]<<std::endl;}
        }
    }
    MCND_solution* sol = new MCND_solution(data.narcs+data.narcs*data.ndemands);
    pump_heur.reset( unfixed, topo);
    int retval = pump_heur.solve(yl, sol->xy, sol->cost);
    
    topo.clear();
    if(retval>=0){
    	has_sol =true;
        lp_mode=LP_Normal;
    	candidates.clear();
    	if(upper_bound() > sol->cost)
    		return sol;
    	delete sol;
    	return 0;
    }
    lp_mode = LP_DiveToFeasibility;
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
    if(ndata){
        std::cout<<"update branch data"<<std::endl;
        //int index = ndata->pos_neg*data.narcs+ndata->branch_var;
        
        //std::cout<<index<<" new phi "<<branch_data.phi[index]<<" "<<ndata->branch_var;
        //delete ndata;
    }
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
    std::deque<const Cover *> track;
    candidates.clear();
    mapd.clear();
    for(int i=track.size();i--;){
    	delete track[i];
    }
    track.clear();
}



