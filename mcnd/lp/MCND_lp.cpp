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
// helper functions

MCND_solution * 
MCND_lp::SolveMinCostFLow(const std::deque<int> & topo){
    
    MinCostFLow lp;
    lp.create_model(topo, &data);
    int retval = lp.solve();
    if(retval < 0){
        std::cout<<"no good"<<std::endl;
        return 0;
    }else if(retval == 0){
        MCND_solution* solution = new MCND_solution(data.narcs+data.narcs*data.ndemands);
        lp.getSolution(topo, &data, solution);
        if(solution->cost<best_soln.cost){
            best_soln = *solution;
        }
        return solution;
    }
    std::cout<<"no good"<<std::endl;
    return 0;
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::update_branch_data(MCND_node_branch_data * ndata){
    if(ndata){
        std::cout<<"update branch data"<<std::endl;
        //int index = ndata->pos_neg*data.narcs+ndata->branch_var;
        
        //std::cout<<index<<" new phi "<<branch_data.phi[index]<<" "<<ndata->branch_var;
        //delete ndata;
    }
}

//-------------------------------------------------------------------------------------------

bool 
MCND_lp::process_flow_solution(const BCP_vec<BCP_var*>& vars, const MCND_solution * sol){
    bool cutoff = true;
    double f=0;
    for(int a=data.narcs;a--;){
        if(vars[a]->lb()==0 && vars[a]->ub()==1 && sol->xy[a]>0){
            cutoff= false;
            //std::cout<<"nop in a0: "<<a<<" "<<sol->xy[a]<<std::endl;
            break;
        }else if(vars[a]->lb()==1){
            f+=data.arcs[a].f;
        }
    }
    if(sol->cost_flow + f >= best_soln.cost){
        //std::cout<<"yes cost "<<sol->cost_flow + f<<" >= "<<best_soln.cost<<std::endl;
        cutoff = true;
    }//else std::cout<<"nop cost "<<sol->cost_flow + f<<" < "<<best_soln.cost<<std::endl;
    return cutoff;
}


//#############################################################################

void
MCND_lp::unpack_module_data(BCP_buffer & buf){
    std::cout<<"try unpack to lp "<<std::endl;
    data.unpack(buf);
    y.resize(data.narcs,0);
    freq.resize(data.narcs,0);
    freq_cand.resize(data.narcs,0);
    
    ss_manager.initialize(&data);
    cover_manager.initialize(&data, 10);
    
    AppVolData.data = &data;
    AppVolData.ss_manager = &ss_manager;
    AppVolData.cover_manager = &cover_manager;
    srand(10);
}

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
    
    
    std::cout<<"initialize_new_search_tree_node "<<cuts.size()<<std::endl;
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
MCND_lp::initialize_int_and_sos_list(std::vector<OsiObject *>& intAndSosObjects){
    std::cout<<"SOS"<<std::endl;
}

//-------------------------------------------------------------------------------------------

OsiSolverInterface *
MCND_lp::initialize_solver_interface(){
    OsiVolSolverInterface* volsolver = new OsiVolSolverInterface("volmcnd.par");
    setOsiBabSolver(volsolver);
    MaxIt = volsolver->volprob_.parm.maxsgriters;
    std::cout<<"vol instantiated "<<MaxIt<<std::endl;
    return volsolver;
    
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
  				if(cut->check_viol(vars)){ cover_manager.covers.insert_front(cut->get_cover());}
  				else{cuts[i]->make_to_be_removed();}
  				//else{ std::cout<<"out id_vi: "<<cut->get_cover()->id_vi<<" , "<<cut->get_cover()->serial_nmbr<<std::endl; cut->get_cover()->print();}
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
	
    par.ubinit = best_soln.cost;
    if(!vollp->HotStartSet){
    	if(current_level()==0)par.maxsgriters = 1000;
    	else par.maxsgriters = 500;
    }else par.maxsgriters = 250;
     
    if(in_strong_branching){ vollp->mode=0;  if(current_level()>0){ AppVolData.intvlVI = 100; vollp->mode=0; }}
    else if(changeType==1 && !in_strong_branching){ vollp->mode=-1;}
    else{ vollp->mode=1;  if(current_level()>0){ AppVolData.intvlVI = 100; vollp->mode=0; }}
    
   
    //if(in_strong_branching ){
    //  std::cout<<"setAuxiliaryInfo "<<lp->getAuxiliaryInfo()<<std::endl;
    //lp->setAuxiliaryInfo(new MCND_parent_branch_data());
    //}
}

//#############################################################################
double
MCND_lp::compute_lower_bound(const double old_lower_bound,
                             const BCP_lp_result& lpres,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts){
    
    //std::cout<<"compute lower bound"<<std::endl;
    const int tc = lpres.termcode();
    if (tc & BCP_ProvenOptimal)
        return lpres.objval();
    
    
    if (tc & BCP_DualObjLimReached)
        return best_soln.cost + 1e-5;
    
    
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
    /*
     getLpProblemPointer()->user_has_lp_result_processing = true;
     NetDesign mip(topo, data);
     mip.create_model(0,0,false);
     mip.solve();
     double * xy = new double[data->narcs+data->narcs*data->ndemands];
     std::cout<<"solution test: "<<mip.getSolution(xy)<<std::endl;
     delete [] xy;*/
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

//-------------------------------------------------------------------------------------------

void
MCND_lp::generate_cuts_in_lp(const BCP_lp_result& lpres,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts,
                             BCP_vec<BCP_cut*>& new_cuts,
                             BCP_vec<BCP_row*>& new_rows){
    
    int sz = cover_manager.covers.sizeOfCollection;
    int newly_added = cover_manager.gend > sz? sz : cover_manager.gend;
    Cover *vi = cover_manager.covers.end;
    std::cout<<"generate cuts: "<<newly_added<<" "<<sz<<" cutsvecsz: "<<cuts.size()<<std::endl;
    for(int i=newly_added;i--;){
        new_cuts.push_back(new CoverCut(vi));
        vi = vi->prev;
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
    std::cout<<"cuts_to_rows: "<<cutnum<<std::endl;
    for (int i=cutnum; i--;) {
        BCP_row * row = new BCP_row();
        row->LowerBound(1);
        row->UpperBound(1e40);
        rows.push_back(row);
    }
}
//#############################################################################
//#############################################################################
//=============================================================================

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
	std::cout<<"slect cuts "<<before_fathom<<std::endl;
	cover_manager.gend=0;
    BCP_lp_user::select_cuts_to_delete(lpres,vars,cuts,before_fathom, deletable);
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
//=============================================================================

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
    // if not feas. Check it for real.
    if((lp_result.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
        std::deque<int> topo;
        for (int a=data.narcs; a--;)
            if(!(vars[a]->lb()==0 && vars[a]->ub()==0)) topo.push_back(a);
        MCND_solution* sol_ = SolveMinCostFLow(topo);
        topo.clear();
        if(sol_){
            std::cout<<"vol stopped but there is sol: "<<sol_->objective_value()<<std::endl;
            //update_branch_data(ndata);
            //cut_off = process_flow_solution(vars, sol_);
            return sol_;
        }
        return sol_ ;
    }else{ // if very feas. Check it for real.
        update_branch_data(ndata);
        double viol=0;
        const double * lhs = lp_result.lhs();
        for(int d=data.nnodes*data.ndemands;d--;)
            viol += std::abs(lhs[d]);
        
        if(viol<0.1 ){
            std::deque<int> topo;
            for (int a=data.narcs; a--;)
                if(!(vars[a]->lb()==0 && vars[a]->ub()==0)) topo.push_back(a);
            MCND_solution* sol_ = SolveMinCostFLow(topo);
            topo.clear();
            std::cout<<"little viol: SolveMinCostFLow"<<std::endl;
            if(sol_){
                cut_off = process_flow_solution(vars, sol_);
                return sol_;
            }
        }
        
        return 0;
    }
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::restore_feasibility(const BCP_lp_result& lpres,
                             const std::vector<double*> dual_rays,
                             const BCP_vec<BCP_var*>& vars,
                             const BCP_vec<BCP_cut*>& cuts,
                             BCP_vec<BCP_var*>& vars_to_add,
                             BCP_vec<BCP_col*>& cols_to_add){
    
    std::cout<<"restore_feasibility"<<std::endl;
}

//#############################################################################

BCP_solution*
MCND_lp::generate_heuristic_solution(const BCP_lp_result& lpres,
                                     const BCP_vec<BCP_var*>& vars,
                                     const BCP_vec<BCP_cut*>& cuts){
    return 0;
    std::cout<<"try heuristic "<<std::endl;
    OsiVolSolverInterface * s = getOsiVolBabSolver();
    const double * x = lpres.x();
    const double * u = lpres.pi();
    CoinWarmStartDual * hs = new CoinWarmStartDual(data.ndemands*data.nnodes, u);
    std::deque<int> test_topo;
    for (int a=data.narcs; a--;){
        if(vars[a]->lb()==1.0){
            test_topo.push_back(a);
            std::cout<<test_topo.back()<<" ";
        }else if(vars[a]->ub()==1.0){
            double r = rand()%101/100.0;
            if(r <= 1){ test_topo.push_back(a);
                std::cout<<test_topo.back()<<"* ";}
        }
    }
    std::cout<<std::endl;
    s->direct_solve(test_topo,hs);
    
    delete hs;
    std::cout<<"ub: "<<s->getObjValue()<<"  UB:"<<upper_bound()<<std::endl;
    if(s->getObjValue() < upper_bound()){
        MCND_solution* sol_;
        sol_ = SolveMinCostFLow(test_topo);
        
        if(sol_!=0){
            std::cout<<"new ub: "<<sol_->cost<<std::endl;
            for(int a=test_topo.size();a--;){
                bool flag=false;
                int arc =  test_topo[a];
                for (int k = 0; k < data.ndemands; ++k){
                    if(sol_->xy[data.narcs+k*data.narcs+arc]){
                        flag=true; break;
                    }
                }
                if(!flag){
                    sol_->xy[arc] = 0.0;
                    sol_->cost -= data.arcs[arc].f;
                    std::cout<<"cut: "<<arc<<std::endl;
                }else sol_->xy[arc] = 1.0;
            }
            std::cout<<"oficial new ub: "<<sol_->cost<<std::endl;
        }
        test_topo.clear();
        return sol_;
    }
    
    return 0;
}


//#############################################################################

void
MCND_lp::pack_feasible_solution(BCP_buffer& buf, const BCP_solution* sol){
    //std::cout<<"pack feas sol: "<<std::endl;
    const MCND_solution* mcndsol = dynamic_cast<const MCND_solution*>(sol);
    
    mcndsol->pack(buf);
    
}

//-------------------------------------------------------------------------

OsiVolSolverInterface*
MCND_lp::getOsiVolBabSolver(){
    return dynamic_cast<OsiVolSolverInterface*>(getLpProblemPointer()->lp_solver);
}


//#############################################################################
