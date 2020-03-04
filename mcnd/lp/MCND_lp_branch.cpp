// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cstdio>
#include <algorithm>

#include "CoinTime.hpp"


#include "MCND_lp.hpp"

#include "BCP_math.hpp"
#include "BCP_lp.hpp"
#include "BCP_lp_node.hpp"
#include "BCP_lp_functions.hpp"

//#############################################################################

BCP_branching_decision
MCND_lp::select_branching_candidates(const BCP_lp_result& lpres,
				   const BCP_vec<BCP_var*>& vars,
				   const BCP_vec<BCP_cut*>& cuts,
				   const BCP_lp_var_pool& local_var_pool,
				   const BCP_lp_cut_pool& local_cut_pool,
				   BCP_vec<BCP_lp_branching_object*>& cands,
				   bool force_branch){
    
    std::cout<<"select_branching_candidates "<<lp_mode<<std::endl;

    //return BCP_DoNotBranch_Fathomed;

   /* if(current_level() == 2){
        return BCP_DoNotBranch_Fathomed;
	}*/
	
	if(lp_mode == LP_ForceNodeAbort){
		lp_mode = LP_Normal;
        return BCP_DoNotBranch_Fathomed;
	}
	
	mapd.clear();
    cover_manager.covers.map_collection(mapd);
    OsiVolSolverInterface * s = getOsiVolBabSolver();
    
    if(local_cut_pool.size()>0){
		//lp_mode = LP_CutAdded;
		return BCP_DoNotBranch;
	}
	
    
	BCP_vec<double> vbd(4, 0.0);
    BCP_vec<int> vpos(1, 0);
    const double * psol = lpres.x();
    int max_cand = 3 ;
    if(lp_mode==LP_Normal || candidates.empty()){
		for (int a=data.narcs; a--;) {
			if(vars[a]->lb()==0 && vars[a]->ub()==1){
				candidates.push_back(Pair2(a, -abs(psol[a]-0.7)));
			}
		}
	}else max_cand = 1;
	
	std::stable_sort(candidates.begin(), candidates.end(), compPair2());
	int arc;
	int ncands=0;
	
	while(!candidates.empty() && ncands<max_cand){
		arc = candidates.front().fst;
		std::cout<<"candidate "<<arc<<" "<<y[arc]<<" "<<psol[arc]<<std::endl;
        candidates.pop_front();		
		vpos[0] = arc;
		vbd[0] = 0.0;
		vbd[1] = 0.0;
		vbd[2] = 1.0;
		vbd[3] = 1.0;
		cands.push_back(new  BCP_lp_branching_object(2,
							  0, 0, /* vars/cuts_to_add */
							  &vpos, 0, &vbd, 0, /* forced parts */
							  0, 0, 0, 0 /* implied parts */));
		
		++ncands;
	}
	candidates.clear();
	std::cout<<"do branch "<<ncands<<" force: "<<force_branch<<std::endl;
    if(ncands) return BCP_DoBranch;
    else return BCP_DoNotBranch_Fathomed;
	
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::logical_fixing(const BCP_lp_result& lpres,
		   const BCP_vec<BCP_var*>& vars,
		   const BCP_vec<BCP_cut*>& cuts,
		   const BCP_vec<BCP_obj_status>& var_status,
		   const BCP_vec<BCP_obj_status>& cut_status,
		   const int var_bound_changes_since_logical_fixing,
		   BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd)
{
	std::cout<<"logical_fixing "<<var_bound_changes_since_logical_fixing<<std::endl;
    //return;
    const double* psol = lpres.x();
    const double * rcsol = lpres.dj();
    double * yarcs = new double [data.narcs];
    double lb = lpres.objval();
    double gij;
    
    for (int a=data.narcs; a--;) yarcs[a] = 0;
    const int cutnum = cuts.size();
	for (int i = getLpProblemPointer()->core->cutnum(); i < cutnum; ++i) {
		CoverCut * cut = dynamic_cast<CoverCut*>(cuts[i]);
		cut->check_logical_fix( vars, yarcs);
	}
    
    BCP_vec<int> changed_to0;
    changed_to0.reserve(data.narcs);
	for (int a=data.narcs; a--;){
		if(yarcs[a]==1){
			std::cout<<"logical fix "<<a<<"  "<<psol[a]<<std::endl;
			changed_pos.push_back(a);
			new_bd.push_back(1.0);
			new_bd.push_back(1.0);
		}else if(vars[a]->lb()==0 && vars[a]->ub()==1){
			gij = rcsol[a];
			//std::cout<<"penalty test: "<<a<<" "<<gij<<" lbs: "<<lb<<" "<<LBi<<std::endl;
			if(gij>0 && (lb+gij)>=upper_bound()){
				std::cout<<a<<" WILL FIX 0 ("<<lb<<" + "<<gij<<") ="<<(lb+gij)<<" "<<upper_bound()<<std::endl;
				changed_pos.push_back(a);
				changed_to0.unchecked_push_back(a);
				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				
			}else if(gij<0 && (lb-gij)>=upper_bound()){
				std::cout<<a<<" WILL FIX 1 ("<<lb<<" "<<gij<<") ="<<(lb-gij)<<" "<<upper_bound()<<std::endl;
				changed_pos.push_back(a);
				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
			}/*else if(psol[a]>=(1.0)){
				std::cout<<"logical fix "<<a<<"  "<<psol[a]<<std::endl;
				changed_pos.push_back(a);
				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
				continue;
			}else if(psol[a]<=1e-30 ){
				std::cout<<"logical fix "<<a<<"  "<<psol[a]<<std::endl;
				changed_pos.push_back(a);
				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				continue;
			}*/
		}
	}
	if(changed_to0.size()>0){
		if(!verify_feasibility( changed_to0, changed_to0.size()))
			lp_mode = LP_ForceNodeAbort;
		changed_to0.clear();
	}
	delete [] yarcs;
}

//-------------------------------------------------------------------------------------------

BCP_branching_object_relation
MCND_lp::compare_branching_candidates(BCP_presolved_lp_brobj* newobj,
				    BCP_presolved_lp_brobj* oldobj)
{

	int var_new = newobj->candidate()->forced_var_pos->front();
    std::cout<<"compare_branching_candidates"<<std::endl;
    
    const BCP_lp_result & child0  = newobj->lpres(0);
    const BCP_lp_result & child1  = newobj->lpres(1);

    int sz = data.ndemands*data.nnodes+cover_manager.covers.sizeOfCollection;
   	double score = std::min((child0.objval() - LBi), (child1.objval() - LBi));
    newobj->user_data()[0] = new MCND_node_branch_data(sz, score , var_new, 0, LBi);
    newobj->user_data()[1] = new MCND_node_branch_data(sz, score , var_new, 1, LBi);
    
    std::cout<<" comparing ("<<var_new<<") y: "<<y[var_new]<<" score: "<<score<<" mode: "<<lp_mode<<std::endl;

     if((child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
     	return BCP_NewPresolvedIsBetter_BranchOnIt;
     }
		
    if(oldobj==0){
    	return BCP_NewPresolvedIsBetter;
    }
    int var_old = oldobj->candidate()->forced_var_pos->front();
    MCND_node_branch_data* old_data = dynamic_cast<MCND_node_branch_data*>(oldobj->user_data()[0]);
    double score_old = old_data->score;
    

    if(lp_mode==LP_DiveToFeasibility){
		if(y[var_new] > y[var_old]){
          	return BCP_NewPresolvedIsBetter;
      	}else return BCP_OldPresolvedIsBetter;
  	}else{
		if(score > score_old){
		  return BCP_NewPresolvedIsBetter;
		}else return BCP_OldPresolvedIsBetter;
    }
   
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::set_user_data_for_children(BCP_presolved_lp_brobj* best,
                                    const int selected){
    
    //BCP_lp_branching_object * cand = best->candidate();
    //std::cout<<" cand children: "<<cand->child_num<<std::endl;
    
    //std::cout<<" set child data "<<std::endl;

    BCP_vec< BCP_user_data * >& childs_data = best->user_data();
    //int var_ = *best->candidate()->forced_var_pos->begin();
    
    const BCP_lp_result & child0  = best->lpres(0);
    const BCP_lp_result & child1  = best->lpres(1);
    
    MCND_node_branch_data * cdata;
    //stock  for further investigation WarmStartDual(getNumRows(), dual, actv);
    
	//int sz = data.ndemands*data.nnodes+cover_manager.covers.sizeOfCollection;
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[0]);
    cdata->hs = new WarmStartDual(cdata->dual_size, child0.pi(), mapd); //std::cout<<cdata->dual_size<<" "<<sz<<std::endl;
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[1]);
    cdata->hs = new WarmStartDual(cdata->dual_size, child1.pi(), mapd);//std::cout<<cdata->dual_size<<" "<<sz<<std::endl;
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::set_actions_for_children(BCP_presolved_lp_brobj* best){
	int nvars = best->candidate()->vars_affected();
	BCP_vec< int > & vars = *(best->candidate()->forced_var_pos);
	for(int v=nvars;v--;) std::cout<<"branching variable: "<<vars[v]<<std::endl;

	const BCP_lp_result & child0  = best->lpres(0);
	const BCP_lp_result & child1  = best->lpres(1);
	BCP_vec< BCP_child_action >& childs_action = best->action();
	
	bool zero_fathomed=false;
	if((child0.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
		childs_action[0] = BCP_FathomChild;
		zero_fathomed=true;
	}else{
		if(!verify_feasibility(vars, nvars)){
			childs_action[0] = BCP_FathomChild;
			zero_fathomed=true;
		}else if(lp_mode==LP_DiveToFeasibility){
			childs_action[0] = BCP_KeepChild;
		}else childs_action[0] = BCP_ReturnChild;
	}
	
	if(zero_fathomed && !((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached)){
		childs_action[1] = BCP_KeepChild;
	}else if((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached){
		childs_action[1] = BCP_FathomChild;
	}else childs_action[1] = BCP_ReturnChild;
	std::cout<<" 0-side childs_action "<<childs_action[0]<<std::endl;
	std::cout<<" 1-side childs_action "<<childs_action[1]<<std::endl;
	lp_mode = LP_Normal;
}


//=======================================================================================

bool 
MCND_lp::verify_feasibility(const BCP_vec<int> & vars_chngd, int nvars){
	OsiVolSolverInterface* lpsolver = getOsiVolBabSolver();
	
	for(;nvars--;){
		lpsolver->setColUpper(vars_chngd[nvars], 0.0);
	} 
	
	const double * collb = lpsolver->getColLower();
    const double * colub = lpsolver->getColUpper();
    
    int ret = lpfeaschecker.solve(collb, colub);
    if(ret<0) return false;
    else return true;
}


