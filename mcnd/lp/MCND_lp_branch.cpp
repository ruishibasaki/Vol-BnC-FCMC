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
    
    std::cout<<"select_branching_candidates"<<std::endl;

	//LBi = lpres.objval();
    //return BCP_DoNotBranch_Fathomed;

   /* if(current_level() == 2){
        return BCP_DoNotBranch_Fathomed;
	}*/
	if(abort){
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
    int max_cand = 1;
    if(lp_mode==LP_Normal || candidates.empty()){
		for (int a=data.narcs; a--;) {
			if(vars[a]->lb()==0 && vars[a]->ub()==1){
				candidates.push_back(Pair2(a, abs(psol[a]-0.5)));
			}
		}
	}else max_cand = 1;
	
	std::stable_sort(candidates.begin(), candidates.end(), compPair2());
	int arc;
	int ncands=0;
	
	while(!candidates.empty() && ncands<max_cand){
		arc = candidates.front().fst;
		std::cout<<"candidate "<<arc<<" "<<y[arc]<<std::endl;
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
	std::cout<<"logical_fixing "<<std::endl;
    //return;
    const double* psol = lpres.x();
    const double * rcsol = lpres.dj();
    double lb = lpres.objval();
    double gij;
    
    
	for (int a=data.narcs; a--;){
		//std::cout<<"logical fix "<<a<<" "<<freq[a]<<"  "<<psol[a]<<std::endl;
		if(vars[a]->lb()==0 && vars[a]->ub()==1){
			gij = rcsol[a];
			//std::cout<<"penalty test: "<<a<<" "<<gij<<std::endl;
			if(gij>0 && (lb+gij)>=upper_bound()){
				std::cout<<a<<" WILL FIX 0 ("<<lb<<" + "<<gij<<") ="<<(lb+gij)<<" "<<upper_bound()<<std::endl;
				changed_pos.push_back(a);
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
}

//-------------------------------------------------------------------------------------------

BCP_branching_object_relation
MCND_lp::compare_branching_candidates(BCP_presolved_lp_brobj* newobj,
				    BCP_presolved_lp_brobj* oldobj)
{

    std::cout<<"compare_branching_candidates"<<std::endl;
    int var_new = newobj->candidate()->forced_var_pos->front();
    const BCP_lp_result & child0  = newobj->lpres(0);
    const BCP_lp_result & child1  = newobj->lpres(1);

    std::cout<<" comparing ("<<var_new<<") y: "<<y[var_new]<<
    	" fst: "<<newobj->candidate()->forced_var_pos->front()<<", bd: "<<newobj->candidate()->forced_var_bd->front()<<
    	" snd: "<<newobj->candidate()->forced_var_pos->back()<<", bd: "<<newobj->candidate()->forced_var_bd->back()
    	<<std::endl;
   
    int sz = data.ndemands*data.nnodes+cover_manager.covers.sizeOfCollection;
    newobj->user_data()[0] = new MCND_node_branch_data(sz, y[var_new], 0, var_new, 0, 0);
    newobj->user_data()[1] = new MCND_node_branch_data(sz, y[var_new], 0, var_new, 1, 0);
    
     
     if((child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
     	return BCP_NewPresolvedIsBetter_BranchOnIt;
     }
		
    if(oldobj==0){
    	return BCP_NewPresolvedIsBetter;
    }
    int var_old = *oldobj->candidate()->forced_var_pos->begin();
    if(lp_mode==LP_DiveToFeasibility){
		if(y[var_new] > y[var_old]){
          	return BCP_NewPresolvedIsBetter;
      	}else return BCP_OldPresolvedIsBetter;
  	}else{
		if(y[var_new] > y[var_old]){
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
	int vars = best->candidate()->vars_affected();
	BCP_vec< int > & var = *(best->candidate()->forced_var_pos);
	for(;vars--;) std::cout<<"branching variable: "<<var[vars]<<std::endl;

	const BCP_lp_result & child0  = best->lpres(0);
	const BCP_lp_result & child1  = best->lpres(1);
	BCP_vec< BCP_child_action >& childs_action = best->action();
	
	bool zero_fathomed=false;
	if((child0.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
		childs_action[0] = BCP_FathomChild;
		zero_fathomed=true;
	}else if(lp_mode==LP_DiveToFeasibility){
		 childs_action[0] = BCP_KeepChild;
	}else childs_action[0] = BCP_ReturnChild;
	
	if(zero_fathomed && !((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached)){
		childs_action[1] = BCP_KeepChild;
	}else if((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached){
		childs_action[1] = BCP_FathomChild;
	}else childs_action[1] = BCP_ReturnChild;
	std::cout<<" 0-side childs_action "<<childs_action[0]<<std::endl;
	std::cout<<" 1-side childs_action "<<childs_action[1]<<std::endl;

}


//#############################################################################
//#############################################################################
//#############################################################################


double 
MCND_lp::lag_subproblem(int a, const double * u){
	double kpsack =0;
	double fillUp =0;
	double x, rc;

	std::list<HeapCell> heap;
	//get reduced cost for each commodity in arc e
    Arc * item = &data.arcs[a];
	for(int k=0;k<data.ndemands;++k){
		rc =  item->c[k] - u[k*data.nnodes + item->j-1] + u[k*data.nnodes + item->i-1];
		if(rc<0.0){
			heap.push_back(HeapCell(k,rc));
		}
	}  
	
	heap.sort(comp());
	//std::stable_sort(heap.begin(), heap.end(), comp());
	while(heap.size()>0){
		if(fillUp < item->capa){
			x = std::min((item->capa - fillUp),  data.d_k[heap.back().k].quantity);
			fillUp += x;
			kpsack += heap.back().rc_ * x;
			heap.pop_back();
			//std::cout<<"in"<<std::endl;
		}else{
			heap.pop_back();
		}	
	}
	return kpsack;
}



