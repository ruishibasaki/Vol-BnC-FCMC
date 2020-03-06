// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cstdio>
#include <algorithm>

#include "CoinTime.hpp"


#include "MCND_lp.hpp"

#include "BCP_math.hpp"
#include "BCP_lp.hpp"
#include "BCP_lp_node.hpp"

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

     if((child0.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
     	to_logical_fix.push_back(Pair2(var_new,1.0));
     }else if((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached){
     	to_logical_fix.push_back(Pair2(var_new,0.0));
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
    
   // std::cout<<" set child data "<<std::endl;

    BCP_vec< BCP_user_data * >& childs_data = best->user_data();
    //int var_ = *best->candidate()->forced_var_pos->begin();
    
    const BCP_lp_result & child0  = best->lpres(0);
    const BCP_lp_result & child1  = best->lpres(1);
    
    MCND_node_branch_data * cdata;
    //stock  for further investigation WarmStartDual(getNumRows(), dual, actv);
    
	//int sz = data.ndemands*data.nnodes+cover_manager.covers.sizeOfCollection;
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[0]);
    //cdata->tofix = to_logical_fix;
    cdata->hs = new WarmStartDual(cdata->dual_size, child0.pi(), mapd); //std::cout<<cdata->dual_size<<" "<<sz<<std::endl;
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[1]);
    //cdata->tofix = to_logical_fix;
    cdata->hs = new WarmStartDual(cdata->dual_size, child1.pi(), mapd);//std::cout<<cdata->dual_size<<" "<<sz<<std::endl;
    //to_logical_fix.clear();
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::set_actions_for_children(BCP_presolved_lp_brobj* best){
	int nvars = best->candidate()->vars_affected();
	BCP_vec< int >&  vars = *(best->candidate()->forced_var_pos);	
	std::cout<<"branching variable: "<<vars[0]<<std::endl;
	
	strong_branch_var_logicfix(best->candidate());
	if(best->candidate()->forced_var_pos->size()>1)
		std::cout<<"YES!! "<<best->candidate()->forced_var_pos->size()<<std::endl;
	
	const BCP_lp_result & child0  = best->lpres(0);
	const BCP_lp_result & child1  = best->lpres(1);
	BCP_vec< BCP_child_action >& childs_action = best->action();
	
	bool feas0 = !((child0.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf);
	bool feas1 = !((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child1.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf);
		
	verify_children_feasibility( best->candidate(),  feas0, feas1 );
	
	if(!feas0){
		childs_action[0] = BCP_FathomChild;
	}else{
		if(lp_mode==LP_DiveToFeasibility){
			childs_action[0] = BCP_KeepChild;
		}else childs_action[0] = BCP_ReturnChild;
	}
	
	if(!feas0 && feas1){
		childs_action[1] = BCP_KeepChild;
	}else if(!feas1){
		childs_action[1] = BCP_FathomChild;
	}else childs_action[1] = BCP_ReturnChild;
	
	std::cout<<" 0-side childs_action "<<childs_action[0]<<std::endl;
	std::cout<<" 1-side childs_action "<<childs_action[1]<<std::endl;
	lp_mode = LP_Normal;
}


//========================================================================================

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

//========================================================================================


void 
MCND_lp::verify_children_feasibility(BCP_lp_branching_object * candidate, bool& feas0, bool& feas1 ){
	if(feas0){
		candidate->apply_child_bd(getOsiVolBabSolver(), 0);
		std::cout<<"testfeas 0"<<std::endl;
		feas0 = verify_feasibility(*(candidate->forced_var_pos),0);
	}
	if(!feas1) return;
	const BCP_vec< int >& vars = *candidate->forced_var_pos;
	const BCP_vec< double >& vars_bd = *candidate->forced_var_bd;
	int sz_implied = vars.size();
	for(;sz_implied--;){
		std::cout<<"testfeas 1 "<<vars[sz_implied]<<" "<<vars_bd[sz_implied*4+3]<<" "<<vars.size()<<" "<<vars_bd.size()<<std::endl;
		if(vars[sz_implied]<data.narcs && vars_bd[sz_implied*4+3]<0.5){
			//candidate->apply_child_bd(getOsiVolBabSolver(), 0);
			candidate->apply_child_bd(getOsiVolBabSolver(), 1);
			feas1 = verify_feasibility(*(candidate->forced_var_pos),0);
			break;
		}
	}	
	
}


//========================================================================================

void
MCND_lp::strong_branch_var_logicfix(BCP_lp_branching_object* candidate){
	BCP_vec< int >* extra_vars = candidate->forced_var_pos;
	BCP_vec< double >* bd_extra_vars = candidate->forced_var_bd;
	int cand = (*candidate->forced_var_pos)[0];
	int old_sz = extra_vars->size();
	int added = 0;
	Pair2* p;   
	
	while(!to_logical_fix.empty()){
		p = &to_logical_fix.back();
		if(p->fst != cand){
			std::cout<<"MORE FIX: "<<p->fst<<" "<<p->snd<<std::endl;
			extra_vars->push_back(p->fst);
			++added;
			for(int i=4;i--;)bd_extra_vars->push_back(p->snd);
		}
		to_logical_fix.pop_back();
	}
	rearrange_bd_vec( old_sz, added, *bd_extra_vars );
}

//========================================================================================

void
MCND_lp::rearrange_bd_vec(int old_sz, int added, BCP_vec< double >& bd_vars ){
	if(!added) return;
	int sz = old_sz+added;
	int old_ = old_sz;
	int strt, strt1;
	double lb,ub;
	BCP_vec< double > aux;
	aux.reserve(bd_vars.size());
	 
	strt = 2*old_sz;
	for(int i=0;i<old_sz;++i){
		aux.unchecked_push_back(bd_vars[strt]);
		aux.unchecked_push_back(bd_vars[strt+1]);
		strt+=2;
	}
	
	for(int i=0;i<added;++i){
		strt = 4*old_sz + 4*i;
		strt1 = 2*old_;
		lb = bd_vars[strt];
		ub = bd_vars[strt+1];
		aux.unchecked_push_back(bd_vars[strt+2]);
		aux.unchecked_push_back(bd_vars[strt+3]);
		bd_vars[strt1] = lb;
		bd_vars[strt1+1] = ub;
		++old_;
	}
	
	strt = 2*old_;
	for(int i=0;i<sz;++i){
		strt1 = 2*i;
		std::cout<<"childone lb: "<<aux[strt1]<<" ub: "<<aux[strt1+1]<<std::endl;
		bd_vars[strt] = aux[strt1];
		bd_vars[strt+1] = aux[strt1+1];
		strt+=2;
	}
	aux.clear();
}




