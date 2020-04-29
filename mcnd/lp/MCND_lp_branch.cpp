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
    
    //std::cout<<"select_branching_candidates "<<lp_mode<<std::endl;

    //return BCP_DoNotBranch_Fathomed;

    /*if(current_level() == 2){
        return BCP_DoNotBranch_Fathomed;
	}*/
	lp_mode &= ~LP_LogicalFixed;
	if(lp_mode & LP_ForceNodeAbort){
		//std::cout<<"node abort"<<std::endl;
		//lp_mode = LP_Normal;
        return BCP_DoNotBranch_Fathomed;
	}
	
    
    if(local_cut_pool.size()>0 ){
    	return BCP_DoNotBranch;
	}else if((lp_mode & LP_CutAddedFromHeuristic)){
		return BCP_DoNotBranch;
	} 
	mapd.clear();
    cover_manager.covers.map_collection(mapd);
	localc_manager.locals.map_collection(mapd);
	globalc_manager.globals.map_collection(mapd);

	
    std::deque<Pair2> candidates;
	BCP_vec<double> vbd(4, 0.0);
    BCP_vec<int> vpos(1, 0);
    double psc0, psc1;
    const double * psol = lpres.x();
    const double * rcsol = lpres.dj();

    int max_cand = 5 ;
    if(candidates.empty()){
		for (int a=data.narcs; a--;) {
			//if(psol[a]==0) continue;
			if(vars[a]->lb()==0 && vars[a]->ub()==1 ){
				/*/if(min(ninsp[a].fst, ninsp[a].snd) >= data.narcs*0.1){
					psc0 = psol[a]*psdcost[a].fst/double(ninsp[a].fst);
					psc1 = (1.0-psol[a])*psdcost[a].snd/double(ninsp[a].snd);
				}*/
				if(pump_heur.branch_candidates[a]>0){
					psc0 = psc1 = pump_heur.branch_candidates[a];
				}
				else{
					psc0 =  psol[a];
					psc1 =  (1.0-psol[a]);
				} 
				/*if(has_sol){
					psc0 = abs(best_sol.xy[a] - psol[a]);
					psc1 = abs(best_sol.xy[a] - psol[a]);
				}*/
				candidates.push_back(Pair2(a, fmin(psc0, psc1)));
				//candidates.push_back(Pair2(a, min(psol[a]-0.7)));
			}
		}
	}
	 
	std::stable_sort(candidates.begin(), candidates.end(), compPair2());
	int arc;
	int ncands=0;
	
	while(!candidates.empty() && ncands<max_cand){
		arc = candidates.front().fst;
		std::cout<<"candidate "<<arc<<" "<<psol[arc]<<" "<<std::setprecision(10)<<candidates.front().snd<<std::endl;
		///*<<" test: "<<fmin(psc0,psc1)<<", "<<min(ninsp[arc].fst, ninsp[arc].snd)*/<<std::endl;
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
	pump_heur.branch_candidates.assign(data.narcs, -1);
	
	//std::cout<<"do branch "<<ncands<<" force: "<<force_branch<<std::endl;
    if(ncands){
     	lp_mode |= LP_StrongBranch;
     	return BCP_DoBranch;
    }else{
    	//std::cout<<"NO CAND"<<std::endl;
    	//lp_mode = LP_Normal;
    	return BCP_DoNotBranch_Fathomed;
	}
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::logical_fixing(const BCP_lp_result& lpres,
		   const BCP_vec<BCP_var*>& vars,
		   const BCP_vec<BCP_cut*>& cuts,
		   const BCP_vec<BCP_obj_status>& var_status,
		   const BCP_vec<BCP_obj_status>& cut_status,
		   const int var_bound_changes_since_logical_fixing,
		   BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd){
		   
	if(lp_mode & LP_ForceNodeAbort) return;

	//std::cout<<"logical_fixing "<<vars.size()<<std::endl;
    //return;
    const double* psol = lpres.x();
    const double* dsol = lpres.pi();
    const double * rcsol = lpres.dj();
    double lb = lpres.objval();
    double gij, ckij;
     
    std::fill(yfix,yfix+data.narcs, -1);
 	for (int a=data.narcs; a--;){
		gij = rcsol[a];
		if(vars[a]->lb()==0.0 && vars[a]->ub()==1.0){
			//std::cout<<"penalty test: "<<a<<" "<<gij<<" lbs: "<<lb<<" "<<LBi<<std::endl;
			if(gij>0 && (lb+gij)>=upper_bound()){
				std::cout<<a<<" LOGFIX 0 ("<<lb<<" + "<<gij<<") ="<<(lb+gij)<<" "<<upper_bound()<<std::endl;
				changed_pos.push_back(a);
 				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				yfix[a]=0; 
				
			}else if(gij<0 && (lb-gij)>=upper_bound()){
				std::cout<<a<<" LOGFIX FIX 1 ("<<lb<<" "<<gij<<") ="<<(lb-gij)<<" "<<upper_bound()<<std::endl;
				changed_pos.push_back(a);
				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
				yfix[a]=1;
			}
		}else if(vars[a]->ub()==0.0){ yfix[a]=0;
		}else if(vars[a]->lb()==1.0){ yfix[a]=1;}
	}
	
	if(!cut_varfix_and_updt( vars, cuts, changed_pos, new_bd)) return;

	for(int a=data.narcs; a--;){
		gij = rcsol[a];
		if(yfix[a]==0){ if(vars[a]->ub()==1.0)lp_mode |= LP_TestConnectivity; continue;}
		Arc * item  =  &data.arcs[a];
		double tt=0;
		for (int k=data.ndemands; k--;){
			int id = data.narcs+k*data.narcs+a;
			if(vars[id]->ub()==0.0) continue;
			ckij = item->c[k] - dsol[k*data.nnodes + item->j-1] + dsol[k*data.nnodes + item->i-1];
			if(ckij>1e-4){
				tt = ((gij > 0)&& (yfix[a]<1)) ? (lb+ gij) : lb;
				if(tt + ckij*vars[id]->ub() > (upper_bound()+1e-4) ){
					double bd = (upper_bound() - tt)/ckij;
					if(bd<1e-8)bd=0.0;
					if(bd > vars[id]->ub()){ std::cout<<"flow logical fix enlarged ub "<<bd<<" "<<vars[id]->ub()<<std::endl; abort();}
					else if(vars[id]->ub()-bd <= 1e-4) continue;
					//std::cout<<"flowlogfix: "<<id<<" "<<bd<<" "<<vars[id]->ub()<<" rc: "<<ckij<<std::endl;
					changed_pos.push_back(id);
					new_bd.push_back(0.0);
					new_bd.push_back(bd);
					lp_mode |= LP_TestConnectivity;
				}
			}
		}
	}
	
	
	if(changed_pos.size()>0) lp_mode |= LP_LogicalFixed;
	if(lp_mode & LP_TestConnectivity){
		lp_mode &= ~LP_TestConnectivity;
		if(!verify_feasibility( changed_pos, new_bd, changed_pos.size())){
			lp_mode |= LP_ForceNodeAbort;
		}
 	}

}

//-------------------------------------------------------------------------------------------

bool 
MCND_lp::cut_varfix_and_updt(const BCP_vec<BCP_var*>& vars, 
							const BCP_vec<BCP_cut*>& cuts,
							BCP_vec<int>& var_changed_pos,
                            BCP_vec<double>& var_new_bd){
                                    
	int contfixed=-1;
	bool viol;
	bool zrofx;
	bool conn_arcfix;
	const int cutnum = cuts.size();
	GlobalCut * gloc;
	while(contfixed < int(var_changed_pos.size())){
     	contfixed = int(var_changed_pos.size());
    	zrofx=false;
    	conn_arcfix=false;
    	for (int i = getLpProblemPointer()->core->cutnum(); i < cutnum; ++i) {
			MCND_Cut * cut = dynamic_cast<MCND_Cut*>(cuts[i]);
			if(cut->purgbl()) continue;
			if(!cut->check_viol_updt_fix(vars,  var_changed_pos,  var_new_bd, viol, zrofx,  yfix)){
				//std::cout<<"insatisfeasible: "<<std::endl; cut->print();
				var_changed_pos.clear();
				var_new_bd.clear();
				lp_mode |= LP_ForceNodeAbort;
				return false;
			}
			if(!viol){ cut->mark_purgbl(); }//std::cout<<"no viol: "<<std::endl; cut->print();}
			//else cut->mark_unpurgbl();
			if(contfixed < int(var_changed_pos.size())){
				contfixed = int(var_changed_pos.size());
				i = getLpProblemPointer()->core->cutnum();
			}
		}
		
		gloc = globalc_manager.globals.begin;
		for(int i = globalc_manager.globals.sizeOfCollection; i--;){
			//std::cout<<"gloc_updt "<<gloc->serial_nmbr<<" "<<gloc->purgbl<<std::endl;
			if(gloc->purgbl) continue;
			if(!gloc->check_viol_updt_fix(vars,  var_changed_pos,  var_new_bd, viol, zrofx,  yfix)){
				//std::cout<<"insatisfeasible: "<<std::endl; gloc->print();
				var_changed_pos.clear();
				var_new_bd.clear();
				lp_mode |= LP_ForceNodeAbort;
				return false;
			}
			if(!viol){ gloc->purgbl=true; }//std::cout<<"no viol: "<<std::endl; cut->print();}
			//else cut->mark_unpurgbl();
			if(contfixed < int(var_changed_pos.size())){
				contfixed = int(var_changed_pos.size());
				i = globalc_manager.globals.sizeOfCollection;
				gloc = globalc_manager.globals.begin;
			}else gloc = gloc->next;
		}
		
		if(zrofx){
			if(!flwconnect.check_connectivity(vars, var_changed_pos, var_new_bd, conn_arcfix, yfix)){
				var_changed_pos.clear();
				var_new_bd.clear();
				//std::cout<<"fathom node flwconnect.check_connectivity"<<std::endl;
				lp_mode |= LP_ForceNodeAbort;
				return false;
			}	
		} 
		if(!conn_arcfix) break;
    }
    return true;
}


//-------------------------------------------------------------------------------------------

BCP_branching_object_relation
MCND_lp::compare_branching_candidates(BCP_presolved_lp_brobj* newobj,
				    BCP_presolved_lp_brobj* oldobj)
{

	int var_new = newobj->candidate()->forced_var_pos->front();
    //std::cout<<"compare_branching_candidates"<<std::endl;
    
    const BCP_lp_result & child0  = newobj->lpres(0);
    const BCP_lp_result & child1  = newobj->lpres(1);

    int sz = getOsiVolBabSolver()->getNumRows();
    double diff0 = (child0.objval() - LBi);
    double diff1 = (child1.objval() - LBi);
   	double score = fmin(diff0, diff1);
    newobj->user_data()[0] = new MCND_node_branch_data(sz, score , var_new, 0, LBi, true);
    newobj->user_data()[1] = new MCND_node_branch_data(sz, score , var_new, 1, LBi, false);
    //std::cout<<"branch0 "<<child0.objval() <<" branch1 "<<child1.objval()<<std::endl;
    //std::cout<<" comparing ("<<var_new<<") y: "<<y[var_new]<<" score: "<<score<<" mode: "<<lp_mode<<std::endl;
     int infeas=0;
     if((child0.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
     	to_logical_fix.push_back(Pair2(var_new,1.0));
     	++infeas;
     }
     if((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
     	(child1.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
     	to_logical_fix.push_back(Pair2(var_new,0.0));
     	lp_mode |= LP_TestConnectivity;
      	++infeas;
	 }
	 
	if(infeas==2){
		//std::cout<<"OPOR"<<std::endl; //abort();
		to_logical_fix.clear();
		lp_mode |= LP_ForceNodeAbort ;
		return BCP_NewPresolvedIsBetter_BranchOnIt;
	}
	 if(oldobj==0){
    	return BCP_NewPresolvedIsBetter;
    }
    
    int var_old = oldobj->candidate()->forced_var_pos->front();
    MCND_node_branch_data* old_data = dynamic_cast<MCND_node_branch_data*>(oldobj->user_data()[0]);
    double score_old = old_data->score;
    

     
	if(score > (score_old+1e-4) ){
	  return BCP_NewPresolvedIsBetter;
	}else if(abs(score - score_old) < 1e-4 && y[var_new] > y[var_old]){
		return BCP_NewPresolvedIsBetter;
	}return BCP_OldPresolvedIsBetter;
    
   
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::set_user_data_for_children(BCP_presolved_lp_brobj* best,
                                    const int selected){
    
    //BCP_lp_branching_object * cand = best->candidate();
    //std::cout<<" cand children: "<<cand->child_num<<std::endl;
    if(lp_mode & LP_ForceNodeAbort){
		//lp_mode = LP_Normal;
		return;
	}
	
    //std::cout<<" set child data "<<std::endl;

    BCP_vec< BCP_user_data * >& childs_data = best->user_data();

    //int var_ = *best->candidate()->forced_var_pos->begin();
    
    const BCP_lp_result & child0  = best->lpres(0);
    const BCP_lp_result & child1  = best->lpres(1);
    
    
    MCND_node_branch_data * cdata;
    //stock  for further investigation WarmStartDual(getNumRows(), dual, actv);
     
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[0]);
    //if(!(lp_mode & LP_LogicalFixed) && y[cdata->branch_var]==0) cdata->reduced_run =true;
    cdata->hs = new WarmStartDual(cdata->dual_size, child0.pi(), mapd); //std::cout<<"dualsz: "<<cdata->dual_size<<std::endl;
    cdata->parent = current_index();
    
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[1]);
    if(lp_mode & LP_TestConnectivity) cdata->test_conn=true;
    cdata->hs = new WarmStartDual(cdata->dual_size, child1.pi(), mapd);//std::cout<<"dualsz: "<<cdata->dual_size<<std::endl;
    cdata->parent = current_index();

	lp_mode = LP_Normal;
}

//-------------------------------------------------------------------------------------------

void
MCND_lp::set_actions_for_children(BCP_presolved_lp_brobj* best){
	int nvars = best->candidate()->vars_affected();
	BCP_vec< int >&  vars = *(best->candidate()->forced_var_pos);
	int var_branch = vars[0];
	BCP_vec< BCP_child_action >& childs_action = best->action();
	if(lp_mode & LP_ForceNodeAbort){
		childs_action[0] = BCP_FathomChild;
		childs_action[1] = BCP_FathomChild;
		//std::cout<<"cancel branching "<<std::endl;
		return;
	}
	std::cout<<"branching variable: "<<var_branch<<std::endl;
	
	if(to_logical_fix.size()) lp_mode |= LP_LogicalFixed;
	strong_branch_var_logicfix(best->candidate());
	//if(best->candidate()->forced_var_pos->size()>1)
		//std::cout<<"YES!! "<<best->candidate()->forced_var_pos->size()<<std::endl;
	
	const BCP_lp_result & child0  = best->lpres(0);
	const BCP_lp_result & child1  = best->lpres(1);
	

	bool feas0 = !((child0.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf);
	bool feas1 = !((child1.termcode() & BCP_PrimalObjLimReached) == BCP_PrimalObjLimReached ||
		(child1.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf);
		
	int ret = verify_children_feasibility( best->candidate(),  feas0, feas1 );
	
	double ybar_branch = child0.x()[var_branch];
    double diff0 = (child0.objval() - LBi);
    double diff1 = (child1.objval() - LBi);
    
	if(!feas0){
		childs_action[0] = BCP_FathomChild;
	}else{
		/*++ninsp[var_branch].fst;
		if((ybar_branch)>1e-8) psdcost[var_branch].fst+=diff0/ybar_branch;
		else psdcost[var_branch].fst+=diff0*1e8;*/
		childs_action[0] = BCP_ReturnChild;
	}
	
	double ybar_branch_ = (1.0- ybar_branch);
	if(!feas0 && feas1){
		/*++ninsp[var_branch].snd;
		if(ybar_branch_>1e-8)psdcost[var_branch].snd+=diff1/ybar_branch_;
		else psdcost[var_branch].snd+=diff1*1e8;*/
		childs_action[1] = BCP_KeepChild;
	}else if(!feas1){
		childs_action[1] = BCP_FathomChild;
	}else{
		/*++ninsp[var_branch].snd;
		if(ybar_branch_>1e-8)psdcost[var_branch].snd+=diff1/ybar_branch_;
		else psdcost[var_branch].snd+=diff1*1e8;*/
		childs_action[1] = BCP_ReturnChild;
	}
	//std::cout<<" 0-side childs_action "<<childs_action[0]<<std::endl;
	//std::cout<<" 1-side childs_action "<<childs_action[1]<<std::endl;
	
}


//========================================================================================

bool 
MCND_lp::verify_feasibility(const BCP_vec<int> & vars_chngd, const BCP_vec<double> & chngd_bd, int nvars){
	OsiVolSolverInterface* lpsolver = getOsiVolBabSolver();
	
	for(int n=nvars;n--;){
		lpsolver->setColUpper(vars_chngd[n], chngd_bd[n*2+1]);
	} 

	const double * collb = lpsolver->getColLower();
    const double * colub = lpsolver->getColUpper();
    
    int ret = lpfeaschecker.solve_feas(collb, colub);
    if(ret<0) return false;
    else return true;
}

//========================================================================================


int 
MCND_lp::verify_children_feasibility(BCP_lp_branching_object * candidate, bool& feas0, bool& feas1 ){
	int ret = 0;
	const BCP_vec< int >& vars = *candidate->forced_var_pos;
	const BCP_vec< double >& vars_bd = *candidate->forced_var_bd;
	if(feas0){
		candidate->apply_child_bd(getOsiVolBabSolver(), 0);
		//std::cout<<"testfeas 0"<<std::endl;
		feas0 = verify_feasibility(vars,vars_bd,0);
		ret +=1;
	}
	if(!feas1) return ret;
	
	if(lp_mode & LP_TestConnectivity){
 		candidate->apply_child_bd(getOsiVolBabSolver(), 1);
		feas1 = verify_feasibility(vars,vars_bd,0);
		ret +=2;
	}	
	return ret;
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
			//std::cout<<"FIX: "<<p->fst<<" to "<<p->snd<<std::endl;
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
		//std::cout<<"childone lb: "<<aux[strt1]<<" ub: "<<aux[strt1+1]<<std::endl;
		bd_vars[strt] = aux[strt1];
		bd_vars[strt+1] = aux[strt1+1];
		strt+=2;
	}
	aux.clear();
}




