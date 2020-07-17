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

    
	if(current_level() == 0 && current_iteration()<10 && !no_heur) return BCP_DoNotBranch;
	/*if(current_level() == 0){
        return BCP_DoNotBranch_Fathomed;
	}*/
	
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
	
	/*if(no_gap_reduct >5 && !break_diving){
    	force_dive=true;
    }else{
    	break_diving = false;
    	force_dive=false;
    }*/
    
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

    //if(no_gap_reduct >5) max_cand=1;
    if(reduced_run)max_cand=1;
    else max_cand=1;
 

	if(candidates.empty()){
		for (int a=data.narcs; a--;) {
			if(vars[a]->lb()==0 && vars[a]->ub()==1 ){
				psc0 =  psol[a];
				psc1 =  (1.0-psol[a]);  psc1 = psc1<0 ? 0: psc1;
				/*if(fmin(psc0, psc1)<0.1) candidates.push_back(Pair2(a, fmin(psc0, psc1)));
				else*/ candidates.push_back(Pair2(a, fmin(psc0, psc1)*data.arcs[a].capa));
			}
		}
	}
		
    std::stable_sort(candidates.begin(), candidates.end(), compPair2());
	  
	int arc;
	int ncands=0;
	
	while(!candidates.empty() && ncands<max_cand){
		arc = candidates.front().fst;
		//std::cout<<"candidate "<<arc<<" "<<psol[arc]<<" "<<std::setprecision(10)<<candidates.front().snd<<std::endl;
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
 	lp_mode &= ~LP_LogicalFixed;
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
	if(current_iteration()>1 && !(lp_mode & LP_TighterBounds)) return;

	//std::cout<<"logical_fixing "<<lpchecker.solved<<std::endl;
    //return;
    const double* psol = lpres.x();
    const double* dsol = lpres.pi();
    const double * rcsol = lpres.dj();
    double lb = lpres.objval();
    double gij, ckij;
    bool arcfx=false;
    unfix=0;
    std::fill(yfix,yfix+data.narcs, -1);
    double ub = upper_bound();
    double pres=0;
 	for (int a=data.narcs; a--;){
		if(vars[a]->lb()==0.0 && vars[a]->ub()==1.0 && yfix[a]==-1){
			//std::cout<<"penalty test: "<<a<<" rc: "<<gij<<" y: "<<psol[a]<<std::endl;
			gij = rcsol[a];
			if(gij>0 && (lb+gij)>=ub){
				std::cout<<a<<" LOGFIX 0 ("<<lb<<" + "<<gij<<") ="<<(lb+gij)<<" "<<ub<<std::endl;
				changed_pos.push_back(a);
 				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				yfix[a]=0;
				arcfx=true;
				continue;
			}else if(gij<0 && (lb-gij)>=ub){
				std::cout<<a<<" LOGFIX FIX 1 ("<<lb<<" "<<gij<<") ="<<(lb-gij)<<" "<<ub<<std::endl;
				changed_pos.push_back(a);
				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
				yfix[a]=1;
				arcfx=true;
				continue;
			} 
			if(current_level()==0){pres=0;}
			else{pres=1e-30*(pow(10,current_level())); if(pres>1e-4)pres=1e-4;}
			if(psol[a]<pres){
				std::cout<<a<<" SOLUTION FIX 0 "<<psol[a]<<std::endl;
				changed_pos.push_back(a);
 				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				yfix[a]=0;
				arcfx=true;
			}else if(psol[a]>1.0-pres){
				std::cout<<a<<" SOLUTION FIX 1 "<<psol[a]<<std::endl;
				changed_pos.push_back(a);
 				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
				yfix[a]=1;
				arcfx=true;
			}
 
		}else if(vars[a]->ub()==0.0){ yfix[a]=0;
		}else if(vars[a]->lb()==1.0){ yfix[a]=1;}
	}
	if(lpchecker.solved) 
 		if(lpchecker.logical_yfix(ub, changed_pos, new_bd, yfix)) arcfx=true;

	if(arcfx)
		if(!cut_varfix_and_updt( vars, cuts, changed_pos, new_bd)){ 
			if(flwconnect.redone)flwconnect.clear_memory(); 
			return;
		}
 	
 	if(lpchecker.solved) 
 		lpchecker.logical_xfix(ub, yfix, vars);
 	
	double tt, dk;
	double bdlb, bdub;
	double lpbd;
	int id;
	Arc * item;
	std::vector<double> lbtigh(data.ndemands,0.0);
	
	for(int a=data.narcs; a--;){
		if(yfix[a]==0){ if(vars[a]->ub()==1.0){lp_mode |= LP_TestConnectivity;} continue;}
		else if(yfix[a]==-1) unfix++;
		
		item  =  &data.arcs[a];
		gij = rcsol[a];
		lbtigh.assign(data.ndemands,0.0);
		if(vars[a]->lb()==1.0)
			flow_lb_fixing( lbtigh,  vars, psol, dsol, rcsol,  a,  lb,  ub);

		for (int k=data.ndemands; k--;){
			id = data.narcs+k*data.narcs+a;
			dk = data.d_k[k].quantity;
 			if(vars[id]->ub()==0.0 || flwconnect.bound_red[k*data.narcs+a].trd == 0.0) continue;
 
			ckij = rcsol[id];
			
			bdub =-1;
			bdlb = fmax(vars[id]->lb(), lbtigh[k]); 
			if(ckij>1e-4){
				tt = ((gij > 0)&& (vars[a]->lb()==0.0)) ? (lb+ gij) : lb;
				if(tt + ckij*(vars[id]->ub()/dk - psol[id]) > (ub+1e-4) ){
					bdub = psol[id] + (ub - tt)/ckij;
					bdub*=dk;
				}
			} 
					
			lpbd = lpchecker.bound_red[a*data.ndemands+k];
			if(bdub >=0) bdub = lpbd>=0 ? fmin(lpbd, bdub): bdub;
			else if(lpbd>=0) bdub = lpbd;
			else if(bdub<0){bdub= vars[id]->ub();}
			
			if(bdlb<1e-4 && bdub<1e-4){bdlb = bdub=0;}
			if(flwconnect.bound_red[k*data.narcs+a].snd>bdub || bdlb>(bdub+1e-30)){
				lp_mode |= LP_ForceNodeAbort;
				changed_pos.clear();
				new_bd.clear();
				//std::cout<<"MCND_lp::logical_fixing::flwconnect.bd conflict "<<bdlb<<" "<<bdub<<std::endl;
				if(lp_mode & LP_TestConnectivity || flwconnect.redone) flwconnect.clear_memory();
				return;
			}
			
			if(vars[id]->ub()-bdub <= 1e-4 && bdlb-vars[id]->lb() <= 1e-4) continue;
			changed_pos.push_back(id);
			new_bd.push_back(bdlb);
			new_bd.push_back(bdub);
			flwconnect.bound_red[k*data.narcs+a].fst = changed_pos.size()-1;
			flwconnect.bound_red[k*data.narcs+a].snd = bdlb;
 			flwconnect.bound_red[k*data.narcs+a].trd = bdub;
			lp_mode |= LP_TestConnectivity;
		}
	}
	lbtigh.clear();
	
	if(lp_mode & LP_TestConnectivity || flwconnect.redone){
		flwconnect.reset(vars, yfix);
		int ret = flwconnect.check_flowbounds(vars, changed_pos, new_bd, yfix);
		flwconnect.clear_memory();
		if(ret<0){
			changed_pos.clear();
			new_bd.clear();
			//std::cout<<"MCND_lp::logical_fixing::flwconnect.check_upperbounds NOT"<<std::endl;
			lp_mode |= LP_ForceNodeAbort;
			return;
		}
		
	} 
	
	if(lpchecker.solved) {
		getOsiVolBabSolver()->reset_dualsol(lpchecker.dual_map);
		lpchecker.bound_red.assign(data.narcs*data.ndemands,-1);
	}
	
 	double gap= (ub  - LBi)/ub;
 	if( (gap<0.01) && !getOsiVolBabSolver()->has_checked){
  		//std::cout<<"getOsiVolBabSolver()->test_opposites"<<std::endl;
  		double prev = changed_pos.size();
 		getOsiVolBabSolver()->test_opposites( changed_pos, new_bd,  yfix, vars,  ub);
 		arcfx = (prev < changed_pos.size()) ?  true: arcfx;
 	} 
 	
	if(changed_pos.size()){
		lp_mode |= LP_LogicalFixed;	
		//lp_mode |= LP_TestFeasibility;
	}
	/*for (int a=changed_pos.size(); a--;){
		if(vars[changed_pos[a]]->ub() < new_bd[a*2+1])std::cout<<"ub "<<vars[changed_pos[a]]->ub()<<" < "<<new_bd[a*2+1]<<std::endl;
		if(vars[changed_pos[a]]->lb() > new_bd[a*2])std::cout<<"lb "<<vars[changed_pos[a]]->lb()<<" > "<<new_bd[a*2]<<std::endl;
	}*/
	
}

//-------------------------------------------------------------------------------------------

void 
MCND_lp::flow_lb_fixing(std::vector<double>& lbtigh, const BCP_vec<BCP_var*>& vars, const double * psol, 
						const double * dsol, const double * rcsol, int arc, double lb, double ub){
	int id, id2, idk;
	bool tight;
	double ckij, sum, dk;
	std::deque<Pair2> rcx;
	Arc * item  =  &data.arcs[arc];
 	for (int k=data.ndemands; k--;){
		ckij = rcsol[data.narcs+k*data.narcs+arc];
		if(ckij<0)rcx.push_back(Pair2(k, double(ckij/data.d_k[k].quantity)));
	}
	std::stable_sort(rcx.begin(), rcx.end(), compPair2());
	//std::cout<<"outra"<<std::endl;
	for (int k=rcx.size(); k--;){
		idk = rcx[k].fst;
		dk = data.d_k[idk].quantity;
		id = data.narcs+idk*data.narcs+arc;
		ckij = -rcsol[id];
		tight=false;
		if(lb + ckij*(psol[id]-vars[id]->lb()/dk) > (ub+1e-4)){
			if(k-1 < 0){ /*std::cout<<"OPA1 "<<id<<" "<<ckij<<" "<<psol[id]<<" "<<vars[id]->lb()<<std::endl;*/ tight=true;
			}else{ 
				sum=0; 
				for (int kk=0; kk<k;++kk){
					id2 = data.narcs+rcx[kk].fst*data.narcs+arc;
					sum += rcsol[id2]*(vars[id2]->ub()/dk - psol[id2]);
					if(sum<-1e-30) break;
				}
				if(sum>=0){
					//std::cout<<"OPA2 "<<sum<<std::endl;
					tight=true;
				}
			}
			if(tight){
				lbtigh[idk] = dk*(psol[id] - (ub-lb)/ckij);
			}
		}
	}
	rcx.clear();
	 
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
				//std::cout<<"MCND_lp::cut_varfix_and_updt::insatisfeasible: "<<std::endl; cut->print();
				if(cut->cut_type==1){
					globalc_manager.globalc_generation_main2(0, yfix);
				}
				var_changed_pos.clear();
				var_new_bd.clear();
				lp_mode |= LP_ForceNodeAbort;
				return false;
			}
			if(!viol){ cut->mark_purgbl(); }//std::cout<<"no viol: "<<std::endl; cut->print();}
			//else cut->mark_unpurgbl();
			if(contfixed < int(var_changed_pos.size())){
				contfixed = int(var_changed_pos.size());
				lp_mode |= LP_LogicalFixed;	
				i = getLpProblemPointer()->core->cutnum();
			}
		}
		
		gloc = globalc_manager.globals.begin;
		for(int i = globalc_manager.globals.sizeOfCollection; i--;){
			//std::cout<<"gloc_updt "<<gloc->serial_nmbr<<" "<<gloc->purgbl<<std::endl;
			if(gloc->purgbl) continue;
			if(!gloc->check_viol_updt_fix(vars,  var_changed_pos,  var_new_bd, viol, zrofx,  yfix)){
				//std::cout<<"MCND_lp::cut_varfix_and_updt::insatisfeasible: "<<std::endl; gloc->print();
				globalc_manager.globalc_generation_main2(0, yfix);
				var_changed_pos.clear();
				var_new_bd.clear();
				lp_mode |= LP_ForceNodeAbort;
				return false;
			}
			if(!viol){ gloc->purgbl=true; }//std::cout<<"no viol: "<<std::endl; cut->print();}
			//else cut->mark_unpurgbl();
			if(contfixed < int(var_changed_pos.size())){
				contfixed = int(var_changed_pos.size());
				lp_mode |= LP_LogicalFixed;	
				i = globalc_manager.globals.sizeOfCollection;
				gloc = globalc_manager.globals.begin;
			}else gloc = gloc->next;
		}
		
		if(zrofx){
 			int ret = flwconnect.check_connectivity(vars, var_changed_pos, var_new_bd, conn_arcfix, yfix);
			if(ret<0){
				if(ret==-1){ globalc_manager.globalc_generation_main2(0, yfix); }
				var_changed_pos.clear();
				var_new_bd.clear();
				//std::cout<<"MCND_lp::cut_varfix_and_updt::flwconnect.check_connectivity NOT"<<std::endl;
				lp_mode |= LP_ForceNodeAbort;
				return false;
			}	
		} 
		if(!conn_arcfix) break;
		else lp_mode |= LP_LogicalFixed;	
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
    newobj->user_data()[0] = new MCND_node_branch_data(sz, diff0 , var_new, 0, LBi, true);
    newobj->user_data()[1] = new MCND_node_branch_data(sz, diff1 , var_new, 1, LBi, false);
    //std::cout<<"var "<<var_new<<" branch0 "<<child0.objval() <<" branch1 "<<child1.objval()<<std::endl;
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
		//std::cout<<"MCND_lp::compare_branching_candidates:: BOTH INFEAS OR TOO HIGHT "<<var_new<<std::endl;
		to_logical_fix.clear();
		lp_mode |= LP_ForceNodeAbort ;
		return BCP_NewPresolvedIsBetter_BranchOnIt;
	}
	if(oldobj==0){
    	return BCP_NewPresolvedIsBetter;
    }
    
    double scoremin = fmin(diff0, diff1);
    double scoremax = max(diff0, diff1);
    int var_old = oldobj->candidate()->forced_var_pos->front();
    MCND_node_branch_data* old_data = dynamic_cast<MCND_node_branch_data*>(oldobj->user_data()[0]);
    diff0 = old_data->score;
    old_data = dynamic_cast<MCND_node_branch_data*>(oldobj->user_data()[1]);
    diff1 = old_data->score;
	double scoreminold = fmin(diff0, diff1);
    double scoremaxold = max(diff0, diff1);
    double pres = lower_bound*1e-4;
	if(scoremin > (scoreminold + pres) ){
	  return BCP_NewPresolvedIsBetter;
	}else if(abs(scoremin - scoreminold) < pres &&  scoremax > scoremaxold){
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
    if(best->action()[0]!=BCP_FathomChild){
    	cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[0]);
   		//if(!(lp_mode & LP_LogicalFixed) && y[cdata->branch_var]==0) cdata->reduced_run =true;
    	cdata->hs = new WarmStartDual(data.narcs, child0.x() ,cdata->dual_size, child0.pi(), mapd); //std::cout<<"dualsz: "<<cdata->dual_size<<std::endl;
    	cdata->parent = current_index();
    }
    if(best->action()[1]!=BCP_FathomChild){
    	cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[1]);
    	if(lp_mode & LP_TestConnectivity) cdata->test_conn=true;
    	cdata->hs = new WarmStartDual(data.narcs, child1.x(), cdata->dual_size, child1.pi(), mapd);//std::cout<<"dualsz: "<<cdata->dual_size<<std::endl;
    	cdata->parent = current_index();
    }
    
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
		
	verify_children_feasibility( best->candidate(),  feas0, feas1 );
	
	double ybar_branch = y[var_branch];
    double diff0 = (child0.objval() - LBi); if(diff0<0)diff0=0;
    double diff1 = (child1.objval() - LBi); if(diff1<0)diff1=0;
    double lbdev = (fmin(child0.objval(), child1.objval()) - lower_bound)/lower_bound;
    //std::cout<<"branching variable: "<<var_branch<<" lp: "<<(lpchecker.solved ? lpchecker.y_[var_branch]:-1.0)<<
    //" y: "<<y[var_branch]<<" nbranch: "<<min(ninsp[var_branch].fst, ninsp[var_branch].snd)<<" capa: "<<data.arcs[var_branch].capa<<std::endl;
    //std::cout<<"child0: "<<child0.objval()<<" child1: "<<child1.objval()<<" feas: "<<feas0<<", "<<feas1<<std::endl<<std::endl;
 	
   	if(dive_for_sol){++times_dived; std::cout<<"force dive for sol"<<std::endl;}
 
	if(!feas0){
		childs_action[0] = BCP_FathomChild;
	}else{
		++ninsp[var_branch].fst;
		if((ybar_branch)>1e-8) psdcost[var_branch].fst+=diff0/ybar_branch;
		else psdcost[var_branch].fst+=diff0*1e8;
		
		if(dive_for_sol && (y[var_branch]<=0.5)) childs_action[0] = BCP_KeepChild;
		else childs_action[0] = BCP_ReturnChild;
	}
	
	double ybar_branch_ = (1.0- ybar_branch);
	ybar_branch_ = ybar_branch_< 0 ? 0: ybar_branch_;
	if(!feas0 && feas1){
		++ninsp[var_branch].snd;
		if(ybar_branch_>1e-8)psdcost[var_branch].snd+=diff1/ybar_branch_;
		else psdcost[var_branch].snd+=diff1*1e8;
		
		childs_action[1] = BCP_KeepChild;
	}else if(!feas1){
		childs_action[1] = BCP_FathomChild;
	}else{
		++ninsp[var_branch].snd;
		if(ybar_branch_>1e-8)psdcost[var_branch].snd+=diff1/ybar_branch_;
		else psdcost[var_branch].snd+=diff1*1e8;
		
		if(dive_for_sol && (y[var_branch]>0.5)) childs_action[1] = BCP_KeepChild;
		else childs_action[1] = BCP_ReturnChild;
	}
	if(childs_action[1] == BCP_KeepChild || childs_action[0] == BCP_KeepChild) 
		lp_mode |= LP_Dive;
	//std::cout<<" 0-side childs_action "<<childs_action[0]<<std::endl;
	//std::cout<<" 1-side childs_action "<<childs_action[1]<<std::endl;
	
}


//========================================================================================

int 
MCND_lp::verify_feasibility(const BCP_vec<int> & vars_chngd, const BCP_vec<double> & chngd_bd, int nvars, bool feas_status){
	OsiVolSolverInterface* lpsolver = getOsiVolBabSolver();
	
	for(int n=nvars;n--;){
		lpsolver->setColUpper(vars_chngd[n], chngd_bd[n*2+1]);
	} 

	const double * collb = lpsolver->getColLower();
    const double * colub = lpsolver->getColUpper();

    int ret = lpfeaschecker.solve_feas(collb, colub, feas_status);
    if(ret == -2){
    	globalc_manager.globalc_generation_main2(colub, 0);
    }
    return ret;
}

//========================================================================================


int 
MCND_lp::verify_children_feasibility(BCP_lp_branching_object * candidate, bool& feas0, bool& feas1 ){
 	const BCP_vec< int >& vars = *candidate->forced_var_pos;
	const BCP_vec< double >& vars_bd = *candidate->forced_var_bd;
	int ret;
	if( (lp_mode & LP_TestConnectivity) /*feas1 &&*/){
		candidate->apply_child_bd(getOsiVolBabSolver(), 1);
		ret = verify_feasibility(vars,vars_bd,0, feas1);
		if(ret < 0){
			feas1 = false;
			if(ret<=-2){
				//std::cout<<"1 child infeasible"<<std::endl;
				feas0 = false;
				return 0;
			}
		}
	}
	//std::cout<<"0bracj "<<feas0<<std::endl;
 	//if(feas0){
	candidate->apply_child_bd(getOsiVolBabSolver(), 0);
	ret = verify_feasibility(vars,vars_bd,0, feas0);
	if(ret<0) feas0=false;
 	//}
 	//if(!feas0) std::cout<<"0 child infeasible"<<std::endl;
	return 0;
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
			std::cout<<"branchfix: "<<p->fst<<" to "<<p->snd<<std::endl;
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




