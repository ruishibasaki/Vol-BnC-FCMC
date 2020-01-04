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
	LBi = lpres.objval();
    return BCP_DoNotBranch_Fathomed;

    if(current_level() == 0){
        return BCP_DoNotBranch_Fathomed;

	}
	
	
	if(LBi<0 || LBi>best_soln.cost){
		std::cout<<"LBi<0 || LBi>best_UB"<<std::endl;
		return BCP_DoNotBranch_Fathomed;
	}
	
	if(cut_off){
		std::cout<<"test cut off"<<std::endl;
		return BCP_DoNotBranch_Fathomed;
	}
	

    OsiVolSolverInterface * s = getOsiBabSolver();
    MCND_node_branch_data* nodedata = dynamic_cast<MCND_node_branch_data*>( get_user_data());
    
    freq.assign(data.narcs,0.0);
    double div = s->getIterationCount();
   // std::cout<<"div: "<<div<<" freq: "<<freq[201]<<" "<<freq[139]<<" "<<freq[89]<<" "<<freq[216]<<std::endl;
    if(nodedata!=0){
        //for(int a=data.narcs; a--;){ freq[a] += nodedata->hist[a];}
           // if(nodedata->hist[a])std::cout<<"ahist: "<<a<<" "<<nodedata->hist[a]std::endl; }
        div= 1;
    }
    //std::cout<<"div: "<<div<<" freq: "<<freq[201]<<" "<<freq[139]<<" "<<freq[89]<<" "<<freq[216]<<std::endl;

	BCP_vec<double> vbd(4, 0.0);
    BCP_vec<int> vpos(1, 0);
    const double * psol = lpres.x();
    std::list<HeapCell> order;
	for (int a=data.narcs; a--;) {
        y[a] = psol[a];
        freq[a] /= div;  if(freq[a] > 1.0) freq[a] = 1.0;
		if(vars[a]->lb()==0 && vars[a]->ub()==1){
			order.push_back(HeapCell(a,freq[a]));
			
		}
	}
	order.sort(comp());
	int arc;
	int ncands=0;
	while(!order.empty() && ncands<4){
		arc = order.back().k;
        order.pop_back();
        ++freq_cand[arc];
		std::cout<<"candidate "<<arc<<" "<<freq[arc]<< " "<<y[arc]<<" "<<freq_cand[arc]<<std::endl;
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
	order.clear();
	std::cout<<"do branch "<<ncands<<" force: "<<force_branch<<std::endl;
    if(ncands) return BCP_DoBranch;
    else return BCP_DoNotBranch_Fathomed;
	//for (int a=data.narcs; a--;) {
		//if(vars[a]->lb()==0 && vars[a]->ub()==1){
			//y[a] = psol[a];
			//test =(1.0-psol[a]);
			//if(test>=0 ){
				//vpos[0] = a;
				//vbd[0] = 0.0;
				//vbd[1] = 0.0;
				//vbd[2] = 1.0;
				//vbd[3] = 1.0;
				//cands.push_back(new  BCP_lp_branching_object(2,
									  //0, 0, /* vars/cuts_to_add */
									  //&vpos, 0, &vbd, 0, /* forced parts */
									  //0, 0, 0, 0 /* implied parts */));
			//}
		//}
	//}
	//std::cout<<topo->szunfxd<<" rdo branch "<<cands.size()<<" force: "<<force_branch<<std::endl;
	//if(cands.size()) return BCP_DoBranch;
	//return BCP_DoNotBranch_Fathomed;
}

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
    return;
    const double* psol = lpres.x();
    const double * dsol = lpres.pi();
    double gij;
    
    if(current_level() == 0){
        for (int a=data.narcs; a--;){
            if(psol[a]>=(1.0-1e-30)){
                std::cout<<"logical fix "<<a<<"  "<<psol[a]<<std::endl;
                changed_pos.push_back(a);
                new_bd.push_back(1.0);
                new_bd.push_back(1.0);
                continue;
            }else if(psol[a]<=1e-30){
                std::cout<<"logical fix "<<a<<"  "<<psol[a]<<std::endl;
                changed_pos.push_back(a);
                new_bd.push_back(0.0);
                new_bd.push_back(0.0);
                continue;
            }
        }
    }else
        for (int a=data.narcs; a--;){
            //std::cout<<"logical fix "<<a<<" "<<freq[a]<<"  "<<psol[a]<<std::endl;
            if(vars[a]->lb()==0 && vars[a]->ub()==1){
                gij = lag_subproblem(a, dsol);
                //std::cout<<"penalty test: "<<a<<std::endl;
                if(gij>0 && (LBi+gij)>=best_soln.cost){
                    std::cout<<a<<" WILL FIX 0"<<std::endl;
                    changed_pos.push_back(a);
                    new_bd.push_back(0.0);
                    new_bd.push_back(0.0);
                }else if(gij<0 && (LBi-gij)>=best_soln.cost){
                    std::cout<<a<<" WILL FIX 1"<<std::endl;
                    changed_pos.push_back(a);
                    new_bd.push_back(1.0);
                    new_bd.push_back(1.0);
                /*}else if(psol[a]>=(1.0)){
                    std::cout<<"logical fix "<<a<<"  "<<psol[a]<<std::endl;
                    changed_pos.push_back(a);
                    new_bd.push_back(1.0);
                    new_bd.push_back(1.0);
                    continue;*/
                }else if(psol[a]<=1e-30 && freq_cand[a]>=3){
                    std::cout<<"logical fix "<<a<<"  "<<psol[a]<<" "<<freq[a]<<std::endl;
                    changed_pos.push_back(a);
                    new_bd.push_back(0.0);
                    new_bd.push_back(0.0);
                    continue;
                }
            }
        }
}


BCP_branching_object_relation
MCND_lp::compare_branching_candidates(BCP_presolved_lp_brobj* newobj,
				    BCP_presolved_lp_brobj* oldobj)
{

    std::cout<<"compare_branching_candidates"<<std::endl;
    int var_ = *newobj->candidate()->forced_var_pos->begin();
    double pscore, nscore;

    const BCP_lp_result & child0  = newobj->lpres(0);
    const BCP_lp_result & child1  = newobj->lpres(1);

    pscore = (child1.objval()-LBi); pscore = (pscore>=0)?pscore: 0;
    nscore = (child0.objval()-LBi);
    if(nscore>0 && !((child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf))
        nscore = (nscore>=0)?nscore: 0;
    else nscore=0;
  
    std::cout<<" comparing ("<<var_<<") "<<freq[var_]<<" (+) "<<pscore<<" (-) "<<nscore<<std::endl;

    double fp = ((1.0-y[var_])>1e-10) ? (1.0-y[var_]) : 0;
    double fn = (y[var_]>1e-10)? y[var_] : 0 ;
    
    if(fp>0) pscore /= fp; else pscore = pscore*1e10;
    if(fn>0) nscore /= fn; else nscore = nscore*1e10;

    pscore *= (freq[var_]);
    nscore *= (freq[var_]);

       //stock si for further investigation
  

    double newscore = nscore+pscore; //score(nscore,pscore,0.33);

    std::cout<<" comparing ("<<var_<<") "<<newscore<<" (+) "<<pscore<<" (-) "<<nscore<<std::endl;
    //std::cout<<" <> "<<newscore<<"( "<<pscore<<" , "<<nscore<<" )"<<std::endl;
    
    //OsiSolverInterface * s = getLpProblemPointer()->lp_solver;
    //OsiVolAuxInfo * pdata = static_cast<OsiVolAuxInfo*>(s->getApplicationData());
   
    //std::cout<<"oi "<<pdata<<std::endl;
    newobj->user_data()[0] = new MCND_node_branch_data( newscore, nscore, var_, 0, 0);
    //std::cout<<"oi2"<<pdata<<std::endl;
    newobj->user_data()[1] = new MCND_node_branch_data( newscore, pscore, var_, 1,0);
    
    //here compare phi
    if(oldobj==0){
      return BCP_NewPresolvedIsBetter;
    }else{
       MCND_node_branch_data* cdata= dynamic_cast<MCND_node_branch_data *>(oldobj->user_data()[0]);
        //std::cout<<"cdata: "<<cdata<<std::endl;
      if(newscore > cdata->score_parent){
          return BCP_NewPresolvedIsBetter;
      }else return BCP_OldPresolvedIsBetter;
    }
   
}


void
MCND_lp::set_user_data_for_children(BCP_presolved_lp_brobj* best,
                                    const int selected){
    
    //BCP_lp_branching_object * cand = best->candidate();
    //std::cout<<" cand children: "<<cand->child_num<<std::endl;
    
    std::cout<<" set child data "<<std::endl;

    BCP_vec< BCP_user_data * >& childs_data = best->user_data();
    //int var_ = *best->candidate()->forced_var_pos->begin();
    
    const BCP_lp_result & child0  = best->lpres(0);
    const BCP_lp_result & child1  = best->lpres(1);
    
    MCND_node_branch_data * cdata;
    //stock  for further investigation

    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[0]);
    cdata->hs = new CoinWarmStartDual(data.ndemands*data.nnodes, child0.pi());
    cdata= dynamic_cast<MCND_node_branch_data *>(childs_data[1]);
    cdata->hs = new CoinWarmStartDual(data.ndemands*data.nnodes, child1.pi());
    
    //best->user_data().push_back();

}


void
MCND_lp::set_actions_for_children(BCP_presolved_lp_brobj* best)
{
  int vars = best->candidate()->vars_affected();
  BCP_vec< int > & var = *(best->candidate()->forced_var_pos);
  for(;vars--;) std::cout<<"branching variable: "<<var[vars]<<std::endl;
  
  const BCP_lp_result & child0  = best->lpres(0);
  const BCP_lp_result & child1  = best->lpres(1);
  BCP_vec< BCP_child_action >& childs_action = best->action();
  
  if((child0.termcode() & BCP_ProvenPrimalInf) == BCP_ProvenPrimalInf){
	childs_action[0] = BCP_FathomChild;
  }else childs_action[0] = BCP_ReturnChild;
  
  childs_action[1] = BCP_ReturnChild;
  std::cout<<" 0-side BCP_FathomChild/BCP_ReturnChild "<<childs_action[0]<<std::endl;
}

//#############################################################################
//#############################################################################
//=============================================================================

void
MCND_lp::select_vars_to_delete(const BCP_lp_result& lpres,
				   const BCP_vec<BCP_var*>& vars,
				   const BCP_vec<BCP_cut*>& cuts,
				   const bool before_fathom,
				   BCP_vec<int>& deletable)
{
	//std::cout<<"slect vars"<<std::endl;
    //select_vars_to_delete(lpres,vars,cuts,before_fathom, deletable);
}

//=============================================================================

void
MCND_lp::select_cuts_to_delete(const BCP_lp_result& lpres,
				   const BCP_vec<BCP_var*>& vars,
				   const BCP_vec<BCP_cut*>& cuts,
				   const bool before_fathom,
				   BCP_vec<int>& deletable)
{
	//std::cout<<"slect cuts"<<std::endl;
    //select_cuts_to_delete(lpres,vars,cuts,before_fathom, deletable);
}

void
MCND_lp::purge_slack_pool(const BCP_vec<BCP_cut*>& slack_pool,
		     BCP_vec<int>& to_be_purged){
				 
	//std::cout<<"try to purge"<<std::endl;	 
				 
}

