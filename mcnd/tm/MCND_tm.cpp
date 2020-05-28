// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#include <fstream>

#include "CoinSearchTree.hpp"

#include "BCP_math.hpp"
#include "BCP_tm.hpp"
#include "MCND_tm.hpp"
#include "MCND_init.hpp"
#include "MCND_branch_score.hpp"

//#############################################################################

int main(int argc, char* argv[])
{
	srand(11);
    MCND_initialize mcnd_init;
    return bcp_main(argc, argv, &mcnd_init);
}

//#############################################################################

void
MCND_tm::pack_module_data(BCP_buffer& buf, BCP_process_t ptype)
{

  switch (ptype) {
    case BCP_ProcessType_LP:
    //std::cout<<"try to pack data to lp "<<data.ndemands<<" "<<data.narcs<<" "<<data.nnodes<<std::endl;
          data.pack(buf);
          buf.pack(getTmProblemPointer()->has_ub());
          init_sol.pack(buf);
    break;
  default:
abort();
  }
}

//#############################################################################

BCP_solution*
MCND_tm::unpack_feasible_solution(BCP_buffer& buf)
{
	
	//std::cout<<"unpack feas solu"<<std::endl;
   MCND_solution* new_sol = new MCND_solution;
   new_sol->unpack(buf);
   //if (new_sol->objective_value() > best_soln.objective_value())
   //init_sol = *new_sol;
      	//std::cout<<"done"<<std::endl;

   return new_sol;
}

//#############################################################################

void
MCND_tm::initialize_core(BCP_vec<BCP_var_core*>& vars,
		       BCP_vec<BCP_cut_core*>& cuts,
		       BCP_lp_relax*& matrix){
  //std::cout<<"tm_init_CORE "<<vars.size()<<" "<<cuts.size()<<std::endl;
  int narcs=data.narcs;
  int ndemands=data.ndemands;
  int nnodes=data.nnodes;
  int contvar=0;
  int controw=0;
  int elems=0;
  
  int  corecols = narcs + narcs*ndemands;
  int  corerows = nnodes*ndemands;

  vars.reserve(corecols); 
  for (int a = 0; a < narcs; ++a){
    vars.unchecked_push_back(new BCP_var_core(BCP_BinaryVar,
					      data.arcs[a].f, 0, 1));
  }
  
  for (int k = 0; k < ndemands; ++k)
	for (int a = 0; a < narcs; ++a){
		 vars.unchecked_push_back(new BCP_var_core(BCP_ContinuousVar,
					      data.arcs[a].c[k], 0, data.arcs[a].b[k]));
	}
	     
  
  cuts.reserve(corerows);
  for (int i = 0; i < nnodes; ++i)
	for (int k = 0; k < ndemands; ++k){
		if(i==data.d_k[k].O-1){
			cuts.unchecked_push_back(new BCP_cut_core(-data.d_k[k].quantity,-data.d_k[k].quantity));
		}else if(i==data.d_k[k].D-1){
			cuts.unchecked_push_back(new BCP_cut_core(data.d_k[k].quantity,data.d_k[k].quantity));
			
		}else{
			 cuts.unchecked_push_back(new BCP_cut_core(0,0));			
		 }
		 
	}
 
}

//#############################################################################

void
MCND_tm::create_root(BCP_vec<BCP_var*>& added_vars,
		   BCP_vec<BCP_cut*>& added_cuts,
		   BCP_user_data*& user_data)
{
  //std::cout<<"no CG or VG"<<std::endl;
  return;
}

//-----------------------------------------------------------------------------

void
MCND_tm::init_new_phase(int phase,
		 BCP_column_generation& colgen,
		 CoinSearchTreeBase*& candidates){
	
	//std::cout<<"init phase "<<phase<<std::endl;
	/*if(phase){
		colgen = BCP_DoNotGenerateColumns_Fathom;
		candidates = new CoinSearchTree<CoinSearchTreeCompareBreadth>;
    }else{
		colgen = BCP_DoNotGenerateColumns_Send;
		candidates = new CoinSearchTree<CoinSearchTreeCompareBreadth>;
	}*/
    
	colgen = BCP_DoNotGenerateColumns_Fathom;
    candidates = new CoinSearchTree<CoinSearchTreeCompareLowerBound>;
}

//#############################################################################
// managing
  
void
MCND_tm::change_candidate_heap(CoinSearchTreeManager& candidates,
			const bool new_solution){
	
	//std::cout<<"MCND_tm:: "<<candidates.getTree()->compName()<<std::endl;
	//std::cout<<"cand size "<<candidates.size()<<" "<<candidates.getTree()->compName()<<std::endl;
	//Best_LB=lower_bound();
	BCP_tm_prob *p = getTmProblemPointer();
	//if(!new_solution){
		CoinSearchTreeBase * tree = candidates.getTree();
		CoinSearchTreeBase *t = new CoinSearchTree<CoinSearchTreeCompareLowerBound>;
		//std::cout<<"change candidate heap size:"<<tree->size()<<std::endl;
        MCND_node_branch_data * user_data;
		BCP_tm_node * n;
		CoinTreeNode ** add = new CoinTreeNode*;
		while(tree->size()){
			*add = tree->top();
			n = dynamic_cast<BCP_tm_node*>(*add);
            //user_data =  dynamic_cast<MCND_node_branch_data*>(n->_data._user.GetRawPtr());
			if(n->status > BCP_PrunedNode_Discarded ||
			   n->status < BCP_PrunedNode_OverUB){
				//std::cout<<"tree: "<<n->getQuality()<<" "<<(*add)->getTrueLB()<<" "<<lower_bound()<<std::endl;
				t->push(1, add);
			}else{
				const double oldTrueLB = floor((*add)->getTrueLB()*p->lb_multiplier);
    			p->lower_bounds.erase(oldTrueLB);
			} 
			tree->pop();
		}
		candidates.setTree(t);
		delete add;
		//std::cout<<"final size "<<t->size()<<std::endl;
	//}
	
	//BCP_tm_user::change_candidate_heap(candidates,new_solution);
}
/*
//-----------------------------------------------------------------------------

CoinSearchTreeBase *
MCND_tm::createSearchType(const char * type){
	
	if(strcmp(type, "CoinSearchTreeCompareBest") == 0)
		return new CoinSearchTree<CoinSearchTreeCompareBest>;
	else if(strcmp(type, "CoinSearchTreeCompareBreadth") == 0)
		return new CoinSearchTree<CoinSearchTreeCompareBreadth>;
	else if(strcmp(type, "CoinSearchTreeCompareDepth") == 0)
		return new CoinSearchTree<CoinSearchTreeCompareDepth>;
	else 
		return new CoinSearchTree<CoinSearchTreeComparePreferred>; 
	
}
*/
//#############################################################################

void
MCND_tm::process_message(BCP_buffer& buf){

	if(Best_LB<lower_bound())
		Best_LB=lower_bound();
}

//-----------------------------------------------------------------------------

void
MCND_tm::display_final_information(const BCP_lp_statistics& lp_stat){
	BCP_tm_user::display_final_information(lp_stat);
    if(Best_LB<lower_bound())
        Best_LB=lower_bound();
	
	double ub = upper_bound();
	double t = double( clock() - t_start ) / double( CLOCKS_PER_SEC );
	if(t<time_lim-0.01) Best_LB = ub;
	if(Best_LB > ub) Best_LB = ub;
	//std::cout<<"The global lower bound: "<<Best_LB<<" / "<<lower_bound()<<std::endl;
    std::ofstream file("fileout", std::ios::app);
    file<<std::setprecision(10)<<instance<<" lb: "<<Best_LB<<" ub: "<<ub<<" gap: "<<(ub-Best_LB)/ub*100<<" nodes: "<<getTmProblemPointer()->search_tree.processed()
    <<" t: "<<t<<std::endl;
    file.close();
    
    
}

//-----------------------------------------------------------------------------

void
MCND_tm::display_node_information(BCP_tree& search_tree,
			   const BCP_tm_node& node,
			   bool after_processing_node){
	if(Best_LB<lower_bound())
		Best_LB=lower_bound();
	
	BCP_buffer buf;
	buf.pack(lower_bound());
	std::cout<<"The lower bound: "<<Best_LB<<std::endl;
	broadcast_message(BCP_ProcessType_LP,buf);
	buf.clear();
}

//-----------------------------------------------------------------------------

void
MCND_tm::display_feasible_solution(const BCP_solution* sol)
{
	//if(Best_LB<lower_bound())
		//Best_LB=lower_bound();
	//std::cout<<"The global lower bound: "<<lower_bound()<<std::endl;	
	//Best_LB=lower_bound();
  /*const MCND_solution* mcnd_sol = dynamic_cast<const MCND_solution*>(sol);
  if (! mcnd_sol) {
    throw BCP_fatal_error("\
MCND_tm::display_feasible_solution() invoked with non-MCND_solution.\n");
  }


  if (tm_par.entry(MCND_tm_par::DisplaySolutionSignature)) {
    mcnd_sol->display(tm_par.entry(MCND_tm_par::SolutionFile));
  }*/
  //std::cout<<"try to display solution"<<std::endl;
}
