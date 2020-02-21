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
    MCND_initialize mcnd_init;
    return bcp_main(argc, argv, &mcnd_init);
}

//#############################################################################

void
MCND_tm::pack_module_data(BCP_buffer& buf, BCP_process_t ptype)
{
    
  switch (ptype) {
    case BCP_ProcessType_LP:
    std::cout<<"try to pack data to lp "<<data.ndemands<<" "<<data.narcs<<" "<<data.nnodes<<std::endl;
          data.pack(buf);
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
   if (new_sol->objective_value() > best_soln.objective_value())
      best_soln = *new_sol;
      	//std::cout<<"done"<<std::endl;

   return new_sol;
}

//#############################################################################

void
MCND_tm::initialize_core(BCP_vec<BCP_var_core*>& vars,
		       BCP_vec<BCP_cut_core*>& cuts,
		       BCP_lp_relax*& matrix)
{
  //std::cout<<"tm_init_CORE "<<vars.size()<<" "<<cuts.size()<<std::endl;
  int narcs=data.narcs;
  int ndemands=data.ndemands;
  int nnodes=data.nnodes;
  int contvar=0;
  int controw=0;
  int elems=0;
  
  int  corecols = narcs + narcs*ndemands;
  int  corerows = nnodes*ndemands;
  
  double* OBJ = new double[corecols];
  double* EV= new double[corerows*corecols];
  double* CLB = new double[corecols];
  double* CUB = new double[corecols];
  double* RLB = new double[corerows];
  double* RUB = new double[corerows];
  int* EI= new int[corerows*corecols];
  int* VB = new int[corerows+1]; /* VB- Starting positions of major-dimension vectors. */
  
  vars.reserve(corecols);
  for (int a = 0; a < narcs; ++a){
    vars.unchecked_push_back(new BCP_var_core(BCP_BinaryVar,
					      data.arcs[a].f, 0, 1));
	OBJ[contvar] = data.arcs[a].f;
	CLB[contvar] = 0;
	CUB[contvar] = 1;
	
	++contvar;
  }
  for (int a = 0; a < narcs; ++a)
	for (int k = 0; k < ndemands; ++k){
		 vars.unchecked_push_back(new BCP_var_core(BCP_ContinuousVar,
					      data.arcs[a].c[k], 0, BCP_DBL_MAX));
		 OBJ[contvar] =  data.arcs[a].c[k];
		 CLB[contvar] = 0;
		 CUB[contvar] = BCP_DBL_MAX;
		 ++contvar;
	}
					     
  
  //cuts.reserve(corerows);
  VB[controw]=elems;
  for (int i = 0; i < nnodes; ++i)
	for (int k = 0; k < ndemands; ++k){
		if(i==data.d_k[k].O-1){
			cuts.push_back(new BCP_cut_core(-data.d_k[k].quantity,-data.d_k[k].quantity));
			RLB[controw] = -data.d_k[k].quantity;
			RUB[controw] = -data.d_k[k].quantity;
			
			
		}else if(i==data.d_k[k].D-1){
			cuts.push_back(new BCP_cut_core(data.d_k[k].quantity,data.d_k[k].quantity));
			RLB[controw] = data.d_k[k].quantity;
			RUB[controw] = data.d_k[k].quantity;
			
		}else{
			 cuts.push_back(new BCP_cut_core(0,0));
			 RLB[controw] = 0;
			 RUB[controw] = 0;
			
		 }
		
		 for (int a = 0; a < narcs; ++a)
			if(i==data.arcs[a].i-1){
				EV[elems] = -1.0;
				EI[elems] = narcs+a*ndemands+k;
				++elems;
			}else if(i==data.arcs[a].j-1){
				EV[elems] = 1.0;
				EI[elems] = narcs+a*ndemands+k;
				++elems;
			}
		++controw;
		VB[controw]=elems;
	}
  
  //std::cout<<"done "<<vars.size()<<" "<<cuts.size()<<std::endl;
  matrix = new BCP_lp_relax(false, corerows, corecols, elems,
			   VB,  EI, EV,
			   OBJ, CLB, CUB,
			   RLB,RUB);
  
  
 
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
	
	std::cout<<"init phase "<<phase<<std::endl;
	/*if(phase){
		colgen = BCP_DoNotGenerateColumns_Fathom;
		candidates = new CoinSearchTree<CoinSearchTreeCompareBreadth>;
    }else{
		colgen = BCP_DoNotGenerateColumns_Send;
		candidates = new CoinSearchTree<CoinSearchTreeCompareBreadth>;
	}*/
    
	colgen = BCP_DoNotGenerateColumns_Fathom;
    candidates = new CoinSearchTree<CoinSearchTreeCompareDepth>;
}

//#############################################################################
// managing
  
void
MCND_tm::change_candidate_heap(CoinSearchTreeManager& candidates,
			const bool new_solution){
	
	std::cout<<"new_solution "<<new_solution<<std::endl;
	std::cout<<"cand size "<<candidates.size()<<std::endl;
	/*//Best_LB=lower_bound();
	if(!new_solution){
		CoinSearchTreeBase * tree = candidates.getTree();
		CoinSearchTreeBase *t = createSearchType(tree->compName());
		std::cout<<"change candidate heap size:"<<tree->size()<<std::endl;
        MCND_node_branch_data * user_data;
		BCP_tm_node * n;
		CoinTreeNode ** add = new CoinTreeNode*;
		while(tree->size()){
			*add = tree->top();
			n = dynamic_cast<BCP_tm_node*>(*add);
            user_data =  dynamic_cast<MCND_node_branch_data*>(n->_data._user.GetRawPtr());
            (*add)->setQuality(-user_data->score);
			//std::cout<<n->getQuality()<<std::endl;
			if(n->status > BCP_PrunedNode_Discarded ||
			   n->status < BCP_PrunedNode_OverUB){
				t->push(1, add);
			}
			tree->pop();
		}
		candidates.setTree(t);
		delete add;
		std::cout<<"final size "<<t->size()<<std::endl;
	}
	
  std::cout<<"change candidate heap size:"<<std::endl;*/
  BCP_tm_user::change_candidate_heap(candidates,new_solution);
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
	std::cout<<"The global lower bound: "<<Best_LB<<std::endl;
	
}

//-----------------------------------------------------------------------------

void
MCND_tm::display_node_information(BCP_tree& search_tree,
			   const BCP_tm_node& node,
			   bool after_processing_node){
	if(Best_LB<lower_bound())
		Best_LB=lower_bound();
	
	std::cout<<"The lower bound: "<<lower_bound()<<std::endl;
				   
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
