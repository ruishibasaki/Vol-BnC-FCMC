// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_TM_H
#define _MCND_TM_H

#include "BCP_tm_user.hpp"
#include "BCP_parameters.hpp"

#include "MCND_tm_param.hpp"

#include "MCND_data.hpp"
#include "MCND_solution.hpp"

class MCND_tm : public BCP_tm_user {
public:
 
  Data data;
  double Best_LB;
  std::string instance;
  
  MCND_solution best_soln;
  
  clock_t t_start;
public:
  MCND_tm() : Best_LB(0) {}
  ~MCND_tm() {}

  
   
  //--------------------------------------------------------------------------
  // pack the module data for the appropriate process
  void
  pack_module_data(BCP_buffer& buf, BCP_process_t ptype);
  //--------------------------------------------------------------------------
  // unpack an MIP feasible solution
  BCP_solution*
  unpack_feasible_solution(BCP_buffer& buf);
  
  virtual void
  process_message(BCP_buffer& buf);
  //--------------------------------------------------------------------------
  // feasible solution displaying
  void
  display_solution(const BCP_solution* soln);
  //--------------------------------------------------------------------------
  // setting the core
  void
  initialize_core(BCP_vec<BCP_var_core*>& vars,
		  BCP_vec<BCP_cut_core*>& cuts,
		  BCP_lp_relax*& matrix);
  //--------------------------------------------------------------------------
  // create the root node
  void
  create_root(BCP_vec<BCP_var*>& added_vars,
	      BCP_vec<BCP_cut*>& added_cuts,
	      BCP_user_data*& user_data);
	
  void
  init_new_phase(int phase,
		 BCP_column_generation& colgen,
		 CoinSearchTreeBase*& candidates);

  //-------------------------------------------------------------------------
  // managing
  /*CoinSearchTreeBase *
  createSearchType(const char * type);*/
  
  void
  change_candidate_heap(CoinSearchTreeManager& candidates,
			const bool new_solution);
  //--------------------------------------------------------------------------
  // Display a feasible solution
  void
  display_feasible_solution(const BCP_solution* sol);
  
  virtual void
  display_final_information(const BCP_lp_statistics& lp_stat);
  
  virtual void
  display_node_information(BCP_tree& search_tree,
			   const BCP_tm_node& node,
			   bool after_processing_node);
  //##########################################################################
};

#endif
