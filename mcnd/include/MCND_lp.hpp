// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_LP_H
#define _MCND_LP_H

#include "BCP_lp_user.hpp"
#include "BCP_lp_functions.hpp"
#include "BCP_parameters.hpp"
#include "CoinHelperFunctions.hpp"
#include <stdlib.h>     /* abs */

#include "OsiVolSolverInterface.hpp"

#include <vector>
#include <deque>
#include <list>

#include "MCND_cut.hpp"
#include "MCND_data.hpp"
#include "MCND_osidata.hpp"
#include "MCND_solution.hpp"
#include "MCND_branch_score.hpp"
#include "MCND_pump.hpp"
//#include "MCND_topoheur.hpp"
#include "MCND_checklpfeas.hpp"


class OsiVolSolverInterface;
class OsiSolverInterface;

enum LP_Mode{
   /** Branch is driven by a Unfeasible pump problem.*/
   LP_Normal = 0,
   /** Branch is driven by a Feasible pump problem.*/
   LP_Solved = 2,
   /** Extra Iteration due to cut addition.*/
   LP_HeuristicRunned = 4,
   /** Extra Iteration due to cut addition.*/
   LP_CutAddedFromHeuristic = 8,
	/** Force node Abort.*/
   LP_ForceNodeAbort = 16,
   /** logical fixing happend in strong branch.*/
   LP_LogicalFixed = 32,
	/**test connectivity due to Logical fix.*/
   LP_TestConnectivity = 64,
	/** strong branching phase.*/
   LP_StrongBranch = 128,
  /** Make less iterations.*/
   LP_ReducedRun = 256,
   /** Tighter bounds.*/
   LP_tighterBounds = 512,
   /**test connectivity due to Logical fix.*/
   LP_TestFeasibility = 1024
};


class MCND_lp : public BCP_lp_user {
    
public:
    
    Data data;
    OsiVolAuxInfo AppVolData;
    CoverManager cover_manager;
    CutSetManager ss_manager;
    LocalCutManager localc_manager;
    GlobalCutManager globalc_manager;

    LPFeasChecker lpfeaschecker;
    //TopoHeur topo_heur;
    FlowConnect flwconnect;
    Pump pump_heur;
    
    int lp_mode;
    double LBi;
    double lower_bound;
	double nomgap;
	int no_gap_reduct;
	
    bool aborted;
    bool has_sol;
    MCND_solution best_sol;

    std::vector<double> y;
    int* yfix;

    std::vector<Pair> ninsp;
    std::vector<PairF> psdcost;
	std::vector<int> tabu;


	int num_nodes;
    std::deque<MCND_CutUnit *> track;
    BCP_vec<Pair2> to_logical_fix;
    std::map<int, int> mapd;
    
        
public:
    inline MCND_lp() : LBi(0), has_sol(false), track(0), aborted(false) { 
    	lp_mode = LP_Normal;  no_gap_reduct=num_nodes=0; lower_bound=0;
    	nomgap = 1e30;
    }
    ~MCND_lp();
    
    
    //--------------------------------------------------------------------------
    //-------------------------------------------------------------------------- 
    //--------------------------------------------------------------------------
    
    virtual void unpack_module_data(BCP_buffer & buf);

    virtual OsiSolverInterface * initialize_solver_interface();
    void process_message(BCP_buffer& buf);
    
    //--------------------------------------------------------------------------
    //-------------------------------------------------------------------------- 
    //--------------------------------------------------------------------------
    
    
    virtual void load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
                              BCP_var_set& vars, BCP_cut_set& cuts);
    
    virtual void modify_lp_parameters(OsiSolverInterface* lp, const int changeType,
                                      bool in_strong_branching);
    
    
     
    virtual void
    initialize_new_search_tree_node(const BCP_vec<BCP_var*>& vars,
                                    const BCP_vec<BCP_cut*>& cuts,
                                    const BCP_vec<BCP_obj_status>& var_status,
                                    const BCP_vec<BCP_obj_status>& cut_status,
                                    BCP_vec<int>& var_changed_pos,
                                    BCP_vec<double>& var_new_bd,
                                    BCP_vec<int>& cut_changed_pos,
                                    BCP_vec<double>& cut_new_bd);
                                    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------                             
    //--------------------------------------------------------------------------
    // Solution processing
    virtual double
    compute_lower_bound(const double old_lower_bound,
                        const BCP_lp_result& lpres,
                        const BCP_vec<BCP_var*>& vars,
                        const BCP_vec<BCP_cut*>& cuts);
    
    virtual void
    process_lp_result(const BCP_lp_result& lpres,
                      const BCP_vec<BCP_var*>& vars,
                      const BCP_vec<BCP_cut*>& cuts,
                      const double old_lower_bound,
                      double& true_lower_bound,
                      BCP_solution*& sol,
                      BCP_vec<BCP_cut*>& new_cuts,
                      BCP_vec<BCP_row*>& new_rows,
                      BCP_vec<BCP_var*>& new_vars,
                      BCP_vec<BCP_col*>& new_cols);
    
    virtual BCP_solution*
    test_feasibility(const BCP_lp_result& lp_result,
                     const BCP_vec<BCP_var*>& vars,
                     const BCP_vec<BCP_cut*>& cuts);
    
    virtual void
    restore_feasibility(const BCP_lp_result& lpres,
                        const std::vector<double*> dual_rays,
                        const BCP_vec<BCP_var*>& vars,
                        const BCP_vec<BCP_cut*>& cuts,
                        BCP_vec<BCP_var*>& vars_to_add,
                        BCP_vec<BCP_col*>& cols_to_add);
    
    /** Try to generate a heuristic solution (or return one generated during
     cut/variable generation. */
    virtual BCP_solution*
    generate_heuristic_solution(const BCP_lp_result& lpres,
                                const BCP_vec<BCP_var*>& vars,
                                const BCP_vec<BCP_cut*>& cuts);
    
   
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    //managing
    virtual BCP_branching_decision
    select_branching_candidates(const BCP_lp_result& lpres,
                                const BCP_vec<BCP_var*>& vars,
                                const BCP_vec<BCP_cut*>& cuts,
                                const BCP_lp_var_pool& local_var_pool,
                                const BCP_lp_cut_pool& local_cut_pool,
                                BCP_vec<BCP_lp_branching_object*>& candidates,
                                bool force_branch = false);
    
    virtual void
    logical_fixing(const BCP_lp_result& lpres,
                   const BCP_vec<BCP_var*>& vars,
                   const BCP_vec<BCP_cut*>& cuts,
                   const BCP_vec<BCP_obj_status>& var_status,
                   const BCP_vec<BCP_obj_status>& cut_status,
                   const int var_bound_changes_since_logical_fixing,
                   BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd);
    
    
    virtual BCP_branching_object_relation
    compare_branching_candidates(BCP_presolved_lp_brobj* new_presolved,
                                 BCP_presolved_lp_brobj* old_presolved);
    virtual void
    set_actions_for_children(BCP_presolved_lp_brobj* best);
    
    virtual void
    set_user_data_for_children(BCP_presolved_lp_brobj* best,
                               const int selected);
    
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    virtual void
    generate_cuts_in_lp(const BCP_lp_result& lpres,
                        const BCP_vec<BCP_var*>& vars,
                        const BCP_vec<BCP_cut*>& cuts,
                        BCP_vec<BCP_cut*>& new_cuts,
                        BCP_vec<BCP_row*>& new_rows);
    
    virtual void
    cuts_to_rows(const BCP_vec<BCP_var*>& vars, // on what to expand
                 BCP_vec<BCP_cut*>& cuts,       // what to expand
                 BCP_vec<BCP_row*>& rows,       // the expanded rows
                 // things that the user can use for lifting cuts if allowed
                 const BCP_lp_result& lpres,
                 BCP_object_origin origin, bool allow_multiple);//std::cout<<"cuts to rows no"<<std::endl;}
    
    virtual void
    vars_to_cols(const BCP_vec<BCP_cut*>& cuts, // on what to expand
                 BCP_vec<BCP_var*>& vars,       // what to expand
                 BCP_vec<BCP_col*>& cols,       // the expanded cols
                 // things that the user can use for lifting vars if allowed
                 const BCP_lp_result& lpres,
                 BCP_object_origin origin, bool allow_multiple){}//std::cout<<"vars to cols no"<<std::endl;}
    //--------------------------------------------------------------------------
    
    
    virtual void
    select_vars_to_delete(const BCP_lp_result& lpres,
                          const BCP_vec<BCP_var*>& vars,
                          const BCP_vec<BCP_cut*>& cuts,
                          const bool before_fathom,
                          BCP_vec<int>& deletable);
    virtual void
    select_cuts_to_delete(const BCP_lp_result& lpres,
                          const BCP_vec<BCP_var*>& vars,
                          const BCP_vec<BCP_cut*>& cuts,
                          const bool before_fathom,
                          BCP_vec<int>& deletable);
    virtual void
    purge_slack_pool(const BCP_vec<BCP_cut*>& slack_pool,
                     BCP_vec<int>& to_be_purged);
    
    
    
    void
    display_lp_solution(const BCP_lp_result& lpres,
                        const BCP_vec<BCP_var*>& vars,
                        const BCP_vec<BCP_cut*>& cuts,
                        const bool final_lp_solution){};
                        
	//--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // helper functions
    OsiVolSolverInterface* getOsiVolBabSolver();
    virtual void pack_feasible_solution(BCP_buffer& buf, const BCP_solution* sol);
    
    //--------------------------------------------------------------------------

    bool cut_varfix_and_updt(const BCP_vec<BCP_var*>& vars, 
							const BCP_vec<BCP_cut*>& cuts,
							BCP_vec<int>& var_changed_pos,
                            BCP_vec<double>& var_new_bd);
    //--------------------------------------------------------------------------
    
    void update_branch_data(MCND_node_branch_data * ndata);
    void keep_track( MCND_CutUnit* vi);
    void rearrange_bd_vec(int old_sz, int added, BCP_vec< double >& bd_vars );
    
    //--------------------------------------------------------------------------
    int verify_children_feasibility(BCP_lp_branching_object * candidate, bool& feas0, bool& feas1 );
	bool verify_feasibility(const BCP_vec<int> & vars_chngd, const BCP_vec<double> & chngd_bd, int nvars, bool feas_status);
	void strong_branch_var_logicfix(BCP_lp_branching_object * candidate);    
};

#endif

















