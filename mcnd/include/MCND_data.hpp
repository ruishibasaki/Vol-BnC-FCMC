// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_H
#define _MCND_H

#include <iostream>
#include <vector>
#include <list>
#include "BCP_var.hpp"
#include "BCP_cut.hpp"
#include "BCP_vector.hpp"
#include "Structures.hpp"

class BCP_buffer;


class Stats{
public:
    int num_strongb_fix;
    int num_ls_fix_solve;
    int num_ls_fix;
	int num_rclbflow_fix;
	int num_rcubflow_fix;
    int num_rcarc_fix;
	int num_conn_fix;
    int num_cpflow_fix;
    int num_cparc_fix; 
    int num_solve_holm;
    int num_holm_fathom;
    
  	int num_ttgend_global; 
    int num_ttgend_cover;  
    int num_ttgend_card;  
    int num_ttgend_sellm; 
    int num_ttgend_feas; 
    int num_ttgend_opt; 
    
    bool strongb_fix;
    bool ls_fix_solve;
    bool ls_fix;
    bool rclbflow_fix;
    bool rcubflow_fix;
    bool rcarc_fix;
	bool conn_fix;
    bool cpflow_fix;
    bool cparc_fix; 
    bool solve_holm;
    bool holm_fathom;
    
    Stats();
    void reset();
    void compute_stats();
    void compute_cutstat(int ttgend_global, 
    int ttgend_cover,  
    int ttgend_card,  
    int ttgend_sellm, 
    int ttgend_feas, 
    int ttgend_opt );
};

//-----------------------------------------------------
// DATA
//-----------------------------------------------------

class Data{
public:
    std::vector<Arc> arcs;
    std::vector<Demand> d_k;
    std::vector<int> grid;
    std::vector<double> kgrid;
    int ndemands, narcs, nnodes; // number of demands, number of arcs, number of nodes
    inline ~Data(){
        d_k.clear();
        arcs.clear();
        grid.clear();
        kgrid.clear();
    }
    
    BCP_buffer& pack(BCP_buffer& buf);
    BCP_buffer& unpack(BCP_buffer& buf);
};

//-----------------------------------------------------

class FlowConnect{
public:
    Data * data;
    int nnodes, ndemands, narcs;
    std::vector<int*> Kd;
    std::vector<int*> Ko;
	std::vector<int> labelf;
	std::vector<int> labelb;
    
	std::list<int> *adjf;
	std::list<int> *adjb;

	std::vector<Trio1> bound_red;
 	bool redone;

	inline FlowConnect(){}
	inline FlowConnect(Data* d){ initialize(d);}
    ~FlowConnect();
	void initialize(Data *d);
	
	void set_bd(int id, double lb, double ub, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd);
	void BFS(bool forwback, int s, const std::list<int> * adj, std::vector<int>& label,const int* K );

	int check_connectivity(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, bool& conn_arcfix, int * yfxd);
    int translate_results(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd,bool& conn_arcfix, int * yfxd);
    int check_flowbounds(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd,  int * yfix);
    void reset(const BCP_vec<BCP_var*>& vars, int * yfxd);
    void clear_memory();
};




#endif
