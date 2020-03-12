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

	inline FlowConnect(){}
	inline FlowConnect(Data* d){ initialize(d);}
    ~FlowConnect();
	void initialize(Data *d);
	
    void BFS(bool forwback, int s, const std::list<int> * adj, std::vector<int>& label,const int* K );
	bool check_connectivity(const std::deque<int>& nonzro, BCP_vec<Pair2>& collb, BCP_vec<Pair2>& colub);
    bool translate_results(const std::deque<int>& nonzro, BCP_vec<Pair2>& collb, BCP_vec<Pair2>& colub);
};




#endif
