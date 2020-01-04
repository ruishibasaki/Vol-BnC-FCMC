// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_H
#define _MCND_H

#include <iostream>
#include <vector>
#include <deque>
#include "BCP_var.hpp"
#include "BCP_cut.hpp"
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




#endif
