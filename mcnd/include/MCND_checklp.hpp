#ifndef _MCND_CKECK_H
#define _MCND_CKECK_H

#include <ilcplex/ilocplex.h>

#include "MCND_data.hpp"
#include "MCND_solution.hpp"
#include "covercollection.hpp"
#include <vector>


class LPChecker{

public:
    
    IloEnv env;
    IloCplex * cplex;
    IloModel * model;
    
    
    int ndemands, nnodes, narcs;
    const Data* data;
    
    // construtores
    inline LPChecker(){};
    virtual ~LPChecker();
    
    
    void set_parameters();
    void create_model(const std::vector<int> & topo, const CoverCollection & covers);
    void initialize(const Data* d);
    // resolve o problema
    int solve();

};

#endif 
