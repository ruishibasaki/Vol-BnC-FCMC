#ifndef _FLOW_H
#define _FLOW_H

#include <ilcplex/ilocplex.h>

#include "MCND_data.hpp"
#include <vector>
#include <deque>


class LPFeasChecker{

public:
    IloEnv env;
    IloCplex cplex;
    IloModel model;
    
    IloNumVarArray x;
    

    int ndemands, nnodes, narcs;
    const Data* data;
    
    // construtores
    inline LPFeasChecker():model(env), cplex(env), x(env) {};
    virtual ~LPFeasChecker();
    
    
    void set_parameters();
    void initialize(const Data* d);
    void create_model(const double *collb, const double * colub);
    int solve(const double *collb, const double * colub);

};

#endif 
