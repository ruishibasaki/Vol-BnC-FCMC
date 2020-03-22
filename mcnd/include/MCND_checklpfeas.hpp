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
   
    IloObjective fobj1;
    IloObjective fobj2;

    IloNumVarArray x;
    

    int ndemands, nnodes, narcs;
    const Data* data;
    
    // construtores
    inline LPFeasChecker():model(env), cplex(env), x(env), fobj1(env), fobj2(env) {};
    virtual ~LPFeasChecker();
    
    
    void set_parameters();
    void initialize(const Data* d);
    int solve_opt(const BCP_vec<BCP_var*>& vars, const double * topo);
    int solve_feas(const double *collb, const double * colub);
    int solve();
    int getSolution(const BCP_vec<BCP_var*>& vars, double * sol, double & solvalue, double &fathmval);

};

#endif 
