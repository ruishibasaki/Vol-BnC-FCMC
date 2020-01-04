#ifndef _FLOW_H
#define _FLOW_H

#include <ilcplex/ilocplex.h>

#include "MCND_data.hpp"
#include "MCND_solution.hpp"
#include <vector>


class MinCostFeas{

public:
	IloEnv* env;
	IloCplex* cplex;
	
	IloModel* model;
	// variaveis
	IloNumVarArray x;
    std::vector<int> idx;
	
	// construtores
	MinCostFeas(const Data * data);
    
    int nnodes;
    int ndemands;
    int narcs;

	

	void set_parameters();
	int solve();
	void create_model(const Topology & topo, const Data * data);
    void clear_model();
	~MinCostFeas();

};

#endif 
