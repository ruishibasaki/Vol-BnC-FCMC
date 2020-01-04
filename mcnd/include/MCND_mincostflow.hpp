#ifndef _FLOW_H
#define _FLOW_H

#include <ilcplex/ilocplex.h>

#include "MCND_data.hpp"
#include "MCND_solution.hpp"
#include <vector>
#include <deque>


class MinCostFLow{

public:
	IloEnv env;
	IloCplex cplex;
	
	IloModel* model;
	// variaveis
	IloNumVarArray x;
	
	
	// construtores
	MinCostFLow();
	

	void set_parameters();
	double getSolution(const std::deque<int> & topo, const Data * data, MCND_solution * sol_);
	int solve();
    void create_model(const std::deque<int> & topo, const Data * data);
	
	~MinCostFLow();

};

#endif 
