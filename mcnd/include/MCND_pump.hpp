#ifndef _PUMPBCP_H
#define _PUMPBCP_H

#include <ilcplex/ilocplex.h>

#include "MCND_data.hpp"
#include <vector>
#include <deque>


class Pump{

public:
	IloEnv env;
	IloCplex cplex;
	
	IloModel* model;
	
	IloObjective Obj;
	
	// variaveis
	IloNumVarArray x;
	IloNumVarArray y;
	
	std::vector<int> idy;
	std::vector<int> idx;
	int nx, ny;
	
	
	//contraints
	IloRangeArray flow_row;
	IloRangeArray capa_row;
	IloRangeArray strong_row; 
	
	std::vector<int> idf_row;
	std::vector<int> idc_row;
	std::vector<int> ids_row;
	int nfr, ncr , nsr;
	
	std::vector<int> yint;
	double factorx, factory;
	
	//problem data
	const Data * data;
	const Topology * topo;
	int nnodes;
	int ndemands;
	int narcs;
	int sznz;
	bool feas_init;
	
	
	
	// construtores
	
	Pump(const Topology * t, const Data * d);
	
	
	
	void set_parameters();
	double getSolution(double * xy);
	int solve();
	void create_model(int numcols, const int * init_cols, double alpha, bool feas);
	
	//price and cut methods
	int price(std::deque<int>& cols_to_add,const IloNumArray & y_,const IloNumArray & pi_,const IloNumArray & alpha_ );
	void cut(const IloNumArray & y_, const IloNumArray & x_);
	void add_cols(std::deque<int>& cols_to_add);
	~Pump();
};

#endif 

