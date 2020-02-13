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
	IloModel * model;
	
	IloNumVarArray x;
	IloNumVarArray y;
	
	double factorx, factory;

	//--------------------------
	
	const Data * data;
	int nnodes, ndemands, narcs;
	int sznz, szunfix;
	double alpha, v0;
	std::deque<Pair2> ytopo;

	//--------------------------

	Pump(const Data * d);
	void set_parameters();
	
	//--------------------------
	void reset(double alph, int szunfix_, std::deque<Pair2>& topo );
	void create_model();
	
	//--------------------------
	int solve(double * yl, double * xy, double & val);
	int cut(const IloNumArray & y_, const IloNumArray & x_);
	double getSolution(double * yl, double * xy, const IloNumArray & x_, const IloNumArray & y_ );

	
	void clear();
	~Pump(); //x.endElements(); 
};

#endif 










