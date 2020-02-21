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
	
	double factorxy, factorp;
	double alpha;
	//--------------------------
	
	const Data * data;
	int nnodes, ndemands, narcs;
	int sznz, szunfix;
	std::deque<Pair2> ytopo;

	//--------------------------
	Pump(): cplex(env), x(env), y(env){data =0; model =0;}
	void set_data(const Data * d, double alph_init );
	void set_parameters();
	
	//--------------------------
	void reset( int szunfix_, std::deque<Pair2>& topo );
	void create_model();
	
	//--------------------------
	int solve(double * yl, double * xy, double & val);
	int cut(const IloNumArray & y_, const IloNumArray & x_);
	double getSolution(double * yl, double * xy, const IloNumArray & x_, const IloNumArray & y_ );

	
	void clear();
	~Pump(); //x.endElements(); 
};

#endif 










