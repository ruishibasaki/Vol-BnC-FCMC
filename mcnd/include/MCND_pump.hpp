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
	IloModel model;
	
	IloNumVarArray x;
	IloNumVarArray y;
	
	IloObjective fobj;
	IloRangeArray cutstrong;
	 
	double factorxy, factorp;
	double alpha;
	//--------------------------
	
	const Data * data;
	int nnodes, ndemands, narcs;
	int szunfix, maxunfix;
 	std::deque<int> unfx;
 	std::vector<int> topo;

	//--------------------------
	Pump(): cplex(env), x(env), y(env), model(env), fobj(env), cutstrong(env) {data =0;  }
	void set_data(const Data * d, double alph_init );
	void initialize(const Data * d, double alph_init);
	void set_parameters();
	
	//--------------------------
	int make_topo(int * fixd, int& closed,  const double * x  , const BCP_vec<BCP_var*>& vars);
	void reset(  const BCP_vec<BCP_var*>& vars);
	void create_model( const BCP_vec<BCP_var*>& vars);
	
	//--------------------------
	int solve( const BCP_vec<BCP_var*>& vars,  double * xy, double & val);
	int cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & y_, const IloNumArray & x_);
	double getSolution(  double * xy,   const IloNumArray & x_, const IloNumArray & y_ );

	
	void clear();
	~Pump(); //x.endElements(); 
};

#endif 










