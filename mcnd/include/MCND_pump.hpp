#ifndef _PUMPBCP_H
#define _PUMPBCP_H

#include <ilcplex/ilocplex.h>

#include "MCND_data.hpp"
#include "UtilsMethods.hpp"
#include "MCND_solution.hpp"

#include <vector>
#include <deque>


class Pump{

public:
	int sizeOfIdSeq;
    Pair2 *map;
    
	IloEnv env;
	IloCplex cplex;
	IloModel model;
	
	IloNumVarArray x;
	IloNumVarArray y;
	
	IloObjective fobj;
	IloRangeArray cutstrong;
	 
	double factory;
 	//--------------------------
	
	const Data * data;
	const MCND_solution* best_sol;

	int nnodes, ndemands, narcs;
	int szunfix, maxunfix;
	int tern, tabusz;
 	std::deque<int> unfx;
 	std::vector<int> topo;
 	std::vector<unsigned int *> tabu;
 
	//--------------------------
	Pump(): cplex(env), x(env), y(env), model(env), fobj(env), cutstrong(env) {data =0;  }
	void set_data(const Data * d);
	void initialize(const Data * d, const MCND_solution* best_sol_);
	void set_parameters();
	
	//--------------------------
	int make_topo( const double * x  , const BCP_vec<BCP_var*>& vars);
	int validate_topology( );
	bool check_tabu(unsigned int* seqtopo);
	void try_perturbation(unsigned int* seqtopo);
	//--------------------------

 	void create_model();
	void check_feas_model();
	//--------------------------
	int solve(int *& fixd, int& closed,  const BCP_vec<BCP_var*>& vars,  MCND_solution*& mipsol);
	int cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & y_, const IloNumArray & x_);
	double getSolution(  const IloNumArray & x_, MCND_solution*& mipsol );
	void get_closed(int*& fixd, int& closed, bool onlyx);
	
	~Pump(); //x.endElements(); 
};

#endif 










