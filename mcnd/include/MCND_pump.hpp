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
	 
	double factorxy, factorp;
	double alpha;
	//--------------------------
	
	const Data * data;
	const MCND_solution* best_sol;

	int nnodes, ndemands, narcs;
	int szunfix, maxunfix;
	int tern;
 	std::deque<int> unfx;
 	std::vector<int> topo;
 	std::vector<unsigned int *> tabu;
    std::vector<double> branch_candidates;

	//--------------------------
	Pump(): cplex(env), x(env), y(env), model(env), fobj(env), cutstrong(env) {data =0;  }
	void set_data(const Data * d, double alph_init );
	void initialize(const Data * d, const MCND_solution* best_sol_, double alph_init);
	void set_parameters();
	
	//--------------------------
	int make_topo(int * fixd, int& closed,  const double * x  , const BCP_vec<BCP_var*>& vars);
	bool check_tabu(unsigned int* seqtopo);
	void try_perturbation(unsigned int* seqtopo);
	//--------------------------

	void reset(  const BCP_vec<BCP_var*>& vars);
	void create_model( const BCP_vec<BCP_var*>& vars);
	
	//--------------------------
	int solve( const BCP_vec<BCP_var*>& vars,  double * xy, double & val );
	int cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & y_, const IloNumArray & x_);
	double getSolution(   double * xy,  const IloNumArray & x_, const IloNumArray & y_ );

	
	void clear();
	~Pump(); //x.endElements(); 
};

#endif 










