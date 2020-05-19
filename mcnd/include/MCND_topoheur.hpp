#ifndef _MCND_CKECK_H
#define _MCND_CKECK_H

#include <ilcplex/ilocplex.h>

#include "covermanager.hpp"
#include "localcutmanager.hpp"
#include "globalcutmanager.hpp"

#include "MCND_data.hpp"
#include "MCND_solution.hpp"
#include "UtilsMethods.hpp"
#include <vector>
#include <list>
#include <algorithm>


class TopoHeur{

public:
    IloEnv env;
    IloCplex cplex;
    IloModel model;
   
    IloObjective fobj;
 
    IloNumVarArray x;
    IloNumVarArray xextra;

	std::vector<int> unfx;
 	std::vector<int> topo;
 	std::vector<unsigned int *> tabu;
 	
 	const CoverCollection* covers;
    const LocalCutCollection* localcs;
    const GlobalCutCollection* globalcs;
 	
	int tern, tabusz;
	int szunfixd;
    int ndemands, nnodes, narcs;
    int sizeOfIdSeq;
    Pair2 *map;

    const Data* data;
    const MCND_solution* best_sol;
    // construtores
	//--------------------------
    inline TopoHeur():model(env), cplex(env), x(env), fobj(env), xextra(env){ data =0; best_sol=0; };
    virtual ~TopoHeur();
    
    //--------------------------

    void set_parameters();
    void initialize(const Data* d, const MCND_solution* best_sol_, const CoverCollection* cover_man_,
    				const LocalCutCollection* localc_man_, const GlobalCutCollection* globalc_man_);
    void make_model();
    void remake_model();
	//--------------------------
	bool check_tabu(unsigned int* seqtopo);
	void try_perturbation(unsigned int* seqtopo);
	
	//--------------------------
	int solve(int& closed, int*& fixd0, MCND_solution*& mipsol);
    int make_topo(  const double * x , const BCP_vec<BCP_var*>& vars);
	int getSolution(MCND_solution*& mipsol);
	//--------------------------

	int check_feas(bool restore);
	int final_feas( std::list<Pair2>& heap,  const IloNumArray & xextrasol);
};

#endif 
