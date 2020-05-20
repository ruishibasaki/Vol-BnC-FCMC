//
//  MCND_checklp.hpp
//  
//
//  Created by Rui Shibasaki on 19/05/2020.
//

#ifndef MCND_checklp_hpp
#define MCND_checklp_hpp

#include <stdio.h>
#include <ilcplex/ilocplex.h>

#include "BCP_vector.hpp"

#include "MCND_data.hpp"

#include "covercollection.hpp"
#include "localcutcollection.hpp"
#include "globalcutcollection.hpp"



class LPChecker{

public:
    IloEnv env;
    IloCplex cplex;
    IloModel model;
   
    IloObjective fobj;
 
    IloNumVarArray x;
    IloNumVarArray y;
    IloRangeArray flowconserv;
    
	IloRangeArray add_strong_force;
	IloRangeArray add_covers;
	IloRangeArray add_locals;
	IloRangeArray add_globals;
	
	int szunfx;
	int numaddcov, numaddloc, numaddgloc, numaddstrong;
    std::vector<int> unfixd;
    std::vector<int> map_addcovers;
    std::vector<int> map_addlocals;
    std::vector<int> map_addglobals;

	const CoverCollection* covers;
    const LocalCutCollection* localcs;
    const GlobalCutCollection* globalcs;
	double* dual_sol;

    int ndemands, nnodes, narcs;
    const Data* data;
    
    // construtores
    inline LPChecker():model(env), cplex(env), x(env), y(env), flowconserv(env), 
    					add_strong_force(env), add_covers(env), add_locals(env), add_globals(env) { dual_sol=0;}
     ~LPChecker();
    
    
    void set_parameters();
    void initialize(const Data* d, const CoverCollection* cover_man_,
    				const LocalCutCollection* localc_man_, const GlobalCutCollection* globalc_man_ );
    void reset(const BCP_vec<BCP_var*>& vars);
    bool coverc_viol(const IloNumArray & y_ , const Cover* cover);
	bool localc_viol(const IloNumArray & y_ , const LocalCut* loc);
	bool globalc_viol(const IloNumArray & y_ , const GlobalCut* gloc);

    int cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & y_, const IloNumArray & x_);
    int solve(double ub, const BCP_vec<BCP_var*>& vars,  double& new_lb);
    void get_dual();
 	void clean();
 	
 	void logical_fix(double ub);
 	
    
 };

#endif /* MCND_checklp_hpp */
