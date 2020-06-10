//
//  MCND_checklp.hpp
//  
//
//  Created by Rui Shibasaki on 19/05/2020.
//

#ifndef MCND_checklp_hpp
#define MCND_checklp_hpp

#include <stdio.h>
#include <map>
 
#include <ilcplex/ilocplex.h>

#include "BCP_vector.hpp"

#include "MCND_data.hpp"
#include "MCND_solution.hpp"

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
    
    IloRangeArray topo_ctrnt;
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

    //-------------------------------
    std::map<int, double> dual_map;  
	double solval;
	bool solved;
	std::vector<double> bound_red;
	IloNumArray x_;
	IloNumArray y_;
	IloNumArray rc_y;
	IloNumArray rc_x;
    //----------------------------------

    int ndemands, nnodes, narcs;
    const Data* data;
    const CoverCollection* covers;
    const LocalCutCollection* localcs;
    const GlobalCutCollection* globalcs;
    
    //----------------------------------
	// Methods
    //----------------------------------

    inline LPChecker():model(env), cplex(env), x(env), y(env), flowconserv(env), 
    					add_strong_force(env), add_covers(env), add_locals(env), add_globals(env),
    					x_(env), y_(env), rc_x(env), rc_y(env), topo_ctrnt(env) {numaddstrong = numaddcov= numaddloc= numaddgloc =0; solved=false; }
     ~LPChecker();
    
    //-------------------------------
    //-------------------------------

    void set_parameters();
    void initialize(const Data* d, const CoverCollection* cover_man_,
    				const LocalCutCollection* localc_man_, const GlobalCutCollection* globalc_man_ );
    void reset(const BCP_vec<BCP_var*>& vars, const double *collb, const double *colub);
   
    //-------------------------------
    //-------------------------------

	 bool coverc_viol(const IloNumArray & y_ , const Cover* cover);
	bool localc_viol(const IloNumArray & y_ , const LocalCut* loc);
	bool globalc_viol(const IloNumArray & y_ , const GlobalCut* gloc);
    int cut(const BCP_vec<BCP_var*>& vars, const double *colub);
    int cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & yb, const IloNumArray & xb, IloRangeArray& tempcnstrt, int& addtemp);
    //-------------------------------

    int solve(double ub, const BCP_vec<BCP_var*>& vars,  const double *collb, const double *colub, double& new_lb);
    int test_high_lowerbounds(  BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, int * yfix,
										 const BCP_vec<BCP_var*>& vars,  const double* volsol, double bestsolv);
    void get_dual();
    void get_solution(MCND_solution* mipsol);
 	void clean();
 	
 	//-------------------------------
    //-------------------------------

 	bool logical_yfix(double ub, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, int * yfix);
 	void logical_xfix(double ub, const int * yfix, const BCP_vec<BCP_var*>& vars);
    
 };

#endif /* MCND_checklp_hpp */
