#ifndef _PUMPBCP_H
#define _PUMPBCP_H

#include <ilcplex/ilocplex.h>
#include "VolVolume.hpp"

#include "MCND_data.hpp"
#include "UtilsMethods.hpp"
#include "MCND_solution.hpp"

#include <vector>
#include <deque>


class Pump : public VOL_user_hooks{

public:

	IloEnv env;
	IloCplex cplex;
	IloModel model;
	
	IloNumVarArray x;
	IloNumVarArray y;
	
	IloObjective fobj;
	IloRangeArray cutstrong;
	 
	double factory;
	double roundup;
	double roundwn;
	 
	//--------------------------
	int tern, tabusz;
	int sizeOfIdSeq;
    Pair2 *map;
 	std::vector<unsigned int *> tabu;

 	//--------------------------
	double altsolval;
	double * altsol;
	const Data * data;
	const MCND_solution* best_sol;

	int nnodes, ndemands, narcs;
	int szunfix, maxunfix;
  	std::vector<int> unfx;
 	std::vector<int> topo;
    const BCP_vec<BCP_var*>* varsp;
	//--------------------------
	Pump(): cplex(env), x(env), y(env), model(env), fobj(env), cutstrong(env), volsolver("volmcnd.par") {data =0; altsolval=0; altsol=0; }
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
	
	//====================================================================
	//====================================================================
	//  volume hooks
	//====================================================================
	//====================================================================
	std::vector<int> fxone;
	int szfxone;
	int sznz;
	VOL_problem volsolver;
	
	int volsolve();
	double knapsack(int a, bool none, const double * rc, double* x);
	int compute_rc(const VOL_dvector& u, VOL_dvector& rc, int actvSSz);

  
	int solve_subproblem(const VOL_dvector& xstar,
                const VOL_dvector& dual, VOL_dvector& rc,
				double& lcost, VOL_dvector& x, double& pcost);
    
  	int resolve_subproblem(const VOL_dvector& dual, VOL_dvector& rc,
                           double& lcost, VOL_dvector& x,double& pcost){return 0;};
    
    int additional_settings(int iter, double& lcost, VOL_dvector& dual, VOL_dvector& rc, VOL_dvector& h,
                                    VOL_dvector& x, const VOL_dvector& xhist,  int actvSSz){return 0;};
   
    int heuristics(const VOL_problem& p,
			  const VOL_dvector& x, double& heur_val){return 0;};
    
    int compute_sg(const VOL_dvector& x, int  actvSSz, VOL_dvector& v);
    
    int addVI(int iter,double lcost, const VOL_dvector& xstar,
          const VOL_dvector& x, VOL_dvector& dstar,  VOL_dvector& dualu, VOL_dvector& dual_lb, VOL_dvector& dual_ub,
          VOL_dvector& rc, VOL_dvector& h, int & actvSSz){return 0;};
    
    int removeVI( int & actvSSz, VOL_dvector& pstarv, VOL_dvector& dstaru,  VOL_dvector& dualu){return 0;};
};

#endif 










