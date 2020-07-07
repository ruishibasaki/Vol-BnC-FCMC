//
//  localcutmanager.hpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#ifndef localcutmanager_hpp
#define localcutmanager_hpp

#include <stdio.h>
#include <list>
#include <vector>
#include "localcutcollection.hpp"
#include "MCND_cut.hpp"

class LocalCutManager: virtual public CutManager {

public:
    
    const Data * data;
    const int * arc_map;
    LocalCutCollection locals;
    int num_actv;
    int lim_to_remv;
    int gend;
   
    int ttgend_sellm; //stat can be removed
    int ttgend_feas; //stat can be removed
    int ttgend_opt; //stat can be removed

    std::deque<LocalCut *> purgbl;
    std::vector<int> fixbl_arcs;
    //-------------------------------------------------------------------------------------------
    //  initializing methods
    //-------------------------------------------------------------------------------------------
    
    inline LocalCutManager(): data(0), arc_map(0){ num_actv = lim_to_remv = gend =0;
    												ttgend_opt = ttgend_feas = ttgend_sellm=0; //stat can be removed
    												}
    inline ~LocalCutManager(){fixbl_arcs.clear(); purgbl.clear();}
    inline void set_arc_map(const int * map){ arc_map = map;}
    void initialize(const Data * d, int lim);
    int reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize, bool recheck_collct);
    void clean_collection();
    //-------------------------------------------------------------------------------------------
    //  main methods
    //-------------------------------------------------------------------------------------------
    int localc0_generation_main( const double * ystar,  const int * closed, int sz, int curr_id);
	int make_localcut0(const double * ystar, int sz,  int* vars_, int curr_id);

    //-------------------------------------------------------------------------------------------
    int localc1_generation_main(double lb, double ub, const double * ystar, const int * y, const double * rc, int curr_id, int max);
    int make_localcut1(std::vector<int>& T, const double * ystar, int curr_id, int oneor0);
    int check_fixable(std::list<Pair2>& rc_, std::vector<int>& T, double lb, double ub, int oneor0);
    void form_t(std::list<Pair2>& rc_, std::vector<int>& T, double lb, double ub, double rc_tt);
    
    //-------------------------------------------------------------------------------------------
	int localc2_generation_main( const double * ystar,  const double * topo, int sz, int curr_id);
	int make_localcut2(const double * ystar, int sz,  int* vars_, double* coef_, double rhs_, int curr_id);
    
    //-------------------------------------------------------------------------------------------
    //  auxiliary methods
    //-------------------------------------------------------------------------------------------
    
     void reposition_locals(int added);
    //-------------------------------------------------------------------------------------------
    //  Volume Integration methods
    //-------------------------------------------------------------------------------------------
    
    int compute_localc_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0);
    int compute_localc_sg(const double * x, const int * actvS, int actvSSz,  double * v);
    void add_local_vi(int added, int * actvS, int & actvSSz,  double * dualsol, double *lhsol,
    				 double * h, double * dstar, double * dual_lb, double * dual_ub );

    double arc_dg_imp(int arc, const double * xy, const double * h, const int * actvS, int actvSSz);
};

#endif /* localcutmanager_hpp */
