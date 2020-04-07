//
//  globalcutmanager.hpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#ifndef globalcutmanager_hpp
#define globalcutmanager_hpp

#include <stdio.h>
#include <list>
#include <vector>
#include "globalcutcollection.hpp"
#include "MCND_cut.hpp"

class GlobalCutManager: virtual public CutManager {

public:
    
    const Data * data;
    const int * arc_map;
    GlobalCutCollection globals;
    int num_actv;
    int lim_to_remv;
 	

     //-------------------------------------------------------------------------------------------
    //  initializing methods
    //-------------------------------------------------------------------------------------------
    
    inline GlobalCutManager(): data(0), arc_map(0){ num_actv = lim_to_remv =0;}
    inline ~GlobalCutManager(){}
    inline void set_arc_map(const int * map){ arc_map = map;}
    void initialize(const Data * d, int lim);
    int reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize, bool recheck_collct);
    void clean_collection();
    //-------------------------------------------------------------------------------------------
    //  main methods
    //-------------------------------------------------------------------------------------------
    
    int globalc_generation_main( const double * ystar, const double * ycoef, int cont0, int curr_id, int type);
    int make_globalcut(const double * ystar, int sz,  int* vars_,  double * coef_, double rhs, int curr_id, int type);
     
    //-------------------------------------------------------------------------------------------
    //  auxiliary methods
    //-------------------------------------------------------------------------------------------   
      void reposition_globals(int added);
      void collect_globals(const BCP_vec<BCP_var*>& vbd, int & currnum);
    //-------------------------------------------------------------------------------------------
    //  Volume Integration methods
    //-------------------------------------------------------------------------------------------
    
    int compute_cover_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0);
    int compute_cover_sg(const double * x, const int * actvS, int actvSSz,  double * v);
    void add_global_vi(int added, int * actvS, int & actvSSz,double * h, double * dual, double * dual_lb, double * dual_ub );
    double arc_dg_imp(int arc, const double * xy, const double * h, const int * actvS, int actvSSz);
    
    //-------------------------------------------------------------------------------------------
    //  manager methods
    //-------------------------------------------------------------------------------------------
    void collect_global(double* collb, double* colub, int & curr_numrows);
};

#endif /* localcutmanager_hpp */
