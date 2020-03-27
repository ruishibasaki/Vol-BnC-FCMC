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

class LocalCutManager {

public:
    
    const Data * data;
    const int * arc_map;
    LocalCutCollection locals;
    int num_actv;
    int lim_to_remv;
    int gend;
    int ttgend;
    std::deque<LocalCut *> purgbl;
    std::vector<int> fixbl_arcs;
    //-------------------------------------------------------------------------------------------
    //  initializing methods
    //-------------------------------------------------------------------------------------------
    
    inline LocalCutManager(): data(0), arc_map(0){ num_actv = lim_to_remv = gend =0;}
    inline ~LocalCutManager(){fixbl_arcs.clear(); purgbl.clear();}
    inline void set_arc_map(const int * map){ arc_map = map;}
    void initialize(const Data * d, int lim);
    int reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize);
    void clean_collection();
    //-------------------------------------------------------------------------------------------
    //  main methods
    //-------------------------------------------------------------------------------------------
    
    int localc_generation_main(double lb, double ub, const double * ystar, const double * y, const double * rc, int curr_id, int max);
    int make_localcut(std::vector<int>& T, const double * ystar, const double * y, int curr_id, int oneor0);
    int check_fixable(std::list<Pair2>& rc_, std::vector<int>& T, double lb, double ub, int oneor0);
    void form_t(std::list<Pair2>& rc_, std::vector<int>& T, double lb, double ub, double rc_tt);
    
    //-------------------------------------------------------------------------------------------
    //  auxiliary methods
    //-------------------------------------------------------------------------------------------
    
    double checkViol(const LocalCut * c, const double *y);
    int add_external_localc(const std::deque<Pair2>& c, int maxNumrows_);
    
    //-------------------------------------------------------------------------------------------
    //  Volume Integration methods
    //-------------------------------------------------------------------------------------------
    
    int compute_cover_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0);
    int compute_cover_sg(const double * x, const int * actvS, int actvSSz,  double * v);
    void add_cover_vi(int added, int * actvS, int & actvSSz,double * h, double * dual, double * dual_lb, double * dual_ub );
      
};

#endif /* localcutmanager_hpp */
