//
//  cutsetGen.hpp
//  
//
//  Created by Rui Shibasaki on 20/05/2019.
//

#ifndef cutsetmanager_hpp
#define cutsetmanager_hpp

#include <stdio.h>
#include "Structures.hpp"
#include "cutsetcollection.hpp"
#include <deque>


class CutSetManager {
    
public:
    
    const Data * data;
    CutSetCollection sets;
    
    //--------------------------------------
    //  initializing methods
    //---------------------------------------
    inline CutSetManager(const Data * d): data(d){}
    inline CutSetManager(): data(0){}
    void initialize(const Data * d);
    
    //--------------------------------------
    //  main methods
    //---------------------------------------
    int cutset_generation_main( const double * ystar, const double * y, bool prepro);
    int compute_cutset(std::vector<int>& s, const double * ystar, const double * y, int initial, bool stosb, bool prepro);
    
    
    //--------------------------------------
    //  Cutset methods
    //---------------------------------------
        
    void make_cutset( const std::vector<int> & s, std::deque<int> & ss_, double &uss, double &dss, std::deque<int> & s_s, double &us_s, double &ds_s);
    
    int separate_cutset(int i, const double * y,std::vector<int> &s, bool stosb);
    
    int get_jToS(const std::vector<int> & s, const double * y, double &uPlus, double&dssPlus, bool stosb);
    
    //--------------------------------------
    //  Preprocessing methods
    //---------------------------------------
    void make_all_singletons();
    int lin_separate_cutset(int initial, const double * ystar, std::vector<int>& s, bool stosb, double tag);
    void computePowerSet();
    int get_lin_jToS(const std::vector<int> & s, const double * ystar, double &uPlus, bool stosb, double tag);
};

//=================================================================================================
//=================================================================================================
//=================================================================================================


class MinCardCSManager{
public:
    const Data * data;
    
    MinCardCSCollection cardscs;
    
    inline MinCardCSManager(Data * d): data(d){}
    inline MinCardCSManager(): data(0){}
    void initialize(const Data * d);
    
    
    //--------------------------------------
    //  main methods
    //---------------------------------------
    int cardcs_generation_main(const double * ystar, const double * y,const CutSetCollection * sets, int actvSSz );
    int cardcs_generation(int ss_size, const int * SS_arcs, double uss, double dss,
                         const double * ystar, const double * y, int id ,int id_owner_);
    int make_cardcs( std::deque<Pair2> & ss_, double dss );
    
    double cutset_preprocess(int sz, double dss, const int * ss_,  std::deque<Pair2>& ss_deque,
                             const double *y, const double *ystar);
    //--------------------------------------
    //  auxiliary methods
    //---------------------------------------
    bool process_card(const double *y, const double *ystar, MinCardCS * c, double & maxy, double & viol);
    bool checkViol(const MinCardCS * c, const double *y, double &viol);
};


#endif /* cutsetGen_hpp */
