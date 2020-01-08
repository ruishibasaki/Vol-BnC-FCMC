//
//  WarmStartDual.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//

#ifndef WarmStartDual_hpp
#define WarmStartDual_hpp

#include <stdio.h>
#include "CoinWarmStartDual.hpp"

class WarmStartDual: public CoinWarmStartDual{
public:
    
    int * map;
    const double * dual_;
    
    //---------------------
    inline WarmStartDual():map(0), dual_(0){};
    inline WarmStartDual(int size, const double* dual, const int* actv):CoinWarmStartDual(size, dual){
        if(size>0) map = new int [size];
        for(int i=size; i--;) map[i] = actv[i];
        dual_ = CoinWarmStartDual::dual();
    }
    inline ~WarmStartDual(){ if(map) delete [] map;}
    
    inline double at(int n) const{
        if(map==0 || CoinWarmStartDual::size()==0) return 0;
        if(map[n]<0) return 0;
        return dual_[map[n]];
    }
};


#endif /* WarmStartDual_hpp */
