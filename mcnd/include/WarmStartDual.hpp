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
    WarmStartDual(int size, const double* dual, const int* actv);
    WarmStartDual(const CoinWarmStartDual* wsd);
    
    inline ~WarmStartDual(){ if(map) delete [] map;}
    
    //---------------------

    double at(int n) const;
};


#endif /* WarmStartDual_hpp */
