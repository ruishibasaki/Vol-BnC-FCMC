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
#include "covercollection.hpp"

class WarmStartDual: public CoinWarmStartDual{
public:
    
    std::map<int, int> mapd; //mapping for extra core cuts;
    const double * dual_; // corevalues/covervalues matching map in sequence.
    
    //---------------------
    inline WarmStartDual(): dual_(0){};
    WarmStartDual(int size, const double* dual, const CoverCollection* covers);
    WarmStartDual(int size, const double* dual, const std::map<int, int>& map_ );
    WarmStartDual(const CoinWarmStartDual* wsd);
    
    inline ~WarmStartDual(){ map.clear;}
    
    //---------------------

    double get_mapped(int key) const;
};


#endif /* WarmStartDual_hpp */
