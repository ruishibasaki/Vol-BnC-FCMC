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
#include "localcutcollection.hpp"

#include "BCP_buffer.hpp"


class WarmStartDual: public CoinWarmStartDual{
public:
    
    std::map<int, int> mapd; //mapping for extra core cuts;
    const double * dual_; // corevalues/covervalues matching map in sequence.
    
    //---------------------
    inline WarmStartDual(): dual_(0){};
    WarmStartDual(int size, const double* dual, const CoverCollection* covers, const LocalCutCollection* locals);
    WarmStartDual(int size, const double* dual, const std::map<int, int>& map_ );
    WarmStartDual(int size, const double* dual);
    WarmStartDual(const WarmStartDual* wsd);
    
    inline ~WarmStartDual(){ mapd.clear();}
    CoinWarmStart * clone () const{ return new WarmStartDual(this);}
    WarmStartDual * clone_ws () const{ return new WarmStartDual(this);}

	//---------------------

     CoinWarmStartDiff * 	generateDiff (const CoinWarmStart *const oldCWS) const;
     void applyDiff (const CoinWarmStartDiff *const cwsdDiff);

    //---------------------

    double get_mapped(int key) const;
    
    //---------------------
    
    void pack(BCP_buffer& buf) const;
    void unpack(BCP_buffer& buf);
};


#endif /* WarmStartDual_hpp */
