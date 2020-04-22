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
#include "globalcutcollection.hpp"

#include "BCP_buffer.hpp"


class WarmStartDual: public CoinWarmStartDual{
public:
    
    std::map<int, int> mapd; //mapping for extra core cuts;
    double * dual; // corevalues/covervalues matching map in sequence.
    int size;
    
    //---------------------
    inline WarmStartDual(): dual(0){size=0;};
    WarmStartDual(int size, const double* dual, const CoverCollection* covers, const LocalCutCollection* locals, const GlobalCutCollection* globals);
    WarmStartDual(int size, const double* dual, const std::map<int, int>& map_ );
    WarmStartDual(int size, const double* dual);
    WarmStartDual(const WarmStartDual* wsd);
    
    inline ~WarmStartDual(){ mapd.clear(); if(size) delete [] dual;}
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
