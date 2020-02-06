//
//  WarmStartDual.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "WarmStartDual.hpp"


WarmStartDual::WarmStartDual(int size, const double* dual, const int* actv):CoinWarmStartDual(size, dual){
        if(size>0) map = new int [size];
        for(int i=size; i--;) map[i] = actv[i];
        dual_ = CoinWarmStartDual::dual();
}

//-------------------------------------------------------------------------------------------


WarmStartDual::WarmStartDual(const CoinWarmStartDual* wsd) : CoinWarmStartDual(*wsd){
    	int sz =  wsd->size();
        if(sz>0) map = new int [sz];
        for(int i=0; i<sz;++i) map[i] = i;
        dual_ = CoinWarmStartDual::dual();
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

double 
WarmStartDual::at(int n) const{
        if(map==0 || CoinWarmStartDual::size()==0) return 0;
        if(map[n]<0) return 0;
        return dual_[map[n]];
}

//-------------------------------------------------------------------------------------------

