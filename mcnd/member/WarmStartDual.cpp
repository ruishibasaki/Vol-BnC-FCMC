//
//  WarmStartDual.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "WarmStartDual.hpp"


WarmStartDual::WarmStartDual(int size, const double* dual, const CoverCollection* covers):CoinWarmStartDual(size, dual){
        int sz = covers->sizeOfCollection;
        Cover *vi = covers->begin;
        for(;sz--;){
        	mapd.insert(std::pair<char,int>(vi->serial_nmbr,vi->id_vi));
        	vi = vi->next;
        }
        dual_ = CoinWarmStartDual::dual();
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int size, const double* dual, const std::map<int, int>& map_ ):CoinWarmStartDual(size, dual){
		dual_ = CoinWarmStartDual::dual();
		mapd = map_;
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int size, const double* dual):CoinWarmStartDual(size, dual){
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(const CoinWarmStartDual* wsd) : CoinWarmStartDual(*wsd){
    	//mapd = wsd->mapd;
        dual_ = CoinWarmStartDual::dual();
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


CoinWarmStartDiff * 
WarmStartDual::generateDiff (const CoinWarmStart *const oldCWS) const{
	std::cout<<"WarmStartDual::generateDiff"<<std::endl;
	return CoinWarmStartDual::generateDiff(oldCWS);

}

//-------------------------------------------------------------------------------------------

void 
WarmStartDual::applyDiff (const CoinWarmStartDiff *const cwsdDiff){
	std::cout<<"WarmStartDual::applyDiff"<<std::endl;
	return CoinWarmStartDual::applyDiff(cwsdDiff);
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

double 
WarmStartDual::get_mapped(int key) const{ // the key is the serial number of the cover
        if(mapd.size()==0 ){//}|| CoinWarmStartDual::size()==0){ 
        	std::cout<<"Attention: map.size()==0 || CoinWarmStartDual::size()==0" <<std::endl;
        	return dual_[key];
        	//return 0;
        }
        std::map<int,int>::const_iterator it = mapd.find(key);
        if(it == mapd.end()) return 0;
        return dual_[it->second];
}

//-------------------------------------------------------------------------------------------

