//
//  WarmStartDual.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "WarmStartDual.hpp"


WarmStartDual::WarmStartDual(int size, const double* dual, const CoverCollection* covers, const LocalCutCollection* locals, const GlobalCutCollection* globals):CoinWarmStartDual(size, dual){
        int sz = covers->sizeOfCollection;
        Cover *vi = covers->begin;
        for(;sz--;){
        	mapd.insert(std::pair<int,int>(vi->serial_nmbr,vi->id_vi));
        	vi = vi->next;
        }
        sz = locals->sizeOfCollection;
        LocalCut *vilc = locals->begin;
        for(;sz--;){
        	mapd.insert(std::pair<int,int>(vilc->serial_nmbr,vilc->id_vi));
        	vilc = vilc->next;
        }
        sz = globals->sizeOfCollection;
        GlobalCut *vigc = globals->begin;
        for(;sz--;){
        	mapd.insert(std::pair<int,int>(vigc->serial_nmbr,vigc->id_vi));
        	vigc = vigc->next;
        }
        dual_ = CoinWarmStartDual::dual();
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int size, const double* dual, const std::map<int, int>& map_ ):CoinWarmStartDual(size, dual){
		dual_ = CoinWarmStartDual::dual();
		mapd = map_;
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int size, const double* dual):CoinWarmStartDual(size, dual){}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(const WarmStartDual* wsd) : CoinWarmStartDual(*wsd){
    	mapd = wsd->mapd;
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
	//std::cout<<"ok "<<size()<<std::endl;
	std::map<int,int>::const_iterator it = mapd.find(key);
	if(it == mapd.end()) return 0.0;
	//std::cout<<key<<" "<<it->second<<std::endl;
	//std::cout<<dual_[it->second]<<std::endl;
	//std::cout<<key<<" "<<it->second<<std::endl;
	return dual_[it->second];
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void 
WarmStartDual::pack(BCP_buffer& buf) const{
	int sz = size(); //std::cout<<"sz dual "<<sz<<" "<< mapd.size()<<std::endl;
	buf.pack(sz);
	for(int i =sz ; i--;) buf.pack(double(dual_[i]));
	
	sz = mapd.size();
	buf.pack(sz);
	std::map<int,int>::const_iterator it = mapd.begin();
	for( ; it != mapd.end(); ++it) buf.pack(int(it->first)).pack(int(it->second));

	
}

//-------------------------------------------------------------------------------------------

void 
WarmStartDual::unpack(BCP_buffer& buf){
	int sz;
	buf.unpack(sz); //std::cout<<"sz dual "<<sz<<std::endl;
	double* d = new double[sz];
	for(int i =sz ; i--;) buf.unpack(d[i]);
	assignDual(sz,d);
	dual_ = CoinWarmStartDual::dual();

	int fst, snd;
	buf.unpack(sz); //std::cout<<"sz map "<<sz<<std::endl;
	for( ; sz--;){
		buf.unpack(fst).unpack(snd);
		mapd.insert(std::make_pair(fst, snd));
	}

}









