//
//  WarmStartDual.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "WarmStartDual.hpp"


WarmStartDual::WarmStartDual(int size_, const double* dual_, const CoverCollection* covers, 
							const LocalCutCollection* locals, const GlobalCutCollection* globals){
							
		size = size_;
 		dual = new double[size];
 		CoinDisjointCopyN(dual_, size, dual);
 		
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
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int size_, const double* dual_, const std::map<int, int>& map_ ){
 		size = size_;
 		dual = new double[size_];
 		CoinDisjointCopyN(dual_, size, dual);
		mapd = map_;
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int size_, const double* dual_){
	size = size_;
 	dual = new double[size_];
 	CoinDisjointCopyN(dual_, size, dual);
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(const WarmStartDual* wsd){
		size = wsd->size;
 		dual = new double[size];
 		CoinDisjointCopyN(wsd->dual, size, dual);
    	mapd = wsd->mapd;
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
		return dual[key];
		//return 0;
	}
	//std::cout<<"ok "<<size()<<std::endl;
	std::map<int,int>::const_iterator it = mapd.find(key);
	if(it == mapd.end()) return 0.0;
	//std::cout<<key<<" "<<it->second<<std::endl;
	//std::cout<<dual_[it->second]<<std::endl;
	//std::cout<<key<<" "<<it->second<<std::endl;
	return dual[it->second];
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void 
WarmStartDual::pack(BCP_buffer& buf) const{
	buf.pack(size);
	for(int i =size ; i--;) buf.pack(double(dual[i]));
	
	int sz = mapd.size();
	buf.pack(sz);
	std::map<int,int>::const_iterator it = mapd.begin();
	for( ; it != mapd.end(); ++it) buf.pack(int(it->first)).pack(int(it->second));

	
}

//-------------------------------------------------------------------------------------------

void 
WarmStartDual::unpack(BCP_buffer& buf){
	buf.unpack(size); //std::cout<<"sz dual "<<sz<<std::endl;
	dual = new double[size];
	for(int i =size ; i--;) buf.unpack(dual[i]);
	
	int fst, snd, sz;
	buf.unpack(sz); //std::cout<<"sz map "<<sz<<std::endl;
	for( ; sz--;){
		buf.unpack(fst).unpack(snd);
		mapd.insert(std::make_pair(fst, snd));
	}

}









