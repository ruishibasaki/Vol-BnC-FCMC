//
//  WarmStartDual.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "WarmStartDual.hpp"


WarmStartDual::WarmStartDual(int sizei_, const double* primal_, int sized_, const double* dual_, const CoverCollection* covers, 
							const LocalCutCollection* locals, const GlobalCutCollection* globals){
							
		sized = sized_;
 		dual = new double[sized];
 		CoinDisjointCopyN(dual_, sized, dual);
 		sizei = sizei_;
 		primal = new double[sizei];
 		CoinDisjointCopyN(primal_, sizei, primal);
 		
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

WarmStartDual::WarmStartDual(int sizei_, const double* primal_, int sized_, const double* dual_, const std::map<int, int>& map_ ){
 		sized = sized_;
 		dual = new double[sized];
 		CoinDisjointCopyN(dual_, sized, dual);
		mapd = map_;
		sizei = sizei_;
 		primal = new double[sizei];
 		CoinDisjointCopyN(primal_, sizei, primal);
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(int sizei_, const double* primal_, int sized_, const double* dual_){
	sized = sized_;
 	dual = new double[sized];
 	CoinDisjointCopyN(dual_, sized, dual);
 	sizei = sizei_;
 		primal = new double[sizei];
 		CoinDisjointCopyN(primal_, sizei, primal);
}

//-------------------------------------------------------------------------------------------

WarmStartDual::WarmStartDual(const WarmStartDual* wsd){
		sized = wsd->sized;
 		dual = new double[sized];
 		CoinDisjointCopyN(wsd->dual, sized, dual);
    	mapd = wsd->mapd;
    	
    	sizei = wsd->sizei;
 		primal = new double[sizei];
 		CoinDisjointCopyN(wsd->primal, sizei, primal);
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
		//return dual[key];
		return 0;
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
	buf.pack(sized);
	for(int i =sized ; i--;) buf.pack(double(dual[i]));
	
	int sz = mapd.size();
	buf.pack(sz);
	std::map<int,int>::const_iterator it = mapd.begin();
	for( ; it != mapd.end(); ++it) buf.pack(int(it->first)).pack(int(it->second));
	
	buf.pack(sizei);
	for(int i =sizei ; i--;) buf.pack(double(primal[i]));
	
}

//-------------------------------------------------------------------------------------------

void 
WarmStartDual::unpack(BCP_buffer& buf){
	buf.unpack(sized); //std::cout<<"sz dual "<<sz<<std::endl;
	dual = new double[sized];
	for(int i =sized ; i--;) buf.unpack(dual[i]);
	
	int fst, snd, sz;
	buf.unpack(sz); //std::cout<<"sz map "<<sz<<std::endl;
	for( ; sz--;){
		buf.unpack(fst).unpack(snd);
		mapd.insert(std::make_pair(fst, snd));
	}
	
	buf.unpack(sizei); //std::cout<<"sz dual "<<sz<<std::endl;
	primal = new double[sizei];
	for(int i =sizei ; i--;) buf.unpack(primal[i]);

}









