//
//  MCND_cut.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "MCND_cut.hpp"

#include "covercollection.hpp"
#include "localcutcollection.hpp"

int CutManager::ttgend = 0;


bool
CoverCut::check_viol(const BCP_vec<BCP_var*>& vars){
	double sum=0;
    double comp= cover->get_rhs();
    int sz = cover->get_total_sz();
    cover->rhs_dimsh=0;
    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	if(vars[cover->at(a)]->lb() > 0.5)
        	sum+= cover->gamma_at(a);
        if(sum>=comp){return false;}
    }
    cover->rhs_dimsh = sum;
    return true;
}

//-------------------------------------------------------------------------------------------

double 
CoverCut::check_viol(const double* vars){
	return cover->viol(vars);
}

//-------------------------------------------------------------------------------------------
// 2 fix open
// 3 fix zero
bool 
CoverCut::check_logical_fix(const BCP_vec<BCP_var*>& vars, int* yarcs){
	double sum=0;
    double rhs= cover->get_rhs();
    int sz = cover->get_total_sz();
    int dimsh=0;
    int arc;
    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	arc = cover->at(a);
		if(/*vars[arc]->lb() > 0.5 ||*/ yarcs[arc]==1 || yarcs[arc]==2){
        	dimsh+= cover->gamma_at(a);
		}else if(/*vars[arc]->ub() > 0.5 &&*/ yarcs[arc]==-1){
			sum+= cover->gamma_at(a);
		}
    }
	rhs -= dimsh;
	if(sum == rhs && rhs>0){
		for(int a=0;a<sz;++a){
			arc = cover->at(a);
			if(vars[arc]->lb() < 0.5 && vars[arc]->ub() > 0.5){
				 yarcs[arc] = 2;
			}
		}
		std::cout<<"serial: "<<cover->serial_nmbr<<" "; cover->print();
		return true;
	}else if(sum < rhs){std::cout<<"check_logical_fix::PRUNE: "<<sum<<" < "<<rhs<<" "; cover->print(); return false;}
	return true;
}

//======================================================================
//======================================================================
// LocalCCut methods
//======================================================================
//======================================================================

bool
LocalCCut::check_viol(const BCP_vec<BCP_var*>& vars){
	double sum=0;
    double rhs= localc->rhs;
    int sz = localc->size;
    localc->rhs_dimsh=0;
    if(localc->sense==1){
		for(int a=0;a<sz;++a){
			//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
			if(vars[localc->vars[a]]->lb() > 0.5)
				sum+= 1;
			if((rhs-sum) <= 0){ std::cout<<"LocalCCutout id: "<<localc->id_vi<<" srial: "<<localc->serial_nmbr
								<<" : "<<rhs<<" - "<<sum<<std::endl;return false;}
		}
		localc->rhs_dimsh = sum;
    }else{
    	double sumzro =0;
    	for(int a=0;a<sz;++a){
			//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
			if(vars[localc->vars[a]]->lb() > 0.5)
				sum+= 1;
			else if(vars[localc->vars[a]]->ub() < 0.5)
				sumzro+=1;
			if(sum >= rhs || sumzro>=(sz)){ std::cout<<"LocalCCutout id: "<<localc->id_vi<<" srial: "<<localc->serial_nmbr
								<<" : "<<rhs<<" - "<<sum<<", "<<sumzro<<"/"<<(sz)<<std::endl;return false;}
		}
		localc->rhs_dimsh = sum;
    }
    return true;
}

//-------------------------------------------------------------------------------------------

double 
LocalCCut::check_viol(const double* vars){
	return localc->viol(vars);
}

//-------------------------------------------------------------------------------------------
// 2 fix open
// 3 fix zero
bool 
LocalCCut::check_logical_fix(const BCP_vec<BCP_var*>& vars, int* yarcs){
	double sum=0;
    double rhs= localc->rhs;
    int sz = localc->size;
    int dimsh=0;
    int arc;
    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	arc = localc->vars[a];
    	if(/*vars[arc]->lb() > 0.5 ||*/ yarcs[arc]==1 || yarcs[arc]==2){
        	dimsh+= 1.0;
		}else if(/*vars[arc]->ub() > 0.5 &&*/ yarcs[arc]==-1){
			sum+= 1.0;
		}
    }
	rhs -= dimsh;
	if(localc->sense==1){
		if(sum == rhs && rhs>0){
			for(int a=0;a<sz;++a){
				arc = localc->vars[a];
				if(vars[arc]->lb() < 0.5 && vars[arc]->ub() > 0.5){
					 yarcs[arc] = 2;
				}
			}
			std::cout<<"serial: "<<localc->serial_nmbr<<" "; localc->print();
			return true;
		}else if(sum < rhs){std::cout<<"check_logical_fix::PRUNE: "<<sum<<" < "<<rhs<<" "; localc->print(); return false;}
	}else{
		if(rhs==0){
			for(int a=0;a<sz;++a){
				arc = localc->vars[a];
				if(vars[arc]->lb() < 0.5 && vars[arc]->ub() > 0.5){
					 yarcs[arc] = 3;
				}
			}
			std::cout<<"serial: "<<localc->serial_nmbr<<" "; localc->print();
			return true;
		}else if(rhs < 0){std::cout<<"check_logical_fix::PRUNE: "<<sum<<" < "<<rhs<<" "; localc->print(); return false;}
	}
	return true;
}



