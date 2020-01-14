//
//  MCND_cut.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "MCND_cut.hpp"


bool
CoverCut::check_viol(const BCP_vec<BCP_var*>& vars){
	double sum=0;
    double comp= cover->get_total_rhs();
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
