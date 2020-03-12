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

bool 
CoverCut::check_logical_fix(const BCP_vec<BCP_var*>& vars, double* yarcs){
	double sum=0;
    double rhs= cover->get_rhs();
    int sz = cover->get_total_sz();
    int dimsh=0;
    int arc;
    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	arc = cover->at(a);
    	if(vars[arc]->lb() > 0.5){
        	dimsh+= cover->gamma_at(a);
		}else if(vars[arc]->ub() > 0.5){
			sum+= cover->gamma_at(a);
		}
    }
	rhs -= dimsh;
	if(sum <= rhs && rhs>0){
		for(int a=0;a<sz;++a){
			arc = cover->at(a);
			if(vars[arc]->lb() < 0.5 && vars[arc]->ub() > 0.5){
				 yarcs[arc] = 1;
			}
		}
		std::cout<<"serial: "<<cover->serial_nmbr<<" "; cover->print();
		return true;
	}
	return false;
}




