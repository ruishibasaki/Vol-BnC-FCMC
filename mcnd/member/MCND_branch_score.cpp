// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.


#include "BCP_USER.hpp"
#include "MCND_branch_score.hpp"

double 
score(double fphin, double fphip, double mu){
	double min, max;
	if(fphin>fphip){
		min = fphip;
		max = fphin;
	}else{
		max = fphip;
		min = fphin;
	}
	return ((1.0 - mu)*min + mu*max);
}

MCND_node_branch_data& 
MCND_node_branch_data::operator=(const MCND_node_branch_data& rhs){
     hs = rhs.hs; pos_neg = rhs.pos_neg; branch_var = rhs.branch_var;
    score = rhs.score; score_parent = rhs.score_parent;
		return *this;
	}


MCND_node_branch_data::MCND_node_branch_data(const MCND_node_branch_data& rhs){
		  hs = rhs.hs; pos_neg = rhs.pos_neg; branch_var = rhs.branch_var;
        score = rhs.score; score_parent = rhs.score_parent;
	 }


void
MCND_node_branch_data::pack(BCP_buffer& buf) const{
   // std::cout<<"pack MCND_node_branch_data "<<hs<<" size:"<<buf.size()<<std::endl;
    
    buf.pack(int(hs->size()));//std::cout<<"hs size: "<<hs->size()<<std::endl;
    const double * d = hs->dual();
    for(int i=hs->size();i--;){ buf.pack(double(d[i]));}//  std::cout<<"hs: "<<d[i]<<std::endl;}
    //delete hs;
    //std::cout<<branch_var<<" "<<pos_neg<<" "<<score<<" "<<score_parent<<std::endl;
    buf.pack(branch_var).pack(pos_neg).pack(score).pack(score_parent).pack(hiters);
    //std::cout<<"packed size: "<<buf.size()<<std::endl;
    
}


void
MCND_node_branch_data::unpack(BCP_buffer& buf){
    //std::cout<<"unpack MCND_node_branch_data size: "<<buf.size()<<std::endl;
    
    int sz;
    buf.unpack(sz);
    double * d = new double [sz]; //std::cout<<"hs size: "<<sz<<std::endl;
    for(int i=sz;i--;){  buf.unpack(d[i]);}//  std::cout<<"hs: "<<d[i]<<std::endl;}
    hs = new CoinWarmStartDual(sz,d);
    //std::cout<<"hs : "<<hs<<std::endl;
    delete [] d;
    
    buf.unpack(branch_var).unpack(pos_neg).unpack(score).unpack(score_parent).unpack(hiters);
   // std::cout<<branch_var<<" "<<pos_neg<<" "<<score<<" "<<score_parent<<std::endl;
    
}
