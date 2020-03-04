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

//----------------------------------------------------------------------------------

MCND_node_branch_data& 
MCND_node_branch_data::operator=(const MCND_node_branch_data& rhs){
     hs = rhs.hs; pos_neg = rhs.pos_neg; branch_var = rhs.branch_var;
    score = rhs.score;
    dual_size = rhs.dual_size; min_lb = rhs.min_lb;
		return *this;
}

//----------------------------------------------------------------------------------

MCND_node_branch_data::MCND_node_branch_data(const MCND_node_branch_data& rhs){
		  hs = rhs.hs; pos_neg = rhs.pos_neg; branch_var = rhs.branch_var;
        score = rhs.score;
            dual_size = rhs.dual_size;
}

//----------------------------------------------------------------------------------

void
MCND_node_branch_data::pack(BCP_buffer& buf) const{
    //std::cout<<"pack MCND_node_branch_data "<<hs<<" size:"<<buf.size()<<std::endl;
    //std::cout<<branch_var<<" "<<pos_neg<<" "<<score<<" "<<score_parent<<std::endl;
    bool hashs = false;
    if(hs && hs->size()>0){
    	hashs = true;
    	buf.pack(hashs);
    	hs->pack(buf);
    }else buf.pack(hashs);
    //std::cout<<branch_var<<" "<<pos_neg<<" "<<score<<" "<<score_parent<<std::endl;
    buf.pack(branch_var).pack(pos_neg).pack(score);
    buf.pack(min_lb);
    //std::cout<<"packed size: "<<buf.size()<<std::endl;
    //std::cout<<branch_var<<" "<<pos_neg<<" "<<score<<" "<<score_parent<<std::endl;
}

//----------------------------------------------------------------------------------

void
MCND_node_branch_data::unpack(BCP_buffer& buf){
    //std::cout<<"unpack MCND_node_branch_data size: "<<buf.size()<<std::endl;
    
    bool hashs;
    buf.unpack(hashs);
    if(hashs){
    	hs = new WarmStartDual();
    	hs->unpack(buf);
    }
    
    buf.unpack(branch_var).unpack(pos_neg).unpack(score);
    buf.unpack(min_lb);
    //std::cout<<branch_var<<" "<<pos_neg<<" "<<score<<" "<<score_parent<<std::endl;
    
}
