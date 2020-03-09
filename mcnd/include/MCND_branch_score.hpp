// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_branch_H
#define _MCND_branch_H

#include <vector>
#include "BCP_buffer.hpp"
#include "WarmStartDual.hpp"
#include "OsiAuxInfo.hpp"
#include "BCP_USER.hpp"

double 
score(double fphin, double fphip, double mu);



class MCND_node_branch_data: public BCP_user_data{
private:
	MCND_node_branch_data& operator=(const MCND_node_branch_data& rhs);
	
public: 

    WarmStartDual * hs;
    int dual_size;
    double score;
    double min_lb;
	int pos_neg;
	int branch_var;
 
	
	MCND_node_branch_data(const MCND_node_branch_data& rhs);
	
    MCND_node_branch_data(): hs(0), min_lb(0) {}
	
    MCND_node_branch_data(int dual_size_, double score_, int branch_var_, int zrone,   double min_lb_):hs(0){
         score = score_; branch_var = branch_var_; pos_neg = zrone;  
        dual_size = dual_size_;
        min_lb =min_lb_;
    }
    
	inline ~MCND_node_branch_data(){
        //std::cout<<"delete MCND_node_branch_data "<<hs<<std::endl;
         if(hs!=0)delete hs;
        //std::cout<<" MCND_node_branch_data deleted"<<std::endl;

    }
	
    void pack(BCP_buffer& buf) const;
    void unpack(BCP_buffer& buf);


};

#endif
