// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "BCP_buffer.hpp"
#include <iostream>

#include "MCND_data.hpp"
#include <cmath>

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
// DATA
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

BCP_buffer& 
Data::pack(BCP_buffer& buf){
	std::cout<<"pack data"<<std::endl;
    const Arc* item;
    const Demand* itemd;
	buf.pack(nnodes).pack(narcs).pack(ndemands);
	for(int a=narcs;a--;){
        item = &arcs[a];
        buf.pack(item->i).pack(item->j).pack(item->capa).pack(item->f);
        for(int k=ndemands;k--;){
            buf.pack(item->c[k]).pack(item->b[k]);
        }
    }
    for(int k=ndemands;k--;){
        itemd = &d_k[k];
        buf.pack(itemd->O).pack(itemd->D).pack(itemd->quantity);
    }
    for(int i=nnodes*nnodes;i--;){
        buf.pack(grid[i]).pack(kgrid[i]);
    }
	return buf;
}

//-------------------------------------------------------------------------------

BCP_buffer& 
Data::unpack(BCP_buffer& buf){
	std::cout<<"unpack data"<<std::endl;
	buf.unpack(nnodes).unpack(narcs).unpack(ndemands);
    if(arcs.empty()) arcs.resize(narcs);
    if(d_k.empty()) d_k.resize(ndemands);
    if(grid.empty()) grid.resize(nnodes*nnodes);
    if(kgrid.empty()) kgrid.resize(nnodes*nnodes);
        
    Arc* item;
    Demand* itemd;
    for(int a=narcs;a--;){
        item = &arcs[a];
        buf.unpack(item->i).unpack(item->j).unpack(item->capa).unpack(item->f);
        item->c.resize(ndemands);
        item->b.resize(ndemands);
        for(int k=ndemands;k--;){
            buf.unpack(item->c[k]).unpack(item->b[k]);
        }
    }
    for(int k=ndemands;k--;){
        itemd = &d_k[k];
        buf.unpack(itemd->O).unpack(itemd->D).unpack(itemd->quantity);
    }
    for(int i=nnodes*nnodes;i--;){
        buf.unpack(grid[i]).unpack(kgrid[i]);
    }
    //std::cout<<" ok"<<std::endl;

	return buf;
}



	
	
	
	
	
	
	
	
	
	

