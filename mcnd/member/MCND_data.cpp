// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "BCP_buffer.hpp"
#include <iostream>

#include "MCND_data.hpp"
#include <cmath>

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
// FlowConnect
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

void 
FlowConnect::initialize(Data *d){
	data = d;
	narcs =  data->narcs;
	ndemands = data->ndemands;
	nnodes  = data->nnodes;
	int szlbl = narcs*ndemands;
	labelf.resize(szlbl,0);
	labelb.resize(szlbl,0);
	Km.resize(nnodes,0);
	Kp.resize(nnodes,0);
	adjf = new list<int>[nnodes]; 
	adjb = new list<int>[nnodes]; 

	std::deque<int> comm;
	for(int i=nnodes;i--;){
		Km[i] = new int[ndemands+1];
		Kp[i] = new int[ndemands+1];
		Km[i][0]=0;
		Kp[i][0]=0;
	}
	int orig, dest, pos;
	int* setm; int* setp;
	for(int i=ndemands;i--;){
		orig = data->d_k[i].O-1;
		dest = data->d_k[i].D-1;
		setm = Km[dest];
		setp = Kp[orig];
		pos = ++setm[0];
		setm[pos] = i;
		pos = ++setp[0];
		setp[pos] = i;
		
	}
}

//-------------------------------------------------------------------------------

bool 
FlowConnect::check_connectivity(const deque<int>& nonzro){
	int arc, i, j;
	const Arc& a;
	for(int i=nonzro.size();i--;){
		arc = nonzro[i];
		a = data->arcs[arc];
		adjf[a.i-1].push_back(a.j-1);
		adjb[a.j-1].push_back(a.i-1);
	}
	for(int i=nnodes;i--;){
		if(Kp[0]>0){
			BFS(true, i, adjf, labelf, Kp[i]);
		}
		if(Km[0]>0){
			BFS(false, i, adjb, labelb, Km[i]);
		}
	}
}

//-------------------------------------------------------------------------------

void
FlowConnect::BFS(bool forwback, int s, const list<int> * adj, std::vector<int>& label,const int* K ){ 
    
    bool *visited = new bool[nnodes]; 
    for(int i = 0; i < nnodes; i++) 
        visited[i] = false; 
  
    std::list<int> queue; 
    visited[s] = true; 
    queue.push_back(s); 
  
    int arc, comm, nk;
    std::list<int>::iterator i; 
    while(!queue.empty()){ 

        s = queue.front(); 
        queue.pop_front(); 
        
        for (i = adj[s].begin(); i != adj[s].end(); ++i){
        	if(forwback)arc = data->grid[s*nnodes+(*i)];
        	else arc = data->grid[(*i)*nnodes+s];
        	nk = K[0];
        	for(int k=1;k<=nk;++k){
        		comm = K[k];
        		label[comm*ndemands+arc] = 1;
        	}
        	
            if (!visited[*i]) { 
                visited[*i] = true; 
                queue.push_back(*i); 
            } 
        } 
    } 

//-------------------------------------------------------------------------------

FlowConnect::~FlowConnect(){
	
	labelf.clear();
	labelb.clear();
	for(int i=data->ndemands;i--;){
		delete [] Km[i] ;
		delete [] Kp[i] ;	
	}
	delete []adjf;
	delete []adjb;

	Km.clear();
	Kp.clear();	
}

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



	
	
	
	
	
	
	
	
	
	

