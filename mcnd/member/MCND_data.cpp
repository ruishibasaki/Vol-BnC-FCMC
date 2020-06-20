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
    redone=false;
    int szlbl = narcs*ndemands;
    bound_red.resize(szlbl);
    labelf.resize(szlbl,0);
    labelb.resize(szlbl,0);
    Kd.resize(nnodes);
    Ko.resize(nnodes);
    adjf = new std::list<int>[nnodes];
    adjb = new std::list<int>[nnodes];
    
    for(int i=nnodes;i--;){
        Kd[i] = new int[ndemands+1];
        Ko[i] = new int[ndemands+1];
        Kd[i][0]=0;
        Ko[i][0]=0;
    }
    int orig, dest, pos;
    int* setm; int* setp;
    for(int i=0;i<ndemands;++i){
        orig = data->d_k[i].O-1;
        dest = data->d_k[i].D-1;
        setm = Kd[dest];
        setp = Ko[orig];
        pos = ++setm[0];
        setm[pos] = i;
        pos = ++setp[0];
        setp[pos] = i;
        
    }
}

//-------------------------------------------------------------------------------

void
FlowConnect::set_bd(int id, double lb, double ub, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd){
	changed_pos.push_back(id);
	new_bd.push_back(lb);
	new_bd.push_back(ub);
	if(id>=narcs){
		bound_red[id-narcs].fst = changed_pos.size()-1;
		bound_red[id-narcs].snd = lb;
 		bound_red[id-narcs].trd = ub;
	}
}

//-------------------------------------------------------------------------------

int
FlowConnect::check_flowbounds(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, int * yfix){
    std::list<int>::const_iterator it;
    int arc, id, idvar;
 	double sum_ubin, sum_lbin, ubnew, lbnew;
	double sum_ubout, sum_lbout, lb;
	bool modifub=false;
	bool modiflb=false;

    for(int n=nnodes;n--;){
        //std::cout<<"node "<<n+1<<std::endl;
        for(int k = ndemands; k--;){
        	if(n == (data->d_k[k].O-1) || n == (data->d_k[k].D-1)) continue;
        	modifub=false; modiflb=false;
        	sum_ubin=0; sum_ubout=0; sum_lbin=0; sum_lbout=0;
        	for (it = adjb[n].begin(); it != adjb[n].end(); ++it){
				arc = data->grid[(*it)*nnodes+n];
				idvar = k*narcs+arc; 
			 	ubnew = bound_red[idvar].trd;
			 	lbnew = bound_red[idvar].snd;
				sum_ubin += (ubnew >= 0)? fmin(ubnew,  vars[narcs+idvar]->ub()) : vars[narcs+idvar]->ub();
				sum_lbin +=  fmax(lbnew,  vars[narcs+idvar]->lb());
            }
        	
            for (it = adjf[n].begin(); it != adjf[n].end(); ++it){
                arc = data->grid[n*nnodes+(*it)];      
                idvar = k*narcs+arc;          
				ubnew = bound_red[k*narcs+arc].trd;
 				lb = fmax(vars[narcs+idvar]->lb(), bound_red[idvar].snd);
 				sum_lbout += lb; 
				if( lb > (sum_ubin+ 1e-4)) return -1;
                
				ubnew = (ubnew >= 0)? fmin(ubnew, vars[narcs+idvar]->ub()) : vars[narcs+idvar]->ub();
				if(sum_ubin < ubnew-1e-4){
					//std::cout<<"ATENCAO!! "<<sum<<" "<<ubnew<<std::endl;
					id = bound_red[idvar].fst;
					if(id>=0){
 						new_bd[id*2+1] = sum_ubin;
 						bound_red[idvar].trd = sum_ubin;
					}else{
						set_bd(narcs+idvar, lb, sum_ubin,  changed_pos, new_bd);
					}
					ubnew = sum_ubin;
					modifub=true;
				}
 				sum_ubout += ubnew;
            }
            if(sum_ubout < sum_lbin) return -1;
            else if(sum_ubin < sum_lbout) return -1;
            else if(sum_lbin>0){
             	for (it = adjf[n].begin(); it != adjf[n].end(); ++it){
                	arc = data->grid[n*nnodes+(*it)];    
            		idvar = k*narcs+arc;      
            		lb = fmax(vars[narcs+idvar]->lb(), bound_red[idvar].snd);
					ubnew = bound_red[idvar].trd;                
					ubnew = (ubnew >= 0)? fmin(ubnew, vars[narcs+idvar]->ub()) : vars[narcs+idvar]->ub();
					lbnew = 0;
					if((sum_ubout-ubnew) < sum_lbin){
						lbnew = sum_lbin-(sum_ubout-ubnew);
						if(lbnew > lb){
							id = bound_red[idvar].fst;
							if(id>=0){
								new_bd[id*2] = lbnew;
								bound_red[idvar].snd = lbnew;
							}else{
								set_bd(narcs+idvar, lbnew, ubnew,  changed_pos, new_bd);
							}
							if(yfix[arc]==-1){
								set_bd(arc, 1.0, 1.0,  changed_pos, new_bd);
								yfix[arc]=1;
								//std::cout<<"flowbounds fix "<<arc<<" to 1"<<std::endl;
							}
							modiflb=true;
						} 
					}
             	}
            }
            if(!modifub){
            	sum_ubin=0;
				for (it = adjb[n].begin(); it != adjb[n].end(); ++it){
					arc = data->grid[(*it)*nnodes+n];
					idvar = k*narcs+arc;
					lb = fmax(vars[narcs+idvar]->lb(), bound_red[idvar].snd);				
					ubnew = bound_red[idvar].trd;
					ubnew = (ubnew >= 0)? fmin(ubnew, vars[narcs+idvar]->ub()) : vars[narcs+idvar]->ub();
					if(sum_ubout < ubnew-1e-4){
						//std::cout<<"ATENCAO!! "<<sumb<<" "<<ubnew<<std::endl;
						id = bound_red[idvar].fst;
						if(id>=0){
							new_bd[id*2+1] = sum_ubout;
							bound_red[idvar].trd = sum_ubout;
						}else{
							set_bd(narcs+idvar, lb, sum_ubout,  changed_pos, new_bd);
						}
						ubnew = sum_ubout;
					}
					sum_ubin += ubnew ;
				 }
             }
             if(!modiflb){
 				 for (it = adjb[n].begin(); it != adjb[n].end(); ++it){
					 arc = data->grid[(*it)*nnodes+n];    
					 idvar = k*narcs+arc;      
					 lb = fmax(vars[narcs+idvar]->lb(), bound_red[idvar].snd);
					 ubnew = bound_red[idvar].trd;                
					 ubnew = (ubnew >= 0)? fmin(ubnew, vars[narcs+idvar]->ub()) : vars[narcs+idvar]->ub();
					 lbnew = 0;
					 if((sum_ubin-ubnew) < sum_lbout){
						 lbnew = sum_lbout-(sum_ubin-ubnew);
						 if(lbnew > lb){
							 id = bound_red[idvar].fst;
							 if(id>=0){
								 new_bd[id*2] = lbnew;
								 bound_red[idvar].snd = lbnew;
							 }else{
								 set_bd(narcs+idvar, lbnew, ubnew,  changed_pos, new_bd);
							 }
							 if(yfix[arc]==-1){
								 set_bd(arc, 1.0, 1.0,  changed_pos, new_bd);
								 yfix[arc]=1;
								 //std::cout<<"flowbounds fix "<<arc<<" to 1"<<std::endl;
							 }
 						 } 					 
					 }
             	 }
             }
        }   
    }
    return 1;
}

//-------------------------------------------------------------------------------

int
FlowConnect::translate_results(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd,bool& conn_arcfix, int * yfxd){
    std::list<int>::const_iterator it;
    int arc, id;
    int arc_unique;
    int comm_sat = 0;
	int ret =1;
    std::vector<int>count_closed(narcs,0);
    	
    for(int n=nnodes;n--;){
        //std::cout<<"node "<<n+1<<std::endl;
        for(int k = ndemands; k--;){
            arc_unique = -1;
            comm_sat = 0;
            for (it = adjf[n].begin(); it != adjf[n].end(); ++it){
                arc = data->grid[n*nnodes+(*it)];
                id = narcs+k*narcs+arc;
                if(labelf[k*narcs+arc] && labelb[k*narcs+arc] && vars[id]->ub()>0){
                    arc_unique = arc;
                    ++comm_sat;
                }else{
                    if(vars[id]->ub()>0 && vars[id]->lb()==0){
                        //std::cout<<"set x_"<<arc+1<<"^"<<k+1<<" to zero. "<<vars[narcs+k*narcs+arc]->lb()<<" "<<vars[narcs+k*narcs+arc]->ub()<<std::endl;
                    	set_bd(id, 0.0, 0.0,  changed_pos, new_bd);
                    	redone=true;
                    }
                    else if(vars[id]->lb()>0){ /*std::cout<<"FlowConnect conflict1"<<std::endl; */ret= -2;}
                    
                    if(++count_closed[arc] == ndemands){
                    	if(yfxd[arc]==-1 && vars[arc]->ub()>0.5){
                    		set_bd(arc, 0.0, 0.0,  changed_pos, new_bd);
                    		conn_arcfix=true;
                    		yfxd[arc]=0;
							//std::cout<<"connfix "<<arc<<" to 0"<<std::endl;
                    	}else if(vars[arc]->lb()>0 || yfxd[arc]==1){ 
                    		/*std::cout<<"FlowConnect conflict2 "<<arc<<" "<<vars[arc]->lb()<<std::endl;*/ 
                    		ret= -2;}
                    }
                }
            }
            if(n == (data->d_k[k].O-1)){
                if(comm_sat==1){
                	id = narcs+k*narcs+arc_unique;
                	double dk = data->d_k[k].quantity;
                	//std::cout<<vars[id]->ub()<<" "<<dk<<std::endl;
                	if(vars[id]->ub() < dk){ /*std::cout<<"FlowConnect conflict3"<<std::endl;*/ ret= -2;}
                	else if(vars[id]->lb() < dk){
                		set_bd(id, dk, dk,  changed_pos, new_bd);
                 		redone=true;
                	} 
                	if(vars[arc_unique]->lb() < 0.5 && yfxd[arc_unique]==-1){
                		set_bd(arc_unique, 1.0, 1.0,  changed_pos, new_bd);
                    	//std::cout<<"set x_"<<arc_unique+1<<"^"<<k+1<<" to the flow capacity"<<std::endl;
                    	conn_arcfix=true;
                    	yfxd[arc_unique]=1;
						//std::cout<<"connfix "<<arc_unique<<" to 1"<<std::endl;
                	} 
                }else if(comm_sat==0){
                	//std::cout<<"FlowConnect conflict4"<<std::endl;  
                    return -1;
                }
            }else if(n == (data->d_k[k].D-1)){
                comm_sat = 0;
                for (it = adjb[n].begin(); it != adjb[n].end(); ++it){
                    arc = data->grid[(*it)*nnodes+n];
                    if(labelf[k*narcs+arc] && labelb[k*narcs+arc]  && vars[narcs+k*narcs+arc]->ub()>0){
                        ++comm_sat;
                        arc_unique = arc;
                        if(comm_sat>1) break;
                    } 
                }
                if(comm_sat==1){
                	int id = narcs+k*narcs+arc_unique;
                	double dk = data->d_k[k].quantity;
                	//std::cout<<vars[id]->ub()<<" "<<dk<<std::endl;
                	if(vars[id]->ub() < dk){ /*std::cout<<"FlowConnect conflict5"<<std::endl;*/ ret= -2;}
                	else if(vars[id]->lb() < dk){
                		set_bd(id, dk, dk,  changed_pos, new_bd);
                 		redone=true;
                	}
                	if(vars[arc_unique]->lb() < 0.5 && yfxd[arc_unique]==-1){
                		set_bd(arc_unique, 1.0, 1.0,  changed_pos, new_bd);
                    	//std::cout<<"set x_"<<arc_unique+1<<"^"<<k+1<<" to the flow capacity"<<std::endl;
                    	conn_arcfix=true;
                    	yfxd[arc_unique]=1;
						//std::cout<<"connfix "<<arc_unique<<" to 1"<<std::endl;
                	} 
                }else if(comm_sat==0){
                    //std::cout<<"FlowConnect conflict6"<<std::endl;  
                    return -1;
                }
            }
        }
        
    }
    return ret;
}

//-------------------------------------------------------------------------------

int 
FlowConnect::check_connectivity(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, bool& conn_arcfix, int * yfxd){
    
    reset(vars, yfxd);
    
    for(int i=nnodes;i--;){
        if(Ko[i][0]>0){
            BFS(true, i, adjf, labelf, Ko[i]);
        }
        if(Kd[i][0]>0){
            BFS(false, i, adjb, labelb, Kd[i]);
        }
    }
    int ret =  translate_results(vars, changed_pos, new_bd, conn_arcfix,yfxd);
    return ret;
}

//-------------------------------------------------------------------------------

void
FlowConnect::BFS(bool forwback, int s, const std::list<int> * adj, std::vector<int>& label,const int* K ){ 
    
     bool *visited = new bool[nnodes];
    for(int i = 0; i < nnodes; i++)
        visited[i] = false;
    
    std::list<int> queue;
    visited[s] = true;
    queue.push_back(s);
    
    int arc, comm, nk;
    std::list<int>::const_iterator i;
    while(!queue.empty()){
        
        s = queue.front();
        queue.pop_front();
        
        for (i = adj[s].begin(); i != adj[s].end(); ++i){
            if(forwback)arc = data->grid[s*nnodes+(*i)];
            else arc = data->grid[(*i)*nnodes+s];
            //std::cout<<"node i: "<<s+1<<" adj"<<forwback<<": "<<*i+1<<" arc "<<arc<<std::endl;

            nk = K[0];
            for(int k=1;k<=nk;++k){
                comm = K[k];
                label[comm*narcs+arc] = 1;
            }
            
            if (!visited[*i]) {
                visited[*i] = true;
                queue.push_back(*i);
            }
        }
    }
    delete [] visited;
} 

//-------------------------------------------------------------------------------

void 
FlowConnect::reset(const BCP_vec<BCP_var*>& vars, int * yfxd){
	for(int i=nnodes;i--;){
		adjf[i].clear();
		adjb[i].clear();
    }
    int szlbl = narcs*ndemands;
    labelf.assign(szlbl,0);
    labelb.assign(szlbl,0);
    
    for(int arc=narcs; arc--;){
    	if(vars[arc]->ub()<0.5){
    		yfxd[arc] =0;
    		continue;
    	}else if( yfxd[arc] ==0){
    		continue;
    	}else if(vars[arc]->lb()>0.5) yfxd[arc]=1;
    	//std::cout<<"arc: "<<arc<<std::endl;
    	 
        const Arc& a = data->arcs[arc];
        adjf[a.i-1].push_back(a.j-1);
        adjb[a.j-1].push_back(a.i-1);
    }
}

//-------------------------------------------------------------------------------

void 
FlowConnect::clear_memory(){
	//std::cout<<"FlowConnect::clear_memory"<<std::endl;
	redone=false;
 	bound_red.assign(narcs*ndemands,Trio1(-1,-1.0,-1.0));	
}

//-------------------------------------------------------------------------------

FlowConnect::~FlowConnect(){
	
	labelf.clear();
	labelb.clear();
	for(int i=data->nnodes;i--;){
		delete [] Ko[i] ;
		delete [] Kd[i] ;	
	}
	delete []adjf;
	delete []adjb;
	bound_red.clear();
 
	Ko.clear();
	Kd.clear();	
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
// DATA
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

BCP_buffer& 
Data::pack(BCP_buffer& buf){
	//std::cout<<"pack data"<<std::endl;
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
	//std::cout<<"unpack data"<<std::endl;
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



	
	
	
	
	
	
	
	
	
	

