#include "MCND_checklp.hpp"



#define INF 9999999
ILOSTLBEGIN

//---------------------------------------------------------------------------

void 
LPChecker::set_parameters() {
	cplex.setParam(IloCplex::Threads,1);
    //cplex->setParam(IloCplex::RootAlg, 2);
    //cplex->setParam(IloCplex::NodeAlg, 3);
    cplex.setParam(IloCplex::ClockType, 1);
    //cplex->setParam(IloCplex::MIPDisplay, 4);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::DataCheck, CPX_DATACHECK_ASSIST);
    cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds
    
    int sizeOfInt=8*sizeof(unsigned int);
    sizeOfIdSeq = (narcs/sizeOfInt)+1;
    
	map = new Pair2[narcs];
    for(int i=0;i<narcs;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
    }
    
}

//---------------------------------------------------------------------------

void
LPChecker::initialize(const Data* d, const MCND_solution* best_sol_, const CoverCollection* cover_man_,
    				const LocalCutCollection* localc_man_, const GlobalCutCollection* globalc_man_ ){
    data = d;
    best_sol = best_sol_;
    covers = cover_man_;
    localcs = localc_man_;
    globalcs = globalc_man_;
    
    ndemands = d->ndemands;
    nnodes = d->nnodes;
    narcs =d->narcs;
    
    unfx.resize(narcs);
    topo.resize(narcs);
	tabu.resize(1000,0);
	tabusz=0;
	tern=0;
    
    set_parameters();
    
    IloExpr obj(env);
    for(int a=0;a<narcs;++a){
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(env, 0.0, IloInfinity));
             obj += data->arcs[a].c[k]*x[a*ndemands+k];
        }
    }
    for (int k = 0; k < ndemands; ++k){
		xextra.add(IloNumVar(env, 0.0, data->d_k[k].quantity));
		obj += 1e8*xextra[k];
	}
	
    fobj = IloMinimize(env, obj);
    model.add(fobj);
    obj.end();

    for(int a=0;a<narcs;++a){
        IloExpr constraint(env);
        for (int k = 0; k < ndemands; k++) {
            constraint += x[a*ndemands+k];
        }
        constraint -= data->arcs[a].capa;
        model.add(constraint <= 0);
        constraint.end();
    }
    for (int k = 0; k < ndemands; ++k) {
        for (int i = 1; i <= nnodes; i++) {
            IloExpr constraint(env);
            
            for(int a=0;a<narcs;++a){
                if(i == data->arcs[a].i){
                    constraint -= x[a*ndemands+k];
                }else if(i == data->arcs[a].j){
                    constraint += x[a*ndemands+k];
                }
            }
            
            if( i == data->d_k[k].D){
            	constraint += xextra[k];
                constraint -= data->d_k[k].quantity;
            }else if( i ==  data->d_k[k].O){
            	constraint -= xextra[k];
                constraint += data->d_k[k].quantity;
            }
            model.add(constraint == 0);
            constraint.end();
        }
    }
    
    cplex.extract(model);
}

//-------------------------------------------------------------------------------------------

void
LPChecker::make_model(){
	double ub;
	for(int a=0;a<narcs;++a){
		if(topo[a]<=0.5){
			ub = 0.0; 
 		}else{ ub = IloInfinity;}
	
		for(int k=0;k<ndemands;++k){
			x[a*ndemands+k].setUB(ub);
		}
	}
}


//-------------------------------------------------------------------------------------------

void
LPChecker::remake_model(){
	double ub;
	int arc;
	for(int a=0;a<szunfixd;++a){
    	arc = unfx[a];
		if(topo[arc]<=0.5){
			ub = 0.0; 
 		}else{ ub = IloInfinity;}
	
		for(int k=0;k<ndemands;++k){
			x[arc*ndemands+k].setUB(ub);
		}
	}
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
LPChecker::make_topo(  const double * x , const BCP_vec<BCP_var*>& vars){
  	unsigned int* seqtopo = new unsigned int[sizeOfIdSeq];
	std::fill(seqtopo, seqtopo+sizeOfIdSeq, 0);
	Pair2* maptitem=0;
	szunfixd=0;
	int rnd;
	int val;
	for(int a=narcs; a--;){
		maptitem = &map[a];
        if(vars[a]->lb()==1.0){
            topo[a]=1;
            setBit(seqtopo, maptitem->fst, maptitem->snd);
        }else if(vars[a]->ub()==1.0){
			rnd = rand()%2;
			if(rnd==1){
				if(best_sol->xy[a]>0.5){
					topo[a]=1;
					setBit(seqtopo, maptitem->fst, maptitem->snd);
				}else{
					topo[a]=0;
				}
			}else{
				rnd = ((rand()%100)/100.0 <= x[a]) ? 1 : 0;
				if(rnd){
					topo[a]=1;
					setBit(seqtopo, maptitem->fst, maptitem->snd);
				}else{
					topo[a]=0;
				}
			}		 
			unfx[szunfixd++] = a;         
        }else{ 
         	topo[a]=0;
        }
    }
   if(check_tabu(seqtopo)){
    	try_perturbation(seqtopo);
    	if(check_tabu(seqtopo)){
    		delete [] seqtopo;
     		return -1;
     	}
    }
	
    return szunfixd;
}

//-------------------------------------------------------------------------------------------

bool
LPChecker::check_tabu(unsigned int* seqtopo){
	unsigned int* titem;
	bool equal;

	for(int i=tabusz;i--;){
		titem = tabu[i];
		equal=true;
		for(int id=0;id<sizeOfIdSeq;++id){
			if(titem[id] != seqtopo[id]){
				equal=false;
				break;
			}
		}
		if(equal){
			return true;
		}
    }  

    if(tern==tabu.size()){ tern=0;}
    unsigned int*& tabuitem = tabu[tern];
	if(tabuitem) { delete [] tabuitem;}
	tabuitem = seqtopo;
	if(tabusz<tabu.size())++tabusz;
	++tern;

    return false;
}

//-------------------------------------------------------------------------------------------

void
LPChecker::try_perturbation(unsigned int* seqtopo){
	int arc;
	Pair2* maptitem=0;
	std::cout<<"LPChecker::try_perturbation"<<std::endl;
	std::random_shuffle(unfx.begin(), unfx.begin()+szunfixd);
 	for(int a=0;a<szunfixd;++a){
		arc = unfx[a];
		maptitem = &map[arc];
		if(topo[arc]==1){
			topo[arc]=0;
			clearBit(seqtopo, maptitem->fst, maptitem->snd);
		}else if(topo[arc]==0){
			topo[arc]=1;
			setBit(seqtopo, maptitem->fst, maptitem->snd);
		}
		if(rand()%2)++a;
	}

}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int
LPChecker::solve(int& closed, int*& fixd0, MCND_solution*& mipsol){
	
	make_model();
	cplex.solve();
	
	int ret;
	if(cplex.getStatus() == IloAlgorithm::Optimal){
        ret = check_feas(true);
        if(ret<0){
        	fixd0 = new int[narcs];
			closed=0;
			for(int a=0;a<narcs;++a){
				if(topo[a]==0){
					fixd0[closed++] = a;
					//std::cout<<"chekclosed: "<<a<<std::endl;
				}
			}
        	if(ret==-1){
				remake_model();
				cplex.solve();
				if(cplex.getStatus() == IloAlgorithm::Optimal){
					ret = check_feas(false);
					if(ret<0) return -2;
				   
				    //std::cout<<"LPChecker: "<<cplex.getObjValue()<<std::endl;
        			getSolution(mipsol);
        			return 0;
        		}else return -3;
			}else if(ret==-2) return -2;
        }
        //std::cout<<"LPChecker: "<<cplex.getObjValue()<<std::endl;
        getSolution(mipsol);
		return 1; 
	}
	if(cplex.getStatus() == IloAlgorithm::Infeasible){	 
		fixd0 = new int[narcs];
		closed=0;
		for(int a=0;a<narcs;++a){
			if(topo[a]==0){
				fixd0[closed++] = a;
				//std::cout<<"chekclosed: "<<a<<std::endl;
			}
		}
 		std::cout<<"LPChecker::solve::cplex.getStatus() == IloAlgorithm::Infeasible"<<std::endl; abort();
		return -3;
	}else return -3;
	
}

//-------------------------------------------------------------------------------------------

int 
LPChecker::getSolution(MCND_solution*& mipsol){
 	double flow;
	double val;
	double cij;
 	
	mipsol = new MCND_solution(narcs+narcs*ndemands);
	mipsol->cost =0;
	IloNumArray x_(env);
	cplex.getValues(x_,x);
	for(int a=0;a<narcs;++a){
 		if(topo[a]==0){
 			mipsol->xy[a] = 0.0;
 			for (int k = 0; k < ndemands; ++k)
				mipsol->xy[narcs+k*narcs+a] = 0.0;
 		}else{
 			flow = 0.0;
 			for (int k = 0; k < ndemands; ++k){
				val = x_[a*ndemands+k];
				cij = data->arcs[a].c[k]*val;
				mipsol->xy[narcs+k*narcs+a] = val;
				mipsol->cost += cij;
				flow += val;
			}
			if(flow>1e-10){
				mipsol->xy[a] = 1.0;
				mipsol->cost += data->arcs[a].f;
			}
 		} 
	}
	x_.end();
 	//std::cout<<"feasibiliy integer sol value: "<<mipsol->cost<<std::endl;
	return 1;
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
LPChecker::check_feas(bool restore){
	int ret=0;
	double val;
	bool feas=true;
	std::list<Pair2> unrouted;
	IloNumArray x_(env);
	cplex.getValues(x_,xextra);
	for (int k = 0; k < ndemands; ++k){
		val = x_[k];
		if(val>1e-8){ 
			feas=false; 
			unrouted.push_back(Pair2(k, val));
			//std::cout<<"comm: "<<k<<" val: "<<val<<std::endl;
		}
	}
	if(!feas){
		if(restore){
			unrouted.sort(compPair2());
			ret = final_feas(unrouted, x_);
			if(ret < 0) return -2;
			else return -1;
		}else return -2;
	}
  	return 0;
}

//---------------------------------------------------------------------------


int
LPChecker::final_feas( std::list<Pair2>& heap,  const IloNumArray & xextrasol){
     
    IloNumArray xsol(env);
	cplex.getValues(xsol,x);
    //std::vector<double> primal(narcs*ndemands,0.0);
    std::vector<int> grid(nnodes*nnodes,-1);
    std::vector<double> wij(narcs,0.0);
    std::vector<double> uij(narcs,0.0);
    std::vector<Pair1> preced(nnodes);
    int arc;
    int k=0;
    double epsP=0.0;
    
    for(int a=0;a<narcs;++a){
        if(topo[a]==1){
        	uij[a] = data->arcs[a].capa;
            for (int k = 0; k < ndemands; ++k)
                uij[a] -= xsol[a*ndemands+k];
        }
    }

    //------- main loop ----------
    while(!heap.empty()){
        k = heap.front().fst;
        epsP =  heap.front().snd;
        //std::cout<<"k "<<k<<" "<<epsP<<" : "<<data->d_k[k].O<<" , "<<data->d_k[k].D<<std::endl;
       for(int a=0;a<narcs;++a){
            if(epsP <= uij[a]){
                //wij[a]= 1e-30;
                //std::cout<<"arc: "<<a<<" : "<<data->arcs[a].i<<" , "<<data->arcs[a].j<<std::endl;
                grid[(data->arcs[a].i-1)*nnodes+data->arcs[a].j-1]= a;
                if(topo[a]==0) wij[a]= data->arcs[a].c[k]*epsP + data->arcs[a].f*epsP/uij[a];
                else wij[a] = data->arcs[a].c[k]*epsP;
                
            }else grid[(data->arcs[a].i-1)*nnodes+data->arcs[a].j-1]=-1;
        }
        dijkstra(data->d_k[k].O, data->d_k[k].D, grid, preced, wij);//find path
        
        //path retrival
        int i,j;
        int contnodes=1;
        j=data->d_k[k].D;
        while (j != data->d_k[k].O){
            i = preced[j-1].node;
            arc = preced[j-1].pos;
             
            if(arc<0){
                //std::cout<<"ERROR 0/0: infeasible heuristic solution !!!!! "<<"commdoity "<<k<<" node: "<<i<<" "<<arc<<"/"<<narcs<<std::endl;
                return -1; //no feasible solution found
            }else if (uij[arc]<=1e-10 || contnodes>nnodes) { //check if impossible path
                //std::cout<<"ERROR 0/1: infeasible heuristic solution !!!!! "<<"commdoity "<<k<<" node: "<<i<<" "<<arc<<"/"<<narcs<<std::endl;
                return -1; //no feasible solution found
            }
            
            if(topo[arc]==0)topo[arc]=1;
            //primal[k*narcs+arcij] += epsP;
            uij[arc]-=epsP;
            // if(dmand==31)std::cout<<data->arcs[arcij].i<<"-"<<data->arcs[arcij].j<<std::endl;
            if(uij[arc]<1e-8){//arc is saturated
                wij[arc] = 1e30; //block
                grid[(i-1)*nnodes+j-1]=-1;
            }
            j = i;
            contnodes++;
        }
        heap.pop_front();
    }
    
    grid.clear();
    wij.clear();
    uij.clear();
    preced.clear();
    return 1;
}



//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

LPChecker::~LPChecker() {
    
    try {
    	delete [] map;
    	for(int i=tabusz;i--;){
			if(tabu[i]) delete tabu[i];
		}
		tabu.clear();
		unfx.clear();
		topo.clear();
		
        x.endElements();
        xextra.endElements();
 		cplex.clearModel();
		model.end();
		cplex.end();
		fobj.end();
        env.end();
    } catch (IloException& e) {
        std::cerr << "ERROR: " << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Error" << std::endl;
    }
}








