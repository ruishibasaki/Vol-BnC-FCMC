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
                constraint -= data->d_k[k].quantity;
            }else if( i ==  data->d_k[k].O){
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
	check_feas(seqtopo);
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

//-------------------------------------------------------------------------------------------

void
LPChecker::check_feas(unsigned int* seqtopo){
	/*int arc=0;
	bool feas;
	GlobalCut* gloc = globalcs->begin;
	for(int i = globalcs->sizeOfCollection; i--;){
		//std::cout<<"gloc_updt "<<gloc->serial_nmbr<<" "<<gloc->purgbl<<std::endl;
 		feas =false;
		for(int a=gloc->size;a--;){
			arc = gloc->vars[a];
			if(topo[arc]>0.5){
				feas =true;
				break;
			}	
		}
		if(!feas){ std::cout<<"LPChecker::check_feas not"<<std::endl; }
		gloc = gloc->next;
	}
	
	Cover* cov = covers->begin;
	for(int i = covers->sizeOfCollection; i--;){
		//std::cout<<"gloc_updt "<<gloc->serial_nmbr<<" "<<gloc->purgbl<<std::endl;
		feas =false;
		double sum=0;
		double comp= cov->get_rhs();
		int sz = cov->get_total_sz();
		for(int a=0;a<sz;++a){
			sum+= cov->gamma_at(a)*topo[cov->at(a)];
			if(sum>=comp){
				feas =true;
				break;
			}
		}
		if(!feas){ std::cout<<"LPChecker::check_feas not"<<std::endl; }
		cov = cov->next;
	}
	
	LocalCut* loc = localcs->begin;
	for(int i = localcs->sizeOfCollection; i--;){
		//std::cout<<"gloc_updt "<<gloc->serial_nmbr<<" "<<gloc->purgbl<<std::endl;
		feas =false;
		double sum=0;
		double comp= cov->get_rhs();
		int sz = cov->get_total_sz();
		for(int a=0;a<sz;++a){
			sum+= cov->gamma_at(a)*topo[cov->at(a)];
			if(sum>=comp){
				feas =true;
				break;
			}
		}
		if(!feas){ std::cout<<"LPChecker::check_feas not"<<std::endl; }
		cov = cov->next;
	}*/
	
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int
LPChecker::solve(int& closed, int*& fixd0, MCND_solution*& mipsol){
	
	make_model();
	cplex.solve();
	
	if(cplex.getStatus() == IloAlgorithm::Optimal){
		getSolution(mipsol);
        std::cout<<"LPChecker: "<<cplex.getObjValue()<<std::endl;
		return 0;
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
		return -1;
	} 
	else return -3;
	
}

//---------------------------------------------------------------------------


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








