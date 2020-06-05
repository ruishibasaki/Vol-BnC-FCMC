#include "MCND_pump.hpp"



#define INF 9999999
ILOSTLBEGIN


void
Pump::set_data(const Data * d){
	
	data = d;
	nnodes = data->nnodes;
	ndemands = data->ndemands;
	narcs = data->narcs;
 	unfx.resize(narcs);
	topo.resize(narcs);
	tabu.resize(1000,0);
 	
	tabusz = tern=0;
	szunfix=0;
	maxunfix = narcs*0.3;
	set_parameters();
	
    int sizeOfInt=8*sizeof(unsigned int);
    sizeOfIdSeq = (narcs/sizeOfInt)+1;
    
	map = new Pair2[narcs];
    for(int i=0;i<narcs;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
    }
    factory=0;
    for(int a=narcs;a--;){
		 if(factory<data->arcs[a].f)factory = data->arcs[a].f;
	}
 	factory *= 100;
}

//-------------------------------------------------------------------------------------------

void Pump::set_parameters() {
	cplex.setParam(IloCplex::Threads,1);
	//cplex->setParam(IloCplex::RootAlg, 2);
	//cplex->setParam(IloCplex::NodeAlg, 3);
	cplex.setParam(IloCplex::ClockType, 1);
	//cplex->setParam(IloCplex::MIPDisplay, 4);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds

}

//---------------------------------------------------------------------------

void
Pump::initialize(const Data * d,  const MCND_solution* best_sol_){
    set_data(d);
    best_sol = best_sol_;
    IloExpr obj(env);
    //x = IloNumVarArray(env);
    //y = IloNumVarArray(env);
    for(int a=0;a<narcs;++a){
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(env, 0.0, IloInfinity));
            obj += data->arcs[a].c[k]*x[a*ndemands+k];
         }
         y.add(IloNumVar(env, 0.0, 1.0));
         obj+=  data->arcs[a].f*y[a];
    }
    fobj = IloMinimize(env, obj);
    model.add(fobj);
    obj.end();
 
    for(int a=0;a<narcs;++a){
        IloExpr constraint(env);
        for (int k = 0; k < ndemands; k++) {
            constraint += x[a*ndemands+k];
        }
        constraint -= data->arcs[a].capa*y[a];
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
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
Pump::make_topo( const double * x ,const BCP_vec<BCP_var*>& vars){
	int rnd=0;
	int cont=0;
	int pertbd =0;
	double maxpertbd = narcs*0.01;
	maxpertbd = maxpertbd>0? maxpertbd: 1;
 	szunfix=0;
 	for(int a=narcs; a--;){
        if(vars[a]->lb()==1.0){
            topo[a]=-2;
         }else if(vars[a]->ub()==1.0){
            if(x[a]>=0.9){ 
            	topo[a]=-2;
            }else if(x[a]>=0.1){
            	rnd = rand()%2;
				if(rnd==1 && pertbd<maxpertbd){
					rnd = ((rand()%101)/100.0 <= x[a]) ? 1 : 0;
					if(rnd){
						topo[a]=-1;
 					}else{
						topo[a]=1;
					}
					++pertbd;
				}else{
					if(best_sol->xy[a]>0.5){
						topo[a]=-1;
 					}else{
						topo[a]=1;
					}
				}	 
				unfx[szunfix++] = a;   
            }else{
             	topo[a]=0;
            }//std::cout<<"cand: "<<a<<" "<<x[a]<<std::endl;}
        	++cont;
        }else{ 
         	topo[a]=0;
        }
    }
    return cont;
    
}

//-------------------------------------------------------------------------------------------

int
Pump::validate_topology( ){
	unsigned int* seqtopo = new unsigned int[sizeOfIdSeq];
	std::fill(seqtopo, seqtopo+sizeOfIdSeq, 0);
	Pair2* maptitem=0;
 	for(int a=narcs; a--;){
		maptitem = &map[a];
        if(topo[a]<0.0){
             setBit(seqtopo, maptitem->fst, maptitem->snd);
        }
    }
     
    if(check_tabu(seqtopo)){
    	try_perturbation(seqtopo);
    	if(check_tabu(seqtopo)){
    		delete [] seqtopo;
     		return -1;
     	}
    }
    return 1;
}

//-------------------------------------------------------------------------------------------

bool
Pump::check_tabu(unsigned int* seqtopo){
	unsigned int* titem;
	bool equal;
	Pair2* maptitem=0;

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
Pump::try_perturbation(unsigned int* seqtopo){
	int arc;
	Pair2* maptitem=0;
	std::cout<<"Pump::try_perturbation"<<std::endl;
	std::random_shuffle(unfx.begin(), unfx.begin()+szunfix);
 	for(int a=0;a<szunfix;++a){
		arc = unfx[a];
		//std::cout<<"arc: "<<arc<<std::endl;
		if(topo[arc]==-1){
			topo[arc]=1;
			maptitem = &map[arc];
			clearBit(seqtopo, maptitem->fst, maptitem->snd);
		}else if(topo[arc]==1){
			topo[arc]=-1;
			maptitem = &map[arc];
			setBit(seqtopo, maptitem->fst, maptitem->snd);
		}
		if(rand()%2)++a;
	}

}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void Pump::create_model( ) {	
	Pair2 item;
	int arc;
	double c;
	
	IloExpr obj(env);
	for(int a=narcs ; a--; ){
		if(topo[a]==0){
			y[a].setUB(0.0);
			y[a].setLB(0.0);
			continue;
		} 
		for (int k = ndemands; k-- ;){
 			obj += data->arcs[a].c[k]*x[a*ndemands+k];
		}
		if(topo[a] == -2){
			y[a].setUB(1.0);
			y[a].setLB(1.0);
			continue;
		}
			
		y[a].setUB(1.0);
		y[a].setLB(0.0);
		 
		obj += factory*topo[a]*y[a];
	}

    cplex.getObjective().setExpr(IloMinimize(env, obj));
	obj.end();
	
}

//-------------------------------------------------------------------------------------------

void Pump::check_feas_model() {	
	Pair2 item;
	int arc;
	double c;
	
	IloExpr obj(env);
	for(int a=narcs ; a--; ){
		if(topo[a]<0){
			//std::cout<<"out: "<<a<<std::endl;
			y[a].setUB(1.0);
			y[a].setLB(1.0);
		}else{
			y[a].setUB(0.0);
			y[a].setLB(0.0);
		} 
	}
    cplex.getObjective().setExpr(fobj);	
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


int 
Pump::cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & y_, const IloNumArray & x_){
	int arc;
	int cont=0;
	double ub;
	for(int a=0;a<szunfix;++a){
		arc = unfx[a];
		for (int k = 0; k < ndemands; ++k){
			ub = data->arcs[arc].b[k];
			if(x_[arc*ndemands+k] - ub*y_[arc]> 1e-2){
				IloExpr constraint(env);
				constraint -= x[arc*ndemands+k];
				constraint+= ub*y[arc];
				++cont;
				cutstrong.add((constraint >= 0));
				model.add(cutstrong[cutstrong.getSize()-1]);
				constraint.end();
			}
		}
	}
	return cont;
}



//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


int 
Pump::solve( int*& fixd0, int& closed,  const BCP_vec<BCP_var*>& vars,  MCND_solution*& mipsol){
	bool feas=true;
	
	/*check_feas_model();
	cplex.solve();
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		std::cout<<"pump:: cplex.getStatus() == IloAlgorithm::Infeasible1"<<std::endl;
		feas=false;
	} 
	if(feas){
		IloNumArray x_(env);
		cplex.getValues(x_,x);
		getSolution(  x_, mipsol );
		x_.end(); 
 		return 1;
	}*/
	std::cout<<"Pump::solve trypump? "<<szunfix<<" "<<maxunfix<<std::endl;
	/*if(szunfix > maxunfix ){
		get_closed(fixd0, closed, true);
     	return 0;
    }*/
    create_model();
    IloNumArray x_(env);
	IloNumArray y_(env);
	bool dontbreak = true;
    while(dontbreak){
		//cplex.exportModel("t.lp");
		cplex.solve();
	
		if(cplex.getStatus() == IloAlgorithm::Infeasible){
			std::cout<<"pump:: cplex.getStatus() == IloAlgorithm::Infeasible2"<<std::endl;
			get_closed(fixd0, closed, false);
			return -1;
		} 
	
		cplex.getValues(x_,x);
		cplex.getValues(y_,y);
		while(cut(vars, y_,x_)){
			cplex.solve();
			//std::cout<<"after cut "<<cplex.getObjValue()<<std::endl;
			if(cplex.getStatus() == IloAlgorithm::Infeasible){
				std::cout<<"pump:: cplex.getStatus() == IloAlgorithm::Infeasible2"<<std::endl;
				get_closed(fixd0, closed, false);
				return -2;
			} 
			cplex.getValues(y_,y);
			cplex.getValues(x_,x);
		}
		dontbreak = false;
		double solpump=0;
		for(int a=0;a<szunfix;++a){
			int arc = unfx[a];
			if(topo[arc]<0)solpump += factory - factory*y_[arc];
			else solpump += factory*y_[arc];
			if(y_[arc] < 0.5 && topo[arc]<0){
				std::cout<<"y "<<arc<<" : "<<y_[arc]<<" topo: "<<topo[arc]<<std::endl;
				cplex.getObjective().setLinearCoef(y[arc], factory);
				topo[arc] = 1;
				dontbreak = true;
			}else if(y_[arc] > 0.5 && topo[arc]>0){
				std::cout<<"y "<<arc<<" : "<<y_[arc]<<" topo: "<<topo[arc]<<std::endl;
				cplex.getObjective().setLinearCoef(y[arc], -factory);
				topo[arc] = -1;
				dontbreak = true;
			}
		}
		std::cout<<"final pump "<<solpump<<std::endl;
    }
 	
	getSolution( x_, mipsol );
	//get_closed(fixd0, closed, true);
	x_.end(); 
	y_.end();
 	return 1;
}

//-------------------------------------------------------------------------------------------

double 
Pump::getSolution(const IloNumArray& x_, MCND_solution*& mipsol ){
	double flow;
	double val;
	double cij;
 	
	mipsol = new MCND_solution(narcs+narcs*ndemands);
	mipsol->cost =0;
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
  	//std::cout<<"feasibiliy integer sol value: "<<mipsol->cost<<std::endl;
	return 1;
}

//-------------------------------------------------------------------------------------------

void 
Pump::get_closed(int*& fixd0, int& closed, bool onlyx){
	fixd0 = new int[narcs];
	closed=0;
	if(onlyx){
		for(int a=0;a<narcs;++a)
			if(topo[a]>=0)
				fixd0[closed++] = a;
	}else{
		for(int a=0;a<narcs;++a)
			if(topo[a]==0)
				fixd0[closed++] = a;
	}
	
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

Pump::~Pump() {

	try {
		cutstrong.endElements();
		x.endElements();
		y.endElements();
		cplex.clearModel();
		model.end();
		cplex.end();
		fobj.end();

		env.end();
		unfx.clear();
		topo.clear();
 		for(int i=tabu.size() ; i--; ){
			if(tabu[i]!=0) delete [] tabu[i];
		}
		tabu.clear();
        delete [] map;
	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







