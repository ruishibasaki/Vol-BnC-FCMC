#include "MCND_pump.hpp"



#define INF 9999999
ILOSTLBEGIN


void
Pump::set_data(const Data * d, double alph_init){
	
	data = d;
	nnodes = data->nnodes;
	ndemands = data->ndemands;
	narcs = data->narcs;
	alpha = alph_init;
	topo.resize(narcs);
	tabu.resize(1000,0);
	branch_candidates.resize(narcs);
	
	tern=0;
	maxunfix = narcs*0.1;
	set_parameters();
	
    int sizeOfInt=8*sizeof(unsigned int);
    sizeOfIdSeq = (narcs/sizeOfInt)+1;
    
	map = new Pair2[narcs];
    for(int i=0;i<narcs;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
    }
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
Pump::initialize(const Data * d,  const MCND_solution* best_sol_, double alph_init){
    set_data(d, alph_init);
    best_sol = best_sol_;
    IloExpr obj(env);
    //x = IloNumVarArray(env);
    //y = IloNumVarArray(env);
    for(int a=0;a<narcs;++a){
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(env, 0.0, IloInfinity));
            obj += x[a*ndemands+k];
         }
         y.add(IloNumVar(env, 0.0, 1.0));
         obj+= y[a];
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
Pump::make_topo(int * fixd, int& closed,  const double * x ,const BCP_vec<BCP_var*>& vars){
	int cont=0;
	closed = 0;
	szunfix=0;
	
	for(int a=narcs; a--;){
        if(vars[a]->lb()==1.0){
            topo[a]=2;
        }else if(vars[a]->ub()==1.0){
            if(x[a]>=0.9){ 
            	topo[a]=2;
            }else if(x[a]>=0.1){
            	if(rand()%2==1){
            		if(best_sol->xy[a]==1){
            			topo[a]=-1;
            		}else{
            			topo[a]=1;
            		}
            	}else if(x[a]>0.5){
            		topo[a]=-1;
            	}else{
            		topo[a]=1;
            	} 
            	++szunfix;
            	unfx.push_back(a);
            }else{
            	fixd[closed++]=a; 
            	topo[a]=0;
            }//std::cout<<"cand: "<<a<<" "<<x[a]<<std::endl;}
        	++cont;
        }else{ 
        	fixd[closed++]=a;  
        	topo[a]=0;
        }
    }
    if(szunfix > maxunfix){
    	return -1;
    }
    unsigned int* seqtopo = new unsigned int[sizeOfIdSeq];
	std::fill(seqtopo, seqtopo+sizeOfIdSeq, 0);
    if(check_tabu(seqtopo)){
    	try_perturbation(seqtopo);
    	if(check_tabu(seqtopo)){
    		delete [] seqtopo;
     		return -1;
     	}
    }
    
    return cont;
}

//-------------------------------------------------------------------------------------------

bool
Pump::check_tabu(unsigned int* seqtopo){
	unsigned int* titem;
	bool equal;
	Pair2* maptitem=0;

	for(int a=narcs; a--;){
		maptitem = &map[a];
		if(topo[a]==2 || topo[a]==-1){
			setBit(seqtopo, maptitem->fst, maptitem->snd);
		}
	}

	for(int i=tern;i--;){
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
	
    return false;
}

//-------------------------------------------------------------------------------------------

void
Pump::try_perturbation(unsigned int* seqtopo){
	int arc;
	Pair2* maptitem=0;
	std::cout<<"Pump::try_perturbation"<<std::endl;
 	for(int a=0;a<szunfix;++a){
		arc = unfx[a];
		std::cout<<"arc: "<<arc<<std::endl;
		if(best_sol->xy[arc]==0 && topo[arc]==-1){
			topo[arc]=1;
			maptitem = &map[arc];
			clearBit(seqtopo, maptitem->fst, maptitem->snd);
		}else if(best_sol->xy[arc]==1 && topo[arc]==1){
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

void  
Pump::reset( const BCP_vec<BCP_var*>& vars){
	int arc;
	double norm=0;
	alpha = szunfix/double(maxunfix);
	
	if(szunfix ==0 ){
		factorxy=1;
		factorp=0;
	}else{
		for(int a=narcs;a--;){
			if(topo[a]==0) continue;
 			if(topo[a]!=2) norm+= pow(data->arcs[a].f,2);
			//std::cout<<data->arcs[arc].f<<std::endl;
			for (int k = 0; k<ndemands;++k)
				norm+=pow(data->arcs[a].c[k],2);
		}
		norm = sqrt(norm);
		factorxy = (1.0 - alpha)/norm;
		factorp = alpha/(sqrt(double(szunfix)));
	}
}

//-------------------------------------------------------------------------------------------

void Pump::create_model( const BCP_vec<BCP_var*>& vars) {	
	Pair2 item;
	int arc;
	double c;
	
	IloExpr obj(env);
	for(int a=narcs ; a--; ){
		if(topo[a]==0){
			//std::cout<<"out: "<<a<<std::endl;
			y[a].setUB(0.0);
			y[a].setLB(0.0);
			continue;
		} 
		for (int k = ndemands; k-- ;){
			//x[a*ndemands+k].setUB(vars[narcs+k*narcs+a]->ub());
			obj += factorxy*data->arcs[a].c[k]*x[a*ndemands+k];
		}
		if(topo[a] == 2){
			y[a].setUB(1.0);
			y[a].setLB(1.0);
		}else{
			y[a].setUB(1.0);
			y[a].setLB(0.0);
			if(topo[a] == -1) c = factorxy*data->arcs[a].f - factorp;
			else if(topo[a] == 1) c = factorxy*data->arcs[a].f  + factorp;
			obj += c*y[a];
		}
	}

	fobj = IloMinimize(env, obj);
    cplex.getObjective().setExpr(fobj);
	obj.end();
	
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
			ub = vars[narcs+k*narcs+arc]->ub();
			if(x_[arc*ndemands+k] - ub*y_[arc]> 1e-10 ){
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
Pump::solve(  const BCP_vec<BCP_var*>& vars,  double * xy, double & val){
	reset( vars);
	create_model( vars);
	//cplex.exportModel("t.lp");
	cplex.solve();
	
	
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		//std::cout<<"pump:: cplex.getStatus() == IloAlgorithm::Infeasible"<<std::endl;
		clear();
		return -1;
	}else if(cplex.getStatus() == IloAlgorithm::Unbounded){
			//std::cout<<"pump:: cplex.getStatus() == IloAlgorithm::Unbounded"<<std::endl;
		clear();
		return -2;
	}
	
	IloNumArray x_(env);
	IloNumArray y_(env);

	cplex.getValues(x_,x);
	cplex.getValues(y_,y);
	while(cut(vars, y_,x_)){
		cplex.solve();
		std::cout<<"after cut "<<cplex.getObjValue()<<std::endl;
		cplex.getValues(y_,y);
		cplex.getValues(x_,x);
	}
	//std::cout<<"final pump "<<cplex.getObjValue()<<std::endl;

	val = getSolution( xy, x_, y_ );
	
	x_.end(); 
	y_.end();
	clear();
	return 0;
	
}

//-------------------------------------------------------------------------------------------

double 
Pump::getSolution( double * xy,  const IloNumArray & x_, const IloNumArray & y_ ){
	
 	double flow;
	double solvalue=0;
	double val;
	for(int a=0;a<narcs;++a){
 		if(topo[a]==0) continue;
 		flow = 0.0;
 		branch_candidates[a]=-1;
		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k];
			xy[narcs+k*narcs+a] = val;
			solvalue+=data->arcs[a].c[k]*val;
			flow += val;
		}
		if(topo[a]!=2){
			val = y_[a];
			std::cout<<"heurarc: "<<a<<" val: "<<val<<" topo: "<<topo[a]<<std::endl;
			if(val>1e-10){
				xy[a] = 1.0;
 				solvalue+=data->arcs[a].f;
			}
			if(topo[a]==-1 && val < 0.9)  branch_candidates[a]= 1.0-val;
			else if(topo[a]==1 && val > 0.1) branch_candidates[a] = val;
		}else if(flow>1e-10){
 			xy[a] = 1.0;
			solvalue+=data->arcs[a].f;
		} 
	}
	
	
	//std::cout<<"sol value: "<<solvalue<<std::endl;
	return solvalue;
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void 
Pump::clear(){
	fobj.end();
	unfx.clear();
}


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
		branch_candidates.clear();
		for(int i=tabu.size() ; i--; ){
			if(tabu[i]!=0) delete [] tabu[i];
		}
		tabu.clear();
	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







