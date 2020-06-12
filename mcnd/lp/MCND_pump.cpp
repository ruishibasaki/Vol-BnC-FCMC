#include "MCND_pump.hpp"



#define INF 9999999
ILOSTLBEGIN


void
Pump::set_data(const Data * d){
	
	data = d;
	nnodes = data->nnodes;
	ndemands = data->ndemands;
	narcs = data->narcs;
	volsolver.ext_initializer(narcs+narcs*ndemands, nnodes*ndemands);
	volsolver.dual_ub = 1e31;
	volsolver.dual_lb = -1e31;
 
    fxone.resize(narcs);
 	unfx.resize(narcs);
	topo.resize(narcs);
  	
 	szunfix=0;
 	round =0.3;
	maxunfix = narcs*0.3;
	set_parameters();
	
    factory=0;
    for(int a=narcs;a--;){
		 if(factory<data->arcs[a].f)factory = data->arcs[a].f;
	}
 	//factory *= 100;
}

//-------------------------------------------------------------------------------------------

void Pump::set_parameters() {
	cplex.setParam(IloCplex::Threads,1);
	//cplex->setParam(IloCplex::RootAlg, 2);
	//cplex->setParam(IloCplex::NodeAlg, 3);
	cplex.setParam(IloCplex::ClockType, 1);
	//cplex->setParam(IloCplex::MIPDisplay, 4);
	cplex.setOut(env.getNullStream());
	//cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds
	cplex.setParam(IloCplex::Param::Advance, 0);

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
	const double *besttopo=0;
	maxpertbd = maxpertbd>0? maxpertbd: 1;
 	szunfix=0;szfxone=0;
 	if(best_sol==0){
 		 if(altsol!=0) besttopo = altsol;
 	}else if(!best_sol->onlyvalue) besttopo = best_sol->xy;
 	else if(altsol!=0){besttopo = altsol;}
 	
 	for(int a=narcs; a--;){
        if(vars[a]->lb()==1.0){
            topo[a]=-2;
            fxone[szfxone++] = a;
         }else if(vars[a]->ub()==1.0){
            /*if(x[a]>=0.99){ 
            	topo[a]=-2;
            	fxone[szfxone++] = a;
            }else*/ if(x[a]>=0.01){
            	rnd = rand()%2;
				if((rnd==1 && pertbd<maxpertbd) || (besttopo==0)){
					rnd = ((rand()%101)/100.0 <= x[a]) ? 1 : 0;
					if(rnd){
						topo[a]=-1;
 					}else{
						topo[a]=1;
					}
 					++pertbd;
				}else{
					if(besttopo[a]>0.5){
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
	if(best_sol->onlyvalue && altsol==0){
    	altsolval=1e30;
    }
     
	bool dontbreak = true;
	std::vector<int>tabu(narcs,0);
    while(dontbreak){
		if(volsolve()<0){
			get_closed(fixd0, closed, false);
			return -1;
		}
		
		dontbreak = false;
		double solpump=0;
		for(int a=0;a<szunfix;++a){
			int arc = unfx[a];
			if(topo[arc]<0)solpump += factory - factory*volsolver.psol[a];
			else solpump += factory*volsolver.psol[a];
			if(volsolver.psol[a] <= round && topo[arc]<0 && tabu[arc]<3){
				//std::cout<<"y "<<arc<<" : "<<volsolver.psol[a]<<" topo: "<<topo[arc]<<std::endl;
				//cplex.getObjective().setLinearCoef(y[arc], factory);
				tabu[arc]++;
				topo[arc] = 1;
				dontbreak = true;
			}else if(volsolver.psol[a] >= 0.3 && topo[arc]>0 && tabu[arc]<3){
				//std::cout<<"y "<<arc<<" : "<<volsolver.psol[a]<<" topo: "<<topo[arc]<<std::endl;
				//cplex.getObjective().setLinearCoef(y[arc], -factory);
				topo[arc] = -1;
				tabu[arc]++;
				dontbreak = true;
			}
		}
		std::cout<<"final pump "<<solpump<<" volume: "<<volsolver.value<<std::endl;
    }
    
	 
 	check_feas_model();
	cplex.solve();
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		std::cout<<"pump:: cplex.getStatus() == IloAlgorithm::Infeasible1"<<std::endl;
		feas=false;
	} 
	if(feas && cplex.getStatus() != IloAlgorithm::Optimal){
		return -10;
	}
	if(feas){
		IloNumArray x_(env);
		cplex.getValues(x_,x);
		getSolution(  x_, mipsol );
		x_.end(); 
		if(best_sol != 0 && mipsol->cost > best_sol->cost)
			round*=1.8; if(round>0.5) round=0.5;
		else if(best_sol ==0) return 1;
		if(best_sol->onlyvalue && mipsol->cost < altsolval){
			altsolval = mipsol->cost;
			if(altsol==0) altsol = new double [narcs];
			for(int a=narcs;a--;)altsol[a] = mipsol->xy[a];
		}
  		return 1;
	}else{ 
		get_closed(fixd0, closed, true); 
		round*=0.8; if(round<0.3) round=0.3;
		return 0;
	}
	
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
		for(int a=0;a<narcs;++a){
			if(topo[a]>=0)
				fixd0[closed++] = a;
 		}
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
 		fxone.clear();
		unfx.clear();
		topo.clear();
 		if(altsol !=0) delete []  altsol;
 	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

//====================================================================
//====================================================================
//  volume hooks
//====================================================================
//====================================================================

int 
Pump::volsolve(){
	
    volsolver.dsol=0;
    sznz = szunfix+szfxone;
    volsolver.active_size = ndemands*nnodes;
    volsolver.psize = szunfix + ndemands*sznz;
    volsolver.active_psize = szunfix;
    int retval = volsolver.solve(*this, true);
    if(retval<0) return -1;
    return 1;

}

//-------------------------------------------------------------------------------------------

int
Pump::solve_subproblem(const VOL_dvector& xstar,
                                        const VOL_dvector& dualu,  VOL_dvector& rc,
                                        double& lcost, VOL_dvector& x,
                                        double& pcost){
    int arc;
    double cost_a;
    pcost= 0;
    lcost=0;
	for(int k=0; k<data->ndemands; ++k){
        lcost += ( dualu[k*nnodes + data->d_k[k].D-1] - dualu[k*nnodes + data->d_k[k].O-1])*data->d_k[k].quantity;
    }    
    for(int a=szunfix; a--; ){
    	//if(xstar[a]>1)std::cout<<"WHAT? "<<xstar[a]<<std::endl;
        arc = unfx[a];
        cost_a = knapsack(a, true, rc.v, x.v);
        rc[a] += cost_a; 
        if( rc[a] < 0.0){
            x[a] =1.0;
            lcost += rc[a];
        }else{
            x[a]=0.0;             
            for(int k=0; k<ndemands; ++k)
                x[szunfix + k*sznz + a]=0.0;
        }  
    }
    for(int a=szfxone; a--;){
        arc = fxone[a];
        cost_a = knapsack(a, false, rc.v, x.v);
        lcost += cost_a;
    }
		 
    
    return 0;
}

//-------------------------------------------------------------------------------------------

double
Pump::knapsack(int a, bool none, const double * rc, double* x){
    double kpsack =0;
    double fillUp =0;
    double rcost, flow;
    int arc, base;
    if(none){ arc = unfx[a]; base = szunfix;}
    else{  arc = fxone[a];  ; base = 2*szunfix;}
    std::list<HeapCell> heap;
    //get reduced cost for each commodity in arc e
    for(int k=ndemands; k--; ){
    	rcost = rc[base + k*sznz + a];
        if(rcost<0.0){
            heap.push_back(HeapCell(k, rcost));
        }else x[base + k*sznz + a] = 0;
        
    }
    
    if(  heap.empty() ){ 
    	 return 0;
    }
    heap.sort(comp());
    //std::stable_sort(heap.begin(), heap.end(), comp());
    
    
    int comm;
    double capa = data->arcs[arc].capa;
    while(heap.size()>0){ 
        comm = heap.back().k;
        if(fillUp < capa){
        	flow = std::min((capa - fillUp), data->arcs[arc].b[comm]);
            x[base + comm*sznz + a] = flow;
            fillUp += flow;
            kpsack += heap.back().rc_ * flow;
        }else x[base + comm*sznz + a] =0;
        heap.pop_back();
    }
    return kpsack;
}

//-----------------------------------------------------------------------

int
Pump::compute_rc(const VOL_dvector& dualu, VOL_dvector& rc, int actvSSz){
    const Arc* item;
    int arc;
    for(int a=szunfix; a--;){
        arc = unfx[a];
        item = &data->arcs[arc];
        if(topo[arc]>0)rc[a] = item->f; //factory;
        else rc[a] = 0.0;
        for(int k=ndemands; k--; ){
            rc[szunfix + k*sznz + a] = item->c[k]  - dualu[k*nnodes + item->j-1] + dualu[k*nnodes + item->i-1];
        }
    }
    for(int a=szfxone; a--;){
        arc = fxone[a];
        item = &data->arcs[arc];
        for(int k=ndemands; k-- ;){
            rc[2*szunfix + k*sznz + a] = item->c[k] - dualu[k*nnodes + item->j-1] + dualu[k*nnodes + item->i-1];
        }
    }

    
    
    
     return 0;
}

//-----------------------------------------------------------------------

int
Pump::compute_sg(const VOL_dvector& x, int actvSSz, VOL_dvector& v){
    //std::cout<<"compute_sg"<<std::endl;
    for(int n=actvSSz; n--; )
    	v[n] = 0;
    	
    int arc, basek, basef;
    int id;
    const Arc* item; const Demand* itemd;
    for(int k=ndemands; k--; ){
        itemd = &data->d_k[k];
        basek = k*nnodes;
        v[basek+ itemd->O-1] -= data->d_k[k].quantity;
        v[basek + itemd->D-1] +=  data->d_k[k].quantity;
        for(int a=0; a<szunfix; ++a){
            arc = unfx[a];
            item = &data->arcs[arc];
            basef = szunfix+k*sznz;
            v[basek + item->i-1] += x[basef+ a];
            v[basek + item->j-1] -= x[basef+ a];
        }
        for(int a=0; a<szfxone; ++a){
            arc = fxone[a];
            item = &data->arcs[arc];
            basef = 2*szunfix+k*sznz;
            v[basek + item->i-1] += x[basef+ a];
            v[basek + item->j-1] -= x[basef+ a];
        }
    }
   
    return 0;
}





