#include "MCND_checklpfeas.hpp"



#define INF 9999999
ILOSTLBEGIN

//---------------------------------------------------------------------------

void LPFeasChecker::set_parameters() {
	cplex.setParam(IloCplex::Threads,1);
    //cplex->setParam(IloCplex::RootAlg, 2);
    //cplex->setParam(IloCplex::NodeAlg, 3);
    cplex.setParam(IloCplex::ClockType, 1);
    //cplex->setParam(IloCplex::MIPDisplay, 4);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::DataCheck, CPX_DATACHECK_ASSIST);
    cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds
    
}

//---------------------------------------------------------------------------

void
LPFeasChecker::initialize(const Data* d){
    data = d;
    ndemands = d->ndemands;
    nnodes = d->nnodes;
    narcs =d->narcs;
    
    set_parameters();
    
    IloExpr obj1(env);
    IloExpr obj2(env);
    x = IloNumVarArray(env);
    for(int a=0;a<narcs;++a){
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(env, 0.0, IloInfinity));
            obj1 += x[a*ndemands+k];
            obj2 += data->arcs[a].c[k]*x[a*ndemands+k];
        }
    }
    fobj1 = IloMinimize(env, obj1);
    fobj2 = IloMinimize(env, obj2);
    model.add(fobj1);
    obj1.end();
    obj2.end();

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

//---------------------------------------------------------------------------
int
LPFeasChecker::solve_opt(const BCP_vec<BCP_var*>& vars, const double * topo) {
    IloNum ub;
    
    if(topo){
    	for(int a=0;a<narcs;++a){
    		if(vars[a]->ub()<=0.5){ ub = 0.0; //std::cout<<"close: "<<a<<" "<<colub[a]<<std::endl;
			}else if(vars[a]->lb()<=0.5 && topo[a]<=0.5){ ub = 0.0; //std::cout<<"close: "<<a<<" "<<colub[a]<<std::endl;
			}else{ ub = IloInfinity;}
		
			for(int k=0;k<ndemands;++k){
				x[a*ndemands+k].setUB(ub);
			}
		}
    }else{
		for(int a=0;a<narcs;++a){
			if(vars[a]->ub()<=0.5){ ub = 0.0; //std::cout<<"close: "<<a<<" "<<colub[a]<<std::endl;
			}else{ ub = IloInfinity;}
		
			for(int k=0;k<ndemands;++k){
				x[a*ndemands+k].setUB(ub);
			}
		}
    }
    cplex.getObjective().setExpr(fobj2);
	return solve();

}


//---------------------------------------------------------------------------

int LPFeasChecker::solve_feas(const double *collb, const double * colub){
    IloNum ub;
    for(int a=0;a<narcs;++a){
        if(colub[a]<=0.5){ ub = 0.0; //std::cout<<"close: "<<a<<" "<<colub[a]<<std::endl;
        }else{ ub = IloInfinity;}
        
        for(int k=0;k<ndemands;++k){
        	if(ub>0)
            	x[a*ndemands+k].setUB(colub[narcs+k*narcs+a]);
            else x[a*ndemands+k].setUB(0.0);
        }
    }
    
	cplex.getObjective().setExpr(fobj1);
	return solve();
}

//---------------------------------------------------------------------------

int LPFeasChecker::solve(){
	//cplex.exportModel("rl.lp");
	cplex.solve();

	if(cplex.getStatus() == IloAlgorithm::Infeasible){
        std::cout<<"infeasible"<<std::endl;
		return -1;
	}else if(cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){		 
        std::cout<<"InfeasibleOrUnbounded"<<std::endl;
		return -2;
		
	}else if(cplex.getStatus() == IloAlgorithm::Optimal){
		 //std::cout<<"optimal "<<cplex.getObjValue()<<std::endl;
		return 0;
	}
	else return -3;
}

//---------------------------------------------------------------------------


int 
LPFeasChecker::getSolution(const BCP_vec<BCP_var*>& vars, double * sol, double& solvalue, double &fathmval){
 	double flow;
	double val;
	double cij;
	int ret=-1; 
	solvalue=0;
	IloNumArray x_(env);

	cplex.getValues(x_,x);
	for(int a=0;a<narcs;++a){
 		flow = 0.0;
		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k];
			sol[narcs+k*narcs+a] = val;
			cij = data->arcs[a].c[k]*val;
			solvalue+=cij;
			fathmval+=cij;
			flow += val;
		}
		if(flow>1e-10){
			sol[a] = 1.0;
			solvalue+=data->arcs[a].f;
			if(vars[a]->lb()<0.5 &&  vars[a]->ub()>0.5) ++ret;
		}//else if(vars[a]->lb()>=0.5){//std::cout<<"LPFeasChecker::getSolution FATHOM"<<std::endl; ret = -1;
	}
	
 	std::cout<<"feasibiliy integer sol value: "<<solvalue<<std::endl;
	return ret;
}

//---------------------------------------------------------------------------

LPFeasChecker::~LPFeasChecker() {

	try {
		cplex.clearModel();

		x.end();
         
        model.end();
		cplex.end();
		env.end();
		

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







