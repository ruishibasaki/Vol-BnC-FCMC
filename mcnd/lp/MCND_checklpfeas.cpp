#include "MCND_checklpfeas.hpp"



#define INF 9999999
ILOSTLBEGIN

//---------------------------------------------------------------------------

void LPFeasChecker::set_parameters() {
    //cplex->setParam(IloCplex::Threads,0);
    //cplex->setParam(IloCplex::RootAlg, 2);
    //cplex->setParam(IloCplex::NodeAlg, 3);
    cplex.setParam(IloCplex::ClockType, 1);
    //cplex->setParam(IloCplex::MIPDisplay, 4);
    cplex.setOut(env.getNullStream());
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
    
    IloExpr obj(env);
    for(int a=0;a<narcs;++a){
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(env));
            obj += x[a*ndemands+k];
        }
    }
    model.add(IloMinimize(env, obj));
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

//---------------------------------------------------------------------------

void LPFeasChecker::create_model(const double *collb, const double * colub) {
    IloNum ub;
    for(int a=0;a<narcs;++a){
        if(colub[a]<=0.5){ ub = 0.0; std::cout<<"close: "<<a<<std::endl;
        }else{ ub = IloInfinity;}
        
        for(int k=0;k<ndemands;++k){
            x[a*ndemands+k].setUB(ub);
        }
    }
}


//---------------------------------------------------------------------------

int LPFeasChecker::solve(const double *collb, const double * colub){
    //std::cout<<"VERIFY FEAS!!"<<std::endl;
    //cplex.exportModel("rl.lp");
    create_model( collb,   colub);
	cplex.solve();
	
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		 
        std::cout<<"infeasible"<<std::endl;
		return -1;
	}else if(cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
		 
        std::cout<<"InfeasibleOrUnbounded"<<std::endl;
		return -2;
		
	}else if(cplex.getStatus() == IloAlgorithm::Optimal){
		 
		return 0;
	}
	else return -3;
	
}

//---------------------------------------------------------------------------

LPFeasChecker::~LPFeasChecker() {

	try {
		cplex.clearModel();

		x.end();
        flow.end();
        capa.end();
        model.end();
		cplex.end();
		env.end();
		

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







