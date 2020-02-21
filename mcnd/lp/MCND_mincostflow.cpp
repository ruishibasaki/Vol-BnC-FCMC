#include "MCND_mincostflow.hpp"



#define INF 9999999
ILOSTLBEGIN


//---------------------------------------------------------------------------

MinCostFLow::MinCostFLow() :
	cplex(env),
	model(0){
	
}

//---------------------------------------------------------------------------

void MinCostFLow::create_model(const std::deque<int> & topo, const Data * data) {
	const Arc* item;
    const Demand* itemd;
	int nnodes = data->nnodes;
	int ndemands = data->ndemands;
	int narcs = data->narcs;
	
	int sznz = topo.size();
	int arc;
	    
    set_parameters();
    
    model = new IloModel(env);
    x = IloNumVarArray(env, sznz*ndemands);
    IloExpr obj(env);
	for(int a=sznz;a--;){
		for (int k = 0; k < ndemands; ++k){
			x[k*sznz+a] = IloNumVar(env);
			obj += x[k*sznz+a];
		}
	}
    model->add(IloMinimize(env, obj));
	obj.end();
	
	for(int a=sznz;a--;){
		arc =  topo[a];        
		IloExpr constraint(env);
		for (int k = 0; k < ndemands; k++) {
			constraint += x[k*sznz+a];
		}
		constraint -= data->arcs[arc].capa;
		model->add(constraint <= 0);
		constraint.end();

	}
	
	bool flag=false;
	for (int k = 0; k < ndemands; ++k) {
        itemd = &data->d_k[k];
		for (int i = 1; i <= nnodes; i++) {
            flag=false;
            IloExpr constraint(env);
            
            for(int a=sznz;a--;){
                arc = topo[a];
                item = &data->arcs[arc];
                if(i == item->i){
                    constraint -= x[k*sznz+a];
                    flag=true;
                }else if(i == item->j){
                    constraint += x[k*sznz+a];
                    flag=true;
                }
            }
            
            if( i == itemd->D){
                flag=true;
                constraint -= itemd->quantity;
            }else if( i == itemd->O){
                flag=true;
                constraint += itemd->quantity;
            }
            if(flag) model->add(constraint == 0);
            constraint.end();
		}
	}
	
	cplex.extract(*model);
}

//---------------------------------------------------------------------------

void MinCostFLow::set_parameters() {
	//cplex->setParam(IloCplex::Threads,0);
	//cplex->setParam(IloCplex::RootAlg, 2);
	//cplex->setParam(IloCplex::NodeAlg, 3);
	cplex.setParam(IloCplex::ClockType, 1);
	//cplex->setParam(IloCplex::MIPDisplay, 4);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds

}

//---------------------------------------------------------------------------

int MinCostFLow::solve(){
	
    //cplex.exportModel("rl.lp");
	cplex.solve();
	
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		cplex.clearModel();
		x.endElements();
		delete model;
		model=0;
        std::cout<<"infeasible"<<std::endl;
		return -1;
	}else if(cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
		cplex.clearModel();
		x.endElements();
		delete model;
		model=0;
        std::cout<<"InfeasibleOrUnbounded"<<std::endl;
		return -2;
		
	}else if(cplex.getStatus() == IloAlgorithm::Optimal){
		cplex.clearModel();
		x.endElements();
		delete model;
		model=0;
		return 0;
	}
	else return -3;
	
}

//---------------------------------------------------------------------------

MinCostFLow::~MinCostFLow() {

	try {
		cplex.clearModel();
		if(model)delete model;

		x.end();
		cplex.end();
		env.end();
		

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







