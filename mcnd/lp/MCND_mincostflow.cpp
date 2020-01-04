#include "MCND_mincostflow.hpp"



#define INF 9999999
ILOSTLBEGIN


//---------------------------------------------------------------------------

MinCostFLow::MinCostFLow() :
	cplex(env),
	model(0)
{
	
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
	
    IloExpr obj(env);
    
    set_parameters();
    
    model = new IloModel(env);
    x = IloNumVarArray(env, sznz*ndemands);
	for(int a=sznz;a--;){
		for (int k = 0; k < ndemands; ++k)
			x[k*sznz+a] = IloNumVar(env);
	}
    
	for(int a=sznz;a--;){
		arc =  topo[a];
        item = &data->arcs[arc];
        obj += item->f;
        
		IloExpr constraint(env);
		for (int k = 0; k < ndemands; k++) {
			constraint += x[k*sznz+a];
			obj += item->c[k] *x[k*sznz+a];
		}
		
		constraint -= item->capa;
		model->add(constraint <= 0);
		constraint.end();

	}
	model->add(IloMinimize(env, obj));
	obj.end();
	
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
		delete model;
		model=0;
        std::cout<<"infeasible"<<std::endl;
		return -1;
	}else if(cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
		cplex.clearModel();
		delete model;
		model=0;
        std::cout<<"InfeasibleOrUnbounded"<<std::endl;
		return -2;
		
	}else if(cplex.getStatus() == IloAlgorithm::Optimal){
		return 0;
	}
	else return -3;
	
}

//---------------------------------------------------------------------------

double 
MinCostFLow::getSolution(const std::deque<int> & topo, const Data * data, MCND_solution * sol_){
	
	sol_->cost_flow = cplex.getObjValue();
	sol_->cost =sol_->cost_flow ;
	
    int sznz = topo.size();
	int arc;
	IloNumArray xsol(env);
    cplex.getValues(x, xsol);

	bool flag;
	for(int a=sznz;a--;){
		//flag=false;
		arc =  topo[a];
		//std::cout<<"unfxd: "<<arc<<" / "<<a<<std::endl;
		for (int k = 0; k < data->ndemands; ++k){
			sol_->xy[data->narcs+k*data->narcs+arc] = xsol[k*sznz+a];
			//if(xsol[k*sznz+a])flag=true;
		}
        sol_->xy[arc] = 1.0;
		/*if(flag){
			sol_->xy[arc] = 1.0;
			sol_->cost += data->fixc[arc];
		}else sol_->xy[arc] = 0.0;*/
	}
	
	cplex.clearModel();
	delete model;
	model=0;
	return sol_->cost;
}

//---------------------------------------------------------------------------

MinCostFLow::~MinCostFLow() {

	try {
		cplex.clearModel();
		delete model;

		x.end();
		cplex.end();
		env.end();
		

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







