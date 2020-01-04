#include "MCND_checkfeas.hpp"



#define INF 9999999
ILOSTLBEGIN



MinCostFeas::MinCostFeas(const Data * data) :
cplex(0),
model(0), env(0)
{
     nnodes = data->nnodes;
     ndemands = data->ndemands;
     narcs = data->narcs;

}

void
MinCostFeas::clear_model(){
    x.end();
    cplex->clearModel();
    cplex->end();
    env->end();
    
    delete env;
    delete cplex;
    delete model;
}



void MinCostFeas::create_model(const Topology & topo, const Data * data) {
	
		
	set_parameters();
	int sznz = topo.sznz;
	int arc;
	
    int nx=0;
    IloExpr obj(*env);
    obj+=0;
    x = IloNumVarArray(*env);
    idx.assign(narcs,0);
    
    for(int a=sznz;a--;){
        arc = topo.nz_arcs[a];
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(*env));
            idx[k*narcs+a]= nx++;
        }
    }
    
    
    //addicional FlowYibility vars
    model = new IloModel(*env);
    model->add(IloMinimize(*env, obj));
    obj.end();
    //constraints
    bool flag=false;
    for (int k = 0; k < ndemands; ++k) {
        for (int i = 1; i <= nnodes; i++) {
            
            flag=false;
            IloExpr constraint(*env);
            
            for(int a=sznz;a--;){
                arc = topo.nz_arcs[a];
                if(idx[k*narcs+arc]>=0){
                    if(i == data->fromi[a]){
                        constraint -= x[idx[k*narcs+arc]];
                        flag=true;
                    }else if(i == data->toj[a]){
                        constraint += x[idx[k*narcs+arc]];
                        flag=true;
                    }
                }
            }
            
            // std::cout<<"what1 "<<std::endl;
            
            
            if( i == data->D[k]){
                flag=true;
                constraint -=data->qnt[k];
            }
            if( i == data->O[k]){
                flag=true;
                constraint +=data->qnt[k];
            }
            
            if(flag){
                model->add((constraint == 0));
            }
            constraint.end();
            
        }
    }
    
    for(int a=sznz;a--;){
        arc = topo.nz_arcs[a];
        IloExpr constraint(*env);
        for (int k = 0; k < ndemands; ++k)
            if(idx[k*narcs+arc]>=0)constraint -= x[idx[k*narcs+arc]];
        constraint+=data->capa[a];
        
        model->add((constraint >= 0));
        constraint.end();
        
    }
    cplex->extract(*model);
}




void MinCostFeas::set_parameters() {
	//cplex->setParam(IloCplex::Threads,0);
	//cplex->setParam(IloCplex::RootAlg, 2);
	//cplex->setParam(IloCplex::NodeAlg, 3);
	cplex->setParam(IloCplex::ClockType, 1);
	//cplex->setParam(IloCplex::MIPDisplay, 4);
	cplex->setOut(env->getNullStream());
	cplex->setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds

}




int MinCostFeas::solve(){
	
	
	cplex->solve();
	
	if(cplex->getStatus() == IloAlgorithm::Infeasible){
		cplex->clearModel();
		delete model;
		model=0;
		return -1;
	}else if(cplex->getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
		cplex->clearModel();
		delete model;
		model=0;
		return -2;
		
	}else if(cplex->getStatus() == IloAlgorithm::Optimal){
		return 0;
	}
	else return -3;
	
}


MinCostFeas::~MinCostFeas() {
    idx.clear();
}







