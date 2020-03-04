#include "MCND_checklp.hpp"



#define INF 9999999
ILOSTLBEGIN

void
LPChecker::initialize(const Data* d) {
    data= d;
    ndemands = data->ndemands;
    nnodes = data->nnodes;
    narcs = data->narcs;
    cplex = new IloCplex(env);
    model = new IloModel(env);
    set_parameters();
}

//-------------------------------------------------------------------------------------------

void LPChecker::create_model(const std::vector<int> & topo, const CoverCollection & covers) {
	
    IloNumVarArray x(env, ndemands*narcs, 0.0, IloInfinity);
    IloNumVarArray y(env, narcs, 0.0, 1.0);
    IloExpr obj(env);
    for(int a=narcs;a--;){
        if(topo[a]==0) continue;
        else if(topo[a]==1)obj += data->arcs[a].f;
        else obj += data->arcs[a].f * y[a];
        for(int k = ndemands; k-- ;){
            obj += data->arcs[a].c[k]*x[k*narcs+a] ;
        }
        
    }
    
    model->add(IloMinimize(env, obj));
    obj.end();
    
    for (int k = ndemands; k--; ) {
        for (int i = nnodes; i-- ; ) {
            IloExpr constraint(env);
            for(int a =narcs ; a--; ){
                if(topo[a]==0) continue;
                if(i == data->arcs[a].i)
                    constraint += x[k*narcs+a];
                else if(i == data->arcs[a].j)
                    constraint -= x[k*narcs+a];
            }
            if( i == data->d_k[k].D)
                constraint += data->d_k[k].quantity;
            if( i == data->d_k[k].O)
                constraint -= data->d_k[k].quantity;
            model->add(constraint == 0);
            constraint.end();
        }
    }
    for(int a=narcs; a--; ){
        if(topo[a]==0) continue;
        IloExpr constraint(env);
        for (int k = ndemands; k-- ;) {
            constraint +=x[k*narcs+a];
        }
        if(topo[a]==1) constraint -= data->arcs[a].capa;
        else constraint -= data->arcs[a].capa * y[a];
        model->add(constraint <= 0);
        constraint.end();
    }
    for(int a=narcs; a--; ){
        if(topo[a]==0 || topo[a]==1) continue;
        for (int k = ndemands; k-- ;) {
            IloExpr constraint(env);
            constraint +=x[k*narcs+a];
            constraint -= data->arcs[a].b[k] * y[a];
            model->add(constraint <= 0);
            constraint.end();
        }
    }
    /*
    Cover * vi = covers.end;
    for(int i=covers.sizeOfCollection; i--;){
        IloExpr constraint(env);
        constraint -= vi->get_total_rhs();
        for(int a=vi->get_total_sz();a--;){
            if(topo[a]==0) continue;
            else if(topo[a]==1) constraint +=  vi->gamma_at(a);
            else constraint +=  vi->gamma_at(a)*y[vi->at(a)];
        }
        model->add(constraint >= 0);
        constraint.end();
        vi = vi->prev;
    }*/
    cplex->extract(*model);
    x.end();
    y.end();
}

//-------------------------------------------------------------------------------------------

void
LPChecker::set_parameters() {
    //cplex->setParam(IloCplex::Threads,1);
    //cplex->setParam(IloCplex::RootAlg, 1);
    //cplex->setParam(IloCplex::NodeAlg, 3);
    //cplex->setParam(IloCplex::ClockType, 1);
    //cplex->setParam(IloCplex::Param::MIP::Limits::EachCutLimit, 0);
    //cplex->setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
    //cplex->setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
    //cplex->setParam(IloCplex::MIPDisplay, 4);
    cplex->setOut(env.getNullStream());
    //cplex->setParam(IloCplex::TiLim, 10000.0); // Time limit in seconds
    
}

//-------------------------------------------------------------------------------------------

int
LPChecker::solve(){
	
	
	cplex->solve();
	
	if(cplex->getStatus() == IloAlgorithm::Infeasible){
		 
		return -1;
	}else if(cplex->getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
    
		return -2;
		
	}else if(cplex->getStatus() == IloAlgorithm::Optimal){
        std::cout<<"LPChecker: "<<cplex->getObjValue()<<std::endl;
		return 0;
	}
	else return -3;
	
}

//-------------------------------------------------------------------------------------------

LPChecker::~LPChecker() {
    
    try {
        cplex->clearModel();
        cplex->end();
        model->end();
        env.end();
        delete model;
        delete cplex;
    } catch (IloException& e) {
        std::cerr << "ERROR: " << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Error" << std::endl;
    }
}








