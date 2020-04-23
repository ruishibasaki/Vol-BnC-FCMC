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
	set_parameters();
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

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void  
Pump::reset( int szunfix_, std::deque<Pair2>& topo ){
	sznz = topo.size();
	szunfix = szunfix_;
	ytopo = topo;
	
	int arc;
	double norm=0;
	
	if(szunfix_==0) alpha= 1.0;
	else alpha= 0.5; 
	
	for(int a=szunfix;a--;){
		arc = topo[a].fst;
		norm+= pow(data->arcs[arc].f,2);
		//std::cout<<data->arcs[arc].f<<std::endl;
		for (int k = 0; k<ndemands;++k)
			norm+=pow(data->arcs[arc].c[k],2);
	}
	
	if(szunfix ==0 ){
		factorxy=1;
		factorp=0;
	}else{
		norm = sqrt(norm);
		factorxy = alpha/norm;
		factorp = (1.0 - alpha)/(sqrt(double(szunfix)));
	}
}

//-------------------------------------------------------------------------------------------

void Pump::create_model(const BCP_vec<BCP_var*>& vars) {	
	Pair2 item;
	int arc;
	double c;
	
	IloExpr obj(env);
	x = IloNumVarArray(env);
	y = IloNumVarArray(env);
	for(int a=0;a<szunfix;++a){
		y.add(IloNumVar(env,0,1));
		item = ytopo[a];
		arc = item.fst;
		if(item.snd == 1.0) c = factorxy*data->arcs[arc].f - factorp;
		else c = factorxy*data->arcs[arc].f  + factorp;
		obj += c*y[a];
		
		for (int k = 0; k < ndemands; ++k){
			x.add(IloNumVar(env, 0.0, vars[narcs+k*narcs+arc]->ub()));
			obj += factorxy*data->arcs[arc].c[k]*x[a*ndemands+k];
		}
	}

	for(int a=szunfix;a<sznz;++a){
		arc = ytopo[a].fst;
		for (int k = 0; k < ndemands; ++k){
			x.add(IloNumVar(env));
			obj += factorxy*data->arcs[arc].c[k]*x[a*ndemands+k];
		}
	}
	
	model = new IloModel(env);
	model->add(IloMinimize(env, obj));
	obj.end();
	
	//constraints
	bool flag;
	for (int k = 0; k < ndemands; ++k) {
		for (int i = 1; i <= nnodes; i++) {
			flag=false;
			IloExpr constraint(env);
		
			for(int a=0;a<sznz;++a){
				arc = ytopo[a].fst;
				if(i == data->arcs[arc].i){
					constraint += x[a*ndemands+k];
					flag = true;
				}else if(i == data->arcs[arc].j){
					constraint -= x[a*ndemands+k];
					flag = true;
				}
			}

		
			if( i == data->d_k[k].D){
				constraint +=data->d_k[k].quantity;
				flag = true;
			}else if( i == data->d_k[k].O){
				constraint -=data->d_k[k].quantity;
				flag = true;
			}
		
			if(flag){
				model->add((constraint == 0));			
			 }
			constraint.end();
		}
	}
	
	for(int a=0;a<sznz;++a){
		arc = ytopo[a].fst;
		IloExpr constraint(env);
		for (int k = 0; k < ndemands; ++k)
			constraint -= x[a*ndemands+k];
		
		if(a<szunfix)constraint+=data->arcs[arc].capa*y[a];
		else constraint+=data->arcs[arc].capa;

		model->add((constraint >= 0));
		constraint.end();
	}
	
	cplex.extract(*model);
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
		arc = ytopo[a].fst;
		for (int k = 0; k < ndemands; ++k){
			ub = vars[narcs+k*narcs+arc]->ub();
			if(x_[a*ndemands+k] - ub*y_[a]> 1e-10 ){
				IloExpr constraint(env);
				constraint -= x[a*ndemands+k];
				constraint+= ub*y[a];
				++cont;
				model->add((constraint >= 0));
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
Pump::solve(const BCP_vec<BCP_var*>& vars,  double * xy, double & val){
	
	create_model(vars);
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
		//std::cout<<"after cut "<<cplex.getObjValue()<<std::endl;
		cplex.getValues(y_,y);
		cplex.getValues(x_,x);
	}
	//std::cout<<"final pump "<<cplex.getObjValue()<<std::endl;

	val = getSolution(  xy,  x_, y_ );
	
	x_.end(); 
	y_.end();
	return 0;
	
}

//-------------------------------------------------------------------------------------------

double 
Pump::getSolution(  double * xy, const IloNumArray & x_, const IloNumArray & y_ ){
	
	int arc;
	double flow;
	double solvalue=0;
	double val;
	for(int a=0;a<sznz;++a){
		arc = ytopo[a].fst;
		flow = 0.0;
		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k];
			xy[narcs+k*narcs+arc] = val;
			solvalue+=data->arcs[arc].c[k]*val;
			flow += val;
		}
		if(a<szunfix){
			val = y_[a];
			if(val>1e-10){
				xy[arc] = 1.0;
 				solvalue+=data->arcs[arc].f;
			}
		}else if(flow>1e-10){
 			xy[arc] = 1.0;
			solvalue+=data->arcs[arc].f;
		} 
	}
	
	clear();
	//std::cout<<"sol value: "<<solvalue<<std::endl;
	return solvalue;
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void 
Pump::clear(){
	cplex.clearModel();
	x.endElements();
	y.endElements();
	if(model)delete model;
	model=0;
}


//-------------------------------------------------------------------------------------------

Pump::~Pump() {

	try {
		clear();
		cplex.end();
		env.end();
		
		ytopo.clear();

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







