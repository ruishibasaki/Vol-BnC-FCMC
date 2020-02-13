#include "MCND_pump.hpp"



#define INF 9999999
ILOSTLBEGIN



Pump::Pump(const Data * d) :
	cplex(env),
	x(env), y(env),
	model(0){
	
	data = d;
	nnodes = data->nnodes;
	ndemands = data->ndemands;
	narcs = data->narcs;
	set_parameters();
}

//-------------------------------------------------------------------------------------------

void Pump::set_parameters() {
	//cplex->setParam(IloCplex::Threads,0);
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
Pump::reset(double alph, int szunfix_, std::deque<Pair2>& topo ){
	sznz = topo.size();
	szunfix = szunfix_;
	ytopo = topo;
	alpha = alph;
	
	int arc;
	double normf=0;
	double normx=0;
	for(int a=szunfix;a--;){
		arc = topo[a].fst;
		normf+= pow(data->arcs[arc].f,2);
		for (int k = 0; k<ndemands;++k)
			normx+=pow(data->arcs[arc].c[k],2);
	}
	factory = alpha/normf;
	factory = sqrt(factory);
	factorx = alpha/normx;
	factorx = sqrt(factorx);
}

//-------------------------------------------------------------------------------------------

void Pump::create_model() {	
	Pair2 item;
	int arc;
	double c;
	double factorp = (1.0 - alpha)/(sqrt(double(szunfix)));
	
	IloExpr obj(env);
	x = IloNumVarArray(env);
	y = IloNumVarArray(env);
	for(int a=0;a<szunfix;++a){
		y.add(IloNumVar(env,0,1));
		item = ytopo[a];
		arc = item.fst;
		if(item.snd == 1.0) c = factory*data->arcs[arc].f - factorp;
		else c = factory*data->arcs[arc].f  + factorp;
		obj += c*y[a];
		
		for (int k = 0; k < ndemands; ++k){
			x.add(IloNumVar(env));
			obj += factorx*data->arcs[arc].c[k]*x[a*ndemands+k];
		}
	}
	v0=0;
	for(int a=szunfix;a<sznz;++a){
		arc = ytopo[a].fst;
		v0 += data->arcs[arc].f;
		for (int k = 0; k < ndemands; ++k){
			x.add(IloNumVar(env));
			obj += factorx*data->arcs[arc].c[k]*x[a*ndemands+k];
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
					constraint -= x[a*ndemands+k];
					flag = true;
				}else if(i == data->arcs[arc].j){
					constraint += x[a*ndemands+k];
					flag = true;
				}
			}

		
			if( i == data->d_k[k].D){
				constraint -=data->d_k[k].quantity;
				flag = true;
			}if( i == data->d_k[k].O){
				constraint +=data->d_k[k].quantity;
				flag = true;
			}
		
			if(flag){
				model->add((constraint == 0));			
			 }
			constraint.end();
		}
	}
	
	for(int a=0;a<sznz;++a){
		IloExpr constraint(env);
		for (int k = 0; k < ndemands; ++k)
			constraint -= x[a*ndemands+k];
		
		if(a<szunfix)constraint+=data->arcs[a].capa*y[a];
		else constraint+=data->arcs[a].capa;

		model->add((constraint >= 0));
		constraint.end();
	}
	
	cplex.extract(*model);
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


int 
Pump::cut(const IloNumArray & y_, const IloNumArray & x_){
	int arc;
	int cont=0;
	for(int a=0;a<szunfix;++a){
		arc = ytopo[a].fst;
		for (int k = 0; k < ndemands; ++k){
			if(x_[a*ndemands+k] > data->arcs[arc].b[k]*y_[a]){
				IloExpr constraint(env);
				constraint -= x[a*ndemands+k];
				constraint+= data->arcs[arc].b[k]*y[a];
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
Pump::solve(double * yl, double * xy, double & val){
	
	create_model();
	cplex.solve();
	
	if(cplex.getStatus() == IloAlgorithm::Infeasible) return -1;
	else if(cplex.getStatus() == IloAlgorithm::Unbounded)return -2;
	
	IloNumArray x_(env);
	IloNumArray y_(env);

	cplex.getValues(x_,x);
	cplex.getValues(y_,y);
	while(cut(y_,x_)){
		cplex.solve();
		//std::cout<<"after cut "<<cplex.getObjValue()<<std::endl;
		cplex.getValues(y_,y);
		cplex.getValues(x_,x);
	}
	std::cout<<"final pump "<<cplex.getObjValue()<<std::endl;

	val = getSolution( yl, xy,  x_, y_ );
	
	x_.end(); 
	y_.end();
	return 0;
	
}

//-------------------------------------------------------------------------------------------

double 
Pump::getSolution(double * yl, double * xy, const IloNumArray & x_, const IloNumArray & y_ ){
	
	int arc;
	double solvalue = v0;
	double val;
	for(int a=0;a<sznz;++a){
		arc = ytopo[a].fst;
		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k];
			xy[narcs+k*narcs+arc] = val;
			solvalue+=data->arcs[arc].c[k]*val;
		}
		val = y_[a];
		yl[arc] = val;
		if(val>1e-10){
			xy[arc] = 1.0;
			solvalue+=data->arcs[arc].f*val;
		}
	}
	
	clear();
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







