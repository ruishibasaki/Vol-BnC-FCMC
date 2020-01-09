#include "MCND_pump.hpp"



#define INF 9999999
ILOSTLBEGIN



Pump::Pump(const Topology * t, const Data * d) :
	cplex(env),
	x(env), y(env),
	flow_row(env), capa_row(env), strong_row(env),
	model(0){
	
	topo = t;
	data = d;
	nnodes = data->nnodes;
	ndemands = data->ndemands;
	narcs = data->narcs;
	sznz = topo->sznz;
	
	yint.resize(narcs,0);
	int arc;
	for(int a=0;a<sznz;++a){
		arc = topo->nz_arcs[a];
		yint[arc] = 1;
	}
	
	double normf=0;
	factory=narcs;
	factorx=0;
	for(int a=0;a<narcs;++a){
		normf+= pow(data->fixc[a],2);
		for (int k = 0; k<ndemands;++k)
			factorx+=pow(data->cost[a][k],2);
	}
	factory /= (normf);
	factory = sqrt(factory);
	factorx = 1.0/(factorx);
	factorx = sqrt(factorx);
}


void Pump::create_model(int numcols, const int * init_cols, double alpha, bool feas) {
	
	
	
	set_parameters();
	idy.resize(narcs,-1);
	idx.resize(narcs*ndemands,-1);
	
	double c;
	nx=0; ny=0;
	factory *= alpha;
	factorx *= alpha;
	
	IloExpr obj(env);
	for(int a=0;a<narcs;++a){
		y.add(IloNumVar(env,0,1));
		idy[a] = ny++;
		
		if(yint[a]) c = factory*data->fixc[a]  + alpha -1;
		else c = factory*data->fixc[a]  +1 - alpha;
		obj += c*y[a];
	}
	
	//addicional feasibility vars
	feas_init = feas;
	if(!feas){
		for (int k = 0; k < ndemands; ++k){
			x.add(IloNumVar(env));
			obj += 1e10 *x[nx];
			++nx;
		}
	}
	if(numcols){
		int a,k;
		for (int c = 0; c < numcols; ++c){
			k = init_cols[c]/ndemands;
			a = init_cols[c]%ndemands;
			x.add(IloNumVar(env));
			obj += factorx*data->cost[a][k]*x[nx];
			idx[init_cols[c]]= nx++;
		}
	}
	
	
	
	if(model){
		cplex.clearModel();
		delete model;
	}
	model = new IloModel(env);
	Obj = IloMinimize(env, obj);
	model->add(Obj);
	obj.end();
	
	//constraints
	idf_row.resize(nnodes*ndemands,-1);
	idc_row.resize(narcs,-1);
	ids_row.resize(narcs*ndemands,-1);
	nfr=0; ncr=0; nsr=0;
	bool flag=false;
	for (int k = 0; k < ndemands; ++k) {
		for (int i = 1; i <= nnodes; i++) {
			    flag=false;
				IloExpr constraint(env);
				
				for(int a=0;a<narcs;++a){
										
					if(i == data->fromi[a]){
						if(idx[k*narcs+a]>=0)
							constraint -= x[idx[k*narcs+a]];
						flag=true;
					}else if(i == data->toj[a]){
						if(idx[k*narcs+a]>=0)
							constraint += x[idx[k*narcs+a]];
						flag=true;
					}
				}
	
				
				if( i == data->D[k]){
					flag=true;
					if(!feas)constraint += x[k];
					constraint -=data->qnt[k];
				}if( i == data->O[k]){
					flag=true;
					if(!feas)constraint -= x[k];
					constraint +=data->qnt[k];
				}
				
				if(flag){
					flow_row.add((constraint == 0));
					model->add(flow_row[nfr]);
					idf_row[k*nnodes+i-1]= nfr++;
					
				 }
				constraint.end();
		}
	}
	
	for(int a=0;a<narcs;++a){
		IloExpr constraint(env);
		for (int k = 0; k < ndemands; ++k)
			if(idx[k*narcs+a]>=0)
				constraint -= x[idx[k*narcs+a]];
		constraint+=data->capa[a]*y[a];
		
		capa_row.add((constraint >= 0));
		model->add(capa_row[ncr]);
		idc_row[a]= ncr++;
		constraint.end();
	}
	
	cplex.extract(*model);
}


int 
Pump::price(std::deque<int>& cols_to_add, const IloNumArray & y_,const IloNumArray & pi_,const IloNumArray & alpha_ ){
	double rc, rcf, totalrc;
	int i,j;
	std::deque<int> col_cand;
	for(int a=0;a<narcs;++a){
		i = data->fromi[a]-1;
		j = data->toj[a]-1;
		totalrc=0;
		//std::cout<<" arc: "<<arc<<std::endl;
		for (int k = 0; k < ndemands; ++k){
			if(idx[k*narcs+a]<0){
				rc = factorx*data->cost[a][k] + alpha_[idc_row[a]];
				rc += pi_[idf_row[k*nnodes+i]]- pi_[idf_row[k*nnodes+j]];
				if(y_[a]>1e-30 && rc<0){
					cols_to_add.push_back(k*narcs+a);
					//std::cout<<"add col k: "<<k<<" id: "<<cols_to_add.back()<<std::endl;
				}else if(rc<0){
					totalrc += -data->qnt[k]*rc;
					col_cand.push_back(k*narcs+a);
					//std::cout<<"candidate k: "<<k<<" id: "<<col_cand.back()<<std::endl;
				}
			}
		}
		if(y_[a]<=1e-30 && !col_cand.empty() ){
			rcf= factory*data->fixc[a] - data->capa[a]*alpha_[idc_row[a]];
			if(rcf<totalrc){
				//std::cout<<"insert candidates"<<std::endl;
				cols_to_add.insert(cols_to_add.end(),col_cand.begin(), col_cand.end());
			}
		}
		col_cand.clear();
	}
	
	//std::cout<<"num cols to add: "<<cols_to_add.size()<<std::endl;
   //std::cout<<"cols to add: "<<cols_to_add.size()<<std::endl;
	//for(int c=cols_to_add.size();c--;)
			//std::cout<<" id: "<<cols_to_add[c]<<std::endl;
	return cols_to_add.size();
}

void 
Pump::add_cols(std::deque<int>& cols_to_add){
	int k,a,i,j;
	//std::cout<<"num cols to add: "<<cols_to_add.size()<<std::endl;
	for(int c=cols_to_add.size();c--;){
		k = cols_to_add[c]/narcs;
		a = cols_to_add[c]%narcs;
		i = data->fromi[a]-1;
		j = data->toj[a]-1;
		//std::cout<<"add col k: "<<k<<" arc: "<<a<<" id: "<<cols_to_add[c]<<std::endl;
		IloNumColumn col =capa_row[idc_row[a]](-1);
		col+=flow_row[idf_row[k*nnodes+i]](-1);
		col+=flow_row[idf_row[k*nnodes+j]](1);
		col+= Obj(factorx*data->cost[a][k]);
		
		x.add(IloNumVar(col));
		idx[cols_to_add[c]]= nx++;
		col.end();
	}
	cols_to_add.clear();
}

void 
Pump::cut(const IloNumArray & y_, const IloNumArray & x_){
	int index;
	for(int a=0;a<narcs;++a){
		for (int k = 0; k < ndemands; ++k){
			index = idx[k*narcs+a];
			if(index>=0){
				if(x_[index]>data->b[a][k]*y_[a]){
					IloExpr constraint(env);
					constraint -= x[index];
					constraint+=data->b[a][k]*y[a];
					
					strong_row.add((constraint >= 0));
					model->add(strong_row[nsr]);
					ids_row[k*narcs+a]= nsr++;
					constraint.end();
				}
			}
		}
	}
}

void Pump::set_parameters() {
	//cplex->setParam(IloCplex::Threads,0);
	//cplex->setParam(IloCplex::RootAlg, 2);
	//cplex->setParam(IloCplex::NodeAlg, 3);
	cplex.setParam(IloCplex::ClockType, 1);
	//cplex->setParam(IloCplex::MIPDisplay, 4);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds

}




int Pump::solve(){
	
	
	cplex.solve();
	IloNumArray x_(env);
	IloNumArray y_(env);
	IloNumArray pi_(env);
	IloNumArray alpha_(env);
	
	cplex.getValues(y_,y);
	cplex.getDuals(pi_,flow_row);
	cplex.getDuals(alpha_,capa_row);
	
	std::deque<int> cols_to_add;
	while(price(cols_to_add, y_, pi_, alpha_)){
		add_cols(cols_to_add);
		cplex.solve();
		cplex.getValues(x_,x);
		cplex.getValues(y_,y);
		std::cout<<"after price "<<cplex.getObjValue()<<std::endl;
		cut(y_,x_);
		//std::cout<<"done cut"<<std::endl;
		cplex.solve();
		std::cout<<"after cut "<<cplex.getObjValue()<<std::endl;
		cplex.getValues(y_,y);
		cplex.getDuals(pi_,flow_row);
		cplex.getDuals(alpha_,capa_row);
	}
	cplex.getValues(x_,x);
	cut(y_,x_);
	std::cout<<"done cut"<<std::endl;
	cplex.solve();
	std::cout<<"final "<<cplex.getObjValue()<<std::endl;
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		cplex.clearModel();
		delete model;
		model=0;
		return -1;
	}else if(cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){
		cplex.clearModel();
		delete model;
		model=0;
		return -2;
		
	}else if(cplex.getStatus() == IloAlgorithm::Optimal){
		return 0;
	}
	else return -3;
	
}

double 
Pump::getSolution(double * xy){
	
	double solvalue = 0;
		
	IloNumArray xsol(env);
	cplex.getValues(x, xsol);
	bool flag;
	
	if(!feas_init)
		for (int k = 0; k < topo->ndemands; ++k)
			if(xsol[k]>0) return -1;
			
	for(int a=0;a<narcs;++a){
		flag=false;
		//std::cout<<"unfxd: "<<arc<<" / "<<a<<std::endl;
		for (int k = 0; k < topo->ndemands; ++k){
			if(idx[k*narcs+a]>=0){
				xy[narcs+k*narcs+a] = xsol[idx[k*narcs+a]];
				solvalue+=data->cost[a][k]*xy[narcs+k*narcs+a];
				//std::cout<<solvalue<<" var "<<idx[k*narcs+arc]<<": "<<xy[narcs+k*narcs+arc]<<std::endl;
			}else xy[narcs+k*narcs+a] =0;
			
			if(xy[narcs+k*narcs+a]>1e-30)flag=true;
			
		}
		if(flag){
			xy[a] = 1.0;
			solvalue += data->fixc[a];
			//std::cout<<solvalue<<" arc "<<arc<<std::endl;
		}
	}
	
	cplex.clearModel();
	delete model;
	model=0;
	return solvalue;
}

Pump::~Pump() {

	try {
		cplex.clearModel();
		delete model;

		x.end();
		y.end();
		Obj.end();
		flow_row.end();
		capa_row.end();
		strong_row.end(); 
		cplex.end();
		env.end();
		idy.clear();
		idx.clear();
		idf_row.clear();
		idc_row.clear();
		ids_row.clear();
		yint.clear();

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







