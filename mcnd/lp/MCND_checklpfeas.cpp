#include "MCND_checklpfeas.hpp"



#define INF 9999999
ILOSTLBEGIN

//---------------------------------------------------------------------------

void LPFeasChecker::set_parameters() {
	cplex.setParam(IloCplex::Threads,1);
    //cplex.setParam(IloCplex::RootAlg, 2);
    //cplex.setParam(IloCplex::NodeAlg, 3);
    cplex.setParam(IloCplex::ClockType, 1);
    //cplex.setParam(IloCplex::MIPDisplay, 4);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::DataCheck, CPX_DATACHECK_ASSIST);
    //cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds
    //cplex.setParam(IloCplex::Param::Simplex::Limits::Iterations	);
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
    addbd=0;
    cplex.extract(model);
}


//---------------------------------------------------------------------------

void 
LPFeasChecker::apply_bounds(const BCP_vec<BCP_var*>& vars){
	 for(int a=0;a<narcs;++a){
		if(vars[a]->ub()<=0.5){  continue; }
		for(int k=0;k<ndemands;++k){
			IloExpr constraint(env);
    		constraint += x[a*ndemands+k];
    		xbd.add((constraint<=vars[narcs+k*narcs+a]->ub()));
    		model.add(xbd[addbd++]);
    		constraint.end();
		}
	}
}

//---------------------------------------------------------------------------

void 
LPFeasChecker::apply_bounds(const double * colub){
	for(int a=0;a<narcs;++a){
		if(colub[a]<=0.5){  continue; }
		for(int k=0;k<ndemands;++k){
			IloExpr constraint(env);
    		constraint += x[a*ndemands+k];
    		xbd.add((constraint<=colub[narcs+k*narcs+a]));
    		model.add(xbd[addbd++]);
    		constraint.end();
		}
	}
	
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int
LPFeasChecker::solve_opt(int unfix, const BCP_vec<BCP_var*>& vars, const double * topo, 
						int * fixd, int& closed, MCND_solution*& mipsol, double fathmval, double betsolv) {
   reset();
   if(topo){
    	for(int a=0;a<narcs;++a){
    		if(vars[a]->ub()<=0.5 || topo[a]<=0.5){
    			fixd[closed++]=a;
    
    			IloExpr constraint(env);
    			for(int k=0;k<ndemands;++k) constraint += x[a*ndemands+k];
    			xbd.add((constraint<=0));
    			model.add(xbd[addbd++]);
    			constraint.end();
			}
		}
    }else{
		for(int a=0;a<narcs;++a){
			if(vars[a]->ub()<=0.5){ 
				fixd[closed++]=a;
				IloExpr constraint(env);
    			for(int k=0;k<ndemands;++k) constraint += x[a*ndemands+k];
    			xbd.add((constraint<=0));
    			model.add(xbd[addbd++]);
    			constraint.end();
			}
		}
    }
    cplex.setParam(IloCplex::Param::Advance, 0);
    cplex.setParam(IloCplex::RootAlg, 0);
    cplex.getObjective().setExpr(fobj2);
	int ret= solve();
	if(ret<0){ return -1;}
	
	bool bd_sat = true;
	double fthmv = fathmval;
	mipsol = new MCND_solution(narcs+narcs*ndemands);
	ret = getSolution(  vars, mipsol->xy, mipsol->cost, fthmv, bd_sat);
	if(topo || unfix==0){
		return 1;
	}
	if(fthmv >= betsolv){
		std::cout<<"test_feasibility::fthmv >= betsolv"<<std::endl;
		//delete mipsol;
		return 0;
	}else if(ret<0) return 2;
	if(bd_sat) return 1;
	
	apply_bounds(vars);
	cplex.setParam(IloCplex::Param::Advance, 1);
	cplex.setParam(IloCplex::RootAlg, 2);
	ret= solve();
	std:cout<<"test_feasibility::resolve"<<std::endl;
	if(ret<0){
		delete mipsol;
		return -2;
	} 
	fthmv = fathmval + compute_fathmval();
	if(fthmv >= betsolv){
		//delete mipsol;
		return 0;
	}
	

	return 1;	
}


//---------------------------------------------------------------------------

int LPFeasChecker::solve_feas(const double *collb, const double * colub, bool feas_status){
    reset();
    for(int a=0;a<narcs;++a){
        if(colub[a]<0.5){ 
        	IloExpr constraint(env);
    		for(int k=0;k<ndemands;++k) constraint += x[a*ndemands+k];
    		xbd.add((constraint<=0));
    		model.add(xbd[addbd++]);
			constraint.end();
        }
    }
	cplex.getObjective().setExpr(fobj1);
	cplex.setParam(IloCplex::Param::Advance, 0);
	cplex.setParam(IloCplex::RootAlg, 0);
 	int ret = solve();	 

	if(ret<0) return -2;
	if(feas_status==false) return -1;
	
	bool bd_sat=true;
	IloNumArray x_(env);
	cplex.getValues(x_,x);
	double val;
	for(int a=0;a<narcs;++a){
 		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k];
			if(val > colub[narcs+k*narcs+a]+1e-4){
				bd_sat=false;
				break;
				//std::cout<<"solve_feas::atencao var: "<<narcs+k*narcs+a<<" val: "<<val<<" ub: "<< colub[narcs+k*narcs+a]<<std::endl;
			}	
		}
	}
	x_.end();
	if(bd_sat) return 0;
	//cplex.exportModel("p.lp");
	apply_bounds(colub);
	cplex.setParam(IloCplex::Param::Advance, 1);
	cplex.setParam(IloCplex::RootAlg, 2);
	ret= solve();
	if(ret<0){
		return -3;
	}
	return 0;
}

//---------------------------------------------------------------------------

int 
LPFeasChecker::test_high_lowerbounds( BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, int * yfix,
										 const BCP_vec<BCP_var*>& vars,  const double* volsol, double bestsolv){
	reset();
	IloExpr obj(env);
	std::vector<int> unfix;
	for(int a=0;a<narcs;++a){
		if(vars[a]->ub()<=0.5 || (yfix[a]==0)){
			IloExpr constraint(env);
    		for(int k=0;k<ndemands;++k) constraint += x[a*ndemands+k];
    		xbd.add((constraint<=0));
    		model.add(xbd[addbd++]);
			constraint.end();
		} 
		else if(vars[a]->lb()>=0.5 || (yfix[a]==1)){
			for (int k = 0; k < ndemands; ++k){
				obj +=  data->arcs[a].c[k]*x[a*ndemands+k];
			}
			obj+=data->arcs[a].f;
		}else{
			if(volsol[a]<0.1 || volsol[a]>0.9)unfix.push_back(a);
			double cost = data->arcs[a].f/data->arcs[a].capa;
			for (int k = 0; k < ndemands; ++k){
				obj +=  (cost + data->arcs[a].c[k])*x[a*ndemands+k];
			}
		}	 
	} 
    
   	cplex.getObjective().setExpr(IloMinimize(env, obj));
   	obj.end();
   	apply_bounds(vars);
   	int ret=0;
   	cplex.setParam(IloCplex::Param::Advance, 0);
	//cplex.setParam(IloCplex::RootAlg, 2);
	
	if(!lbtested){ 
		ret= solve();
		if(ret<0){ std::cout<<"LPFeasChecker::test_high_lowerbounds: infeasible"<<std::endl; return -1;}
		std::cout<<"LPFeasChecker::test_high_lowerbounds: "<<cplex.getObjValue()<<" ub: "<<bestsolv<<std::endl;
		if(cplex.getObjValue()>bestsolv) return -1;
	}
	//cplex.setParam(IloCplex::Param::Advance, 0);

	int arc;
	int addctrnt =0;
	double additional=0;
	IloRangeArray xbdtemp(env);
	while(!unfix.empty()){

		arc = unfix.back();
		unfix.pop_back();
 		if(volsol[arc]>0.9){
 			IloExpr constraint(env);
    		for(int k=0;k<ndemands;++k) constraint += x[arc*ndemands+k];
    		xbdtemp.add((constraint<=0));
    		model.add(xbdtemp[addctrnt++]);
			constraint.end();
			 
		}else{
			 
			for(int k=0;k<ndemands;++k){
				cplex.getObjective().setLinearCoef(x[arc*ndemands+k], data->arcs[arc].c[k]); 
			}
			additional+=data->arcs[arc].f;
		}
 		
 		ret= solve();
		
 		if(ret>=0)std::cout<<"LPFeasChecker::test_high_lowerbounds "<<arc<<" : "<<cplex.getObjValue()+additional<<" ub: "<<bestsolv<<" "<<volsol[arc]<<std::endl;
 		else std::cout<<"LPFeasChecker::test_high_lowerbounds "<<arc<<": infeasible"<<std::endl;
 		
 		if( ret<0 || cplex.getObjValue()+additional > bestsolv){
			
			if(volsol[arc]>0.9){
				std::cout<<arc<<" LPFeasChecker::test_high_lowerbounds FIX 1 "<<std::endl;
				changed_pos.push_back(arc);
				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
				yfix[arc]=1;
				
 				for(int k=0;k<ndemands;++k){
					cplex.getObjective().setLinearCoef(x[arc*ndemands+k], data->arcs[arc].c[k]); 
				}
				additional+=data->arcs[arc].f;
			}else{
				std::cout<<arc<<" LPFeasChecker::test_high_lowerbounds FIX 0 "<<std::endl;
				changed_pos.push_back(arc);
				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				yfix[arc]=0;
				
				IloExpr constraint(env);
    			for(int k=0;k<ndemands;++k) constraint += x[arc*ndemands+k];
    			xbdtemp.add((constraint<=0));
    			model.add(xbdtemp[addctrnt++]);
				constraint.end();
			}
		}
		if(volsol[arc]>0.9){
			 xbdtemp[addctrnt-1].removeFromAll();
		}else{
			for(int k=0;k<ndemands;++k){
				cplex.getObjective().setLinearCoef(x[arc*ndemands+k], data->arcs[arc].f/data->arcs[arc].capa+data->arcs[arc].c[k]); 
			}
			additional-=data->arcs[arc].f;
		}
		
	}
	xbdtemp.endElements();
	return 0;
}


//---------------------------------------------------------------------------

int LPFeasChecker::solve(){
	//cplex.exportModel("rl.lp");
	cplex.solve();

	if(cplex.getStatus() == IloAlgorithm::Infeasible){
        //std::cout<<"infeasible"<<std::endl;
		return -2;
	}else if(cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded){		 
       // std::cout<<"InfeasibleOrUnbounded"<<std::endl;
		return -3;
		
	}else if(cplex.getStatus() == IloAlgorithm::Optimal){
		 //std::cout<<"optimal "<<cplex.getObjValue()<<std::endl;
		return 0;
	}
	else return -4;
}

//---------------------------------------------------------------------------


int 
LPFeasChecker::getSolution(const BCP_vec<BCP_var*>& vars, double * sol, double& solvalue, double &fathmval, bool& bd_sat){
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
			if(val> vars[narcs+k*narcs+a]->ub()+1e-4){
				bd_sat=false;
				//std::cout<<"atencao var: "<<narcs+k*narcs+a<<" val: "<<val<<" ub: "<< vars[narcs+k*narcs+a]->ub()<<std::endl;
			}	
		}
		if(flow>1e-4){
			sol[a] = 1.0;
			solvalue+=data->arcs[a].f;
			if(vars[a]->lb()<0.5 &&  vars[a]->ub()>0.5){ ++ret; }//std::cout<<"sola: "<<a<<" flow: "<<flow<<std::endl;}
		}//else if(vars[a]->lb()>=0.5){//std::cout<<"LPFeasChecker::getSolution FATHOM"<<std::endl; ret = -1;
		//else std::cout<<"solclosed: "<<a<<std::endl;
	}
	x_.end();
 	std::cout<<"LPFeasChecker::getSolution::feasibiliy integer sol value: "<<solvalue<<std::endl;
	return ret;
}

//---------------------------------------------------------------------------

double 
LPFeasChecker::compute_fathmval(){
 	double val;
	double cij;
	double fathmval=0;
	IloNumArray x_(env);
	cplex.getValues(x_,x);
	for(int a=0;a<narcs;++a){
 		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k];
 			cij = data->arcs[a].c[k]*val;
 			fathmval+=cij;

		}
	}
	x_.end();
 	//std::cout<<"feasibiliy integer sol value: "<<solvalue<<std::endl;
	return fathmval;
}

//---------------------------------------------------------------------------

void
LPFeasChecker::reset(){
	xbd.endElements();
	addbd=0;
}

//---------------------------------------------------------------------------

LPFeasChecker::~LPFeasChecker() {

	try {
		cplex.clearModel();
		xbd.endElements();
		x.endElements();
        fobj1.end();
        fobj2.end();
        
        model.end();
		cplex.end();
		env.end();
		

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







