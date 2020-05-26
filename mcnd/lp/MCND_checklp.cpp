//
//  MCND_checklp.cpp
//  
//
//  Created by Rui Shibasaki on 19/05/2020.
//

#include "MCND_checklp.hpp"

ILOSTLBEGIN

 
void LPChecker::set_parameters() {
	cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::RootAlg, 2);
    cplex.setParam(IloCplex::ClockType, 1);
    //cplex->setParam(IloCplex::MIPDisplay, 4);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::DataCheck, CPX_DATACHECK_ASSIST);
    cplex.setParam(IloCplex::TiLim, 3600.0); // Time limit in seconds
    //cplex.setParam(IloCplex::Param::Simplex::Limits::Iterations	);
}

//---------------------------------------------------------------------------

void
LPChecker::initialize(const Data* d, const CoverCollection* cover_man_,
    				const LocalCutCollection* localc_man_, const GlobalCutCollection* globalc_man_ ){
    data = d;
    covers = cover_man_;
    localcs = localc_man_;
    globalcs = globalc_man_;

    data = d;
    ndemands = d->ndemands;
    nnodes = d->nnodes;
    narcs =d->narcs;
    unfixd.resize(narcs,-1);
    bound_red.resize(narcs*ndemands,-1);
    
    szunfx=0;
    set_parameters();
    
    IloExpr obj(env);
     for(int a=0;a<narcs;++a){
    	y.add(IloNumVar(env, 0.0, 1.0));
    	obj += data->arcs[a].f*y[a];
        for (int k = 0; k < ndemands; ++k){
            x.add(IloNumVar(env, 0.0, IloInfinity));
            obj += data->d_k[k].quantity*data->arcs[a].c[k]*x[a*ndemands+k];
        }
    }
    fobj = IloMinimize(env, obj);
    model.add(fobj);
    obj.end();

    for(int a=0;a<narcs;++a){
        IloExpr constraint(env);
        for (int k = 0; k < ndemands; k++) {
            constraint += data->d_k[k].quantity*x[a*ndemands+k];
        }
        constraint -= data->arcs[a].capa*y[a];
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
                constraint -= 1.0;
            }else if( i ==  data->d_k[k].O){
                constraint += 1.0;
            }
            flowconserv.add((constraint == 0));
            model.add(flowconserv[k*nnodes+i-1]);
            constraint.end();
        }
    }
    
    cplex.extract(model);
}

//---------------------------------------------------------------------------

void 
LPChecker::reset(const BCP_vec<BCP_var*>& vars){
	szunfx=0;
	 for(int a=0;a<narcs;++a){
		if(vars[a]->ub()<0.5){  
			y[a].setUB(0.0);// std::cout<<"LPChecker::close: "<<a<<" "<<vars[a]->ub()<<std::endl;
			y[a].setLB(0.0);
		}else if(vars[a]->lb()>0.5){  
			y[a].setUB(1.0);
			y[a].setLB(1.0);
		}else{
			y[a].setUB(1.0);
			y[a].setLB(0.0);
			unfixd[szunfx++]=a;
		}
	}
	
 	map_addcovers.assign(covers->sizeOfCollection,-1);
    map_addlocals.assign(localcs->sizeOfCollection,-1);
    map_addglobals.assign(globalcs->sizeOfCollection,-1);
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void 
LPChecker::get_dual(){ 
	int fx = ndemands*nnodes ; 
	int sz = covers->sizeOfCollection+ localcs->sizeOfCollection+ globalcs->sizeOfCollection + fx;
	int id;
	double val;
	 
	
	IloNumArray pi(env);
	cplex.getDuals(pi, flowconserv);
	for (int k = 0; k < fx; ++k) 
		dual_map.insert(std::pair<int,double>(k, pi[k]));
    pi.end();
    
    if(numaddcov>0){
    	IloNumArray alph1(env);
    	//std::cout<<"numaddcov "<<numaddcov<<" "<<add_covers.getSize()<<std::endl;
		cplex.getDuals(alph1, add_covers);
		Cover * cov = covers->begin;
		sz = covers->sizeOfCollection;
		for(int i=0;i<sz;++i){
			id = map_addcovers[i];
			if(id>=0 ){
				val = alph1[id];
				dual_map.insert(std::pair<int,double>(fx+cov->serial_nmbr, val));
				//std::cout<<"cpxdual: "<<cov->serial_nmbr<<" "<<alph1[id]<<" "<<dual_map.find(fx+cov->serial_nmbr)->second<<std::endl;
			}
			cov = cov->next;
		}
		alph1.end();
    }
    
	if(numaddloc>0){
		//std::cout<<"numaddloc "<<numaddloc<<std::endl;
		IloNumArray alph2(env);
		cplex.getDuals(alph2, add_locals);
		LocalCut * loc = localcs->begin;
		sz = localcs->sizeOfCollection;
		for(int i=0;i<sz;++i){
			id = map_addlocals[i];
			if(id>=0 ){
				dual_map.insert(std::pair<int,double>(fx+loc->serial_nmbr, alph2[id]));
				//std::cout<<"cpxdual: "<<loc->serial_nmbr<<" "<<alph2[id]<<" "<<dual_map.find(fx+loc->serial_nmbr)->second<<std::endl;

			} 
			loc = loc->next;
		}
		alph2.end();
	}

	if(numaddgloc>0){
		//std::cout<<"numaddgloc "<<numaddgloc<<std::endl;
		IloNumArray alph3(env);
		cplex.getDuals(alph3, add_globals);
		GlobalCut * gloc = globalcs->begin;
		sz = globalcs->sizeOfCollection;
		for(int i=0;i<sz;++i){
			id = map_addglobals[i];
			if(id>=0 ){
				dual_map.insert(std::pair<int,double>(fx+gloc->serial_nmbr, alph3[id]));
			} 
			gloc = gloc->next;
		}
		alph3.end();
	}
}

//---------------------------------------------------------------------------

void 
LPChecker::get_solution(MCND_solution* mipsol){
	double flow;
	double val;
	double cij;
	int ret=-1; 
	mipsol->cost=0;
	 
	for(int a=0;a<narcs;++a){
 		flow = 0.0;
		for (int k = 0; k < ndemands; ++k){
			val = x_[a*ndemands+k]*data->d_k[k].quantity;
			mipsol->xy[narcs+k*narcs+a] = val;
			mipsol->cost += data->arcs[a].c[k]*val;
			flow += val;
		}
		if(flow>1e-4){
			mipsol->xy[a] = 1.0;
			mipsol->cost+=data->arcs[a].f;
		}
	}
 	std::cout<<"LPChecker::getSolution:: integer sol value: "<<mipsol->cost<<std::endl;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


int LPChecker::solve(double ub, const BCP_vec<BCP_var*>& vars,  double& new_lb){
	
	clean();
	reset(vars);
	std::cout<<"LPChecker::solve"<<std::endl;
	cplex.solve();
	cplex.setParam(IloCplex::Param::Advance, 1);

	//std::cout<<"LPChecker first "<<cplex.getObjValue()<<std::endl;
	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		cplex.exportModel("rl.lp");
		std::cout<<"LPChecker::solve:: cplex.getStatus() == IloAlgorithm::Infeasible 1"<<std::endl;
  		return -1;
	} 
	
	cplex.getValues(x_,x);
	cplex.getValues(y_,y);
	while(cut(vars, y_,x_)){
		cplex.solve();
		if(cplex.getStatus() == IloAlgorithm::Infeasible){ 
			std::cout<<"LPChecker::solve:: cplex.getStatus() == IloAlgorithm::Infeasible 2"<<std::endl; 
			return -2; 
		}
		//std::cout<<"after cut "<<cplex.getObjValue()<<std::endl;
		if(cplex.getObjValue()+1e-4>= ub){
			new_lb = cplex.getObjValue();
			std::cout<<"stop lpchecker "<<cplex.getObjValue()<<std::endl;
			return 1;
		}
		cplex.getValues(y_,y);
		cplex.getValues(x_,x);
	}
	solval = cplex.getObjValue();
	std::cout<<"final lpchecker "<<solval<<std::endl;
	get_dual();
	new_lb = solval;
	
	return 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


int 
LPChecker::cut(const BCP_vec<BCP_var*>& vars, const IloNumArray & y_, const IloNumArray & x_){
	int arc, sz;
	int cont=0;
	double ub;
	for(int a=0;a<szunfx;++a){
		arc = unfixd[a];
		for (int k = 0; k < ndemands; ++k){
			ub = vars[narcs+k*narcs+arc]->ub();
			if(data->d_k[k].quantity*x_[arc*ndemands+k] - ub*y_[arc]> 1e-4 ){
				IloExpr constraint(env);
				constraint += x[arc*ndemands+k];
				constraint -= (ub/data->d_k[k].quantity)*y[arc];
 				//std::cout<<"add: "<<arc<<" "<<k<<" viol: "<< data->d_k[k].quantity*x_[arc*ndemands+k] - ub*y_[arc]<<std::endl;
 				add_strong_force.add((constraint <= 0));
				model.add(add_strong_force[numaddstrong++]);
				constraint.end();
				++cont;
			}
		}
	}
 	Cover * cov = covers->begin;
	sz = covers->sizeOfCollection;
	for(int i=0;i<sz;++i){
 		if(map_addcovers[i]<0 && coverc_viol(y_ , cov) ){
			//std::cout<<"LPChecker::inserst cover "<<cov->serial_nmbr<<std::endl;
			IloExpr constraint(env);
			for(int a=cov->get_total_sz(); a--; )
				constraint += cov->gamma_at(a)*y[ cov->at(a) ];
			constraint -= cov->get_rhs();
	
			add_covers.add((constraint >= 0));
			map_addcovers[i] = numaddcov;
			model.add(add_covers[numaddcov++]);
			constraint.end();
			++cont;
		}
		cov = cov->next;
	}
  	LocalCut * loc = localcs->begin;
	sz = localcs->sizeOfCollection;
 	for(int i=0;i<sz;++i){
		if(map_addlocals[i]<0 && localc_viol(y_ , loc) ){
			//std::cout<<"LPChecker::inserst local"<<loc->serial_nmbr<<std::endl;
			IloExpr constraint(env);
			for(int a=loc->size; a--; ){
				 constraint += loc->coef_at(a) * y[ loc->vars[a] ];
			}
			constraint -= loc->rhs;
			
			if(loc->sense==1)
				add_locals.add((constraint >= 0));
			else add_locals.add((constraint <= 0));
			map_addlocals[i] = numaddloc;
			model.add(add_locals[numaddloc++]);
			constraint.end();
			++cont;
		}
		loc = loc->next;
	}
	GlobalCut * gloc = globalcs->begin;
	sz = globalcs->sizeOfCollection;
 	for(int i=0;i<sz;++i){
		if(map_addglobals[i]<0 && globalc_viol(y_ , gloc) ){
			//std::cout<<"LPChecker::inserst global"<<std::endl;
			IloExpr constraint(env);
			for(int a=gloc->size; a--; ){
				 constraint += y[ gloc->vars[a] ];
			}
			constraint -= 1.0;
			
			add_globals.add((constraint >= 0));
			map_addglobals[i] = numaddgloc;
			model.add(add_globals[numaddgloc++]);
			constraint.end();
			++cont;
		}
		gloc = gloc->next;
	}
 	return cont;
}

//---------------------------------------------------------------------------

bool 
LPChecker::globalc_viol(const IloNumArray & y_ , const GlobalCut* gloc) {
 	int arc=0;
 	double sum=0;
    for(int a=gloc->size ;a--; ){
      	arc = gloc->vars[a];
      	sum += y_[arc];
      	if(sum >= 1.0)
      		 return false;
    }
    return true;
}

//---------------------------------------------------------------------------

bool 
LPChecker::coverc_viol(const IloNumArray & y_ , const Cover* cover) {
 	int arc=0;
 	double sum=0;
 	double rhs = cover->get_rhs();
 	int sz = cover->get_total_sz();
    for(int a=sz ;a--; ){
      	sum+= cover->gamma_at(a)*y_[ cover->at(a) ];
        if(sum>=rhs){return false;}   
    }
    return true;
}

//---------------------------------------------------------------------------

bool 
LPChecker::localc_viol(const IloNumArray & y_ , const LocalCut* loc) {
	int arc;
 	double sum=0; 
	double rhs_ = loc->rhs;
	int size = loc->size;

    for(int a=0;a<size;++a){
    	sum+= loc->coef_at(a) * y_[ loc->vars[a] ];
        if(sum>=rhs_ && loc->sense==1){return false;} 
        else if(sum>=rhs_ && loc->sense==-1){return true;}
    }
	if(loc->sense==1)
   		return true;
   	else return false;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void 
LPChecker::logical_yfix(double ub, BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, int * yfix){
	int arc, id;
	double rc, ysol;
	
	cplex.getReducedCosts(rc_y, y);

	for(int a=0;a<szunfx;++a){
		arc = unfixd[a];
		if(yfix[arc] >=0) continue;
		
		rc = rc_y[arc];
		if(solval + abs(rc) >= ub){
			if(ysol >=1.0-1e-10){
				changed_pos.push_back(arc);
 				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
				yfix[arc] = 1.0;
				std::cout<<arc<<" LOGFIX 1 (LP optimal)"<<std::endl;
			}else if(ysol <=1e-10){
				changed_pos.push_back(arc);
 				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
				yfix[arc] = 0.0;
				std::cout<<arc<<" LOGFIX 0 (LP optimal)"<<std::endl;
			} 
		} 
		
	}

}

//---------------------------------------------------------------------------

void 
LPChecker::logical_xfix(double ub, const int * yfix, const BCP_vec<BCP_var*>& vars){
	int arc, id;
	double tt;
	double xrc, dk;
 	cplex.getReducedCosts(rc_x, x);

	for(int a=0;a<szunfx;++a){
		arc = unfixd[a];
		if(yfix[arc] ==0) continue;
		
		for (int k = 0; k < ndemands; k++) {
			id = arc*ndemands+k;
			xrc = rc_x[id];
			
			if(x_[id]<=0 && xrc>1e-10){
				dk = data->d_k[k].quantity;
				tt = solval +rc_y[arc]*(1.0 -  y_[arc]);
				
				if(tt + xrc*(vars[narcs+k*narcs+arc]->ub()/dk) > ub +1e-4){
					bound_red[id] = dk*(ub - tt)/rc_x[id];
					
				}else bound_red[id] = -1;
			}else bound_red[id] = -1;
       	}
       	
	}
	
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void
LPChecker::clean(){
	map_addcovers.clear();
    map_addlocals.clear();
   	map_addglobals.clear();
    numaddstrong = numaddcov= numaddloc= numaddgloc =0;
	add_strong_force.endElements();
	add_covers.endElements();
	add_locals.endElements();
	add_globals.endElements();
	cplex.setParam(IloCplex::Param::Advance, 0);
	dual_map.clear();
}
	

//---------------------------------------------------------------------------


LPChecker::~LPChecker() {

	try {
		cplex.clearModel();

		x.endElements();
		y.endElements();
		x_.end();
		y_.end();
		rc_x.end();
		rc_y.end();
		flowconserv.endElements();
    
		add_strong_force.endElements();
		add_covers.endElements();
		add_locals.endElements();
		add_globals.endElements();
	
		unfixd.clear();
    	map_addcovers.clear();
    	map_addlocals.clear();
   		map_addglobals.clear();
		
        fobj.end();
        
        model.end();
		cplex.end();
		env.end();
 
	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}
