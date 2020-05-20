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
    //cplex->setParam(IloCplex::NodeAlg, 3);
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
			y[a].setUB(0.0);
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
	
	/*add_strong_force = IloRangeArray(env);
	add_covers= IloRangeArray(env);
	add_locals= IloRangeArray(env);
	add_globals= IloRangeArray(env);*/
	map_addcovers.assign(covers->sizeOfCollection,-1);
    map_addlocals.assign(localcs->sizeOfCollection,-1);
    map_addglobals.assign(globalcs->sizeOfCollection,-1);
    numaddstrong = numaddcov= numaddloc= numaddgloc =0;
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
			if(data->d_k[k].quantity*x_[arc*ndemands+k] - ub*y_[arc]> 1e-10 ){
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
 		if(coverc_viol(y_ , cov) ){
			std::cout<<"LPChecker::inserst cover"<<std::endl;
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
		if(localc_viol(y_ , loc) ){
			std::cout<<"LPChecker::inserst local"<<std::endl;
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
		if(globalc_viol(y_ , gloc) ){
			std::cout<<"LPChecker::inserst global"<<std::endl;
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
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


int LPChecker::solve(double ub, const BCP_vec<BCP_var*>& vars,  double& new_lb){
	
	//cplex.exportModel("rl.lp");
	reset(vars);
	std::cout<<"LPChecker::solve"<<std::endl;
	cplex.solve();
	std::cout<<"LPChecker first "<<cplex.getObjValue()<<std::endl;

	if(cplex.getStatus() == IloAlgorithm::Infeasible){
		std::cout<<"LPChecker::solve:: cplex.getStatus() == IloAlgorithm::Infeasible 1"<<std::endl;
  		return -1;
	} 
	IloNumArray x_(env);
	IloNumArray y_(env);
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
	std::cout<<"final lpchecker "<<cplex.getObjValue()<<std::endl;
	get_dual();
	new_lb = cplex.getObjValue();
	
	return 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void 
LPChecker::get_dual(){ 
	int fx = ndemands*nnodes ; 
	int sz = covers->sizeOfCollection+ localcs->sizeOfCollection+ globalcs->sizeOfCollection + fx;
	int id;
	if(dual_sol) delete [] dual_sol;
	dual_sol = new double [sz];
	
	IloNumArray pi(env);
	cplex.getDuals(pi, flowconserv);
	for (int k = 0; k < fx; ++k) 
       	dual_sol[k] = pi[k];
    pi.end();
    
    IloNumArray alph1(env);
	cplex.getDuals(alph1, add_covers);
    Cover * cov = covers->begin;
    sz = covers->sizeOfCollection;
	for(int i=0;i<sz;++i){
		id = map_addcovers[i];
 		if(id>=0 ){
			dual_sol[cov->id_vi] = alph1[id];
		}else dual_sol[cov->id_vi]=0;
		cov = cov->next;
	}
	alph1.end();
	
	IloNumArray alph2(env);
	cplex.getDuals(alph2, add_locals);
    LocalCut * loc = localcs->begin;
	sz = localcs->sizeOfCollection;
	for(int i=0;i<sz;++i){
		id = map_addlocals[i];
 		if(id>=0 ){
			dual_sol[loc->id_vi] = alph2[id];
		}else dual_sol[loc->id_vi]=0;
		loc = loc->next;
	}
	alph2.end();

	IloNumArray alph3(env);
	cplex.getDuals(alph3, add_globals);
    GlobalCut * gloc = globalcs->begin;
	sz = globalcs->sizeOfCollection;
	for(int i=0;i<sz;++i){
		id = map_addglobals[i];
 		if(id>=0 ){
			dual_sol[gloc->id_vi] = alph3[id];
		}else dual_sol[gloc->id_vi]=0;
		gloc = gloc->next;
	}
	alph3.end();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
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
LPChecker::logical_fix(double ub){
	int arc;
	double lb = cplex.getObjValue();
	IloNumArray rc(env);
	cplex.getReducedCosts(rc, y);
	for(int a=0;a<szunfx;++a){
		arc = unfixd[a];
		
	}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void
LPChecker::clean(){
	add_strong_force.endElements();
	add_covers.endElements();
	add_locals.endElements();
	add_globals.endElements();
}
	

//---------------------------------------------------------------------------


LPChecker::~LPChecker() {

	try {
		cplex.clearModel();

		x.endElements();
		y.endElements();
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
		if(dual_sol) delete [] dual_sol;

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}
