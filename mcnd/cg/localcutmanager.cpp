//
//  localcutmanager.cpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#include "localcutmanager.hpp"

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  initializing methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void
LocalCutManager::initialize(const Data * d, int lim) {
    data=d;
    locals.initialize(data->narcs);
    fixbl_arcs.resize(data->narcs,-1);

    lim_to_remv = lim;
    num_actv = gend = 0;
}

//-------------------------------------------------------------------------------------------

int
LocalCutManager::reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize, bool recheck_collct){
    int cont;
    int sz = locals.sizeOfCollection;
    LocalCut* vi = locals.begin;
    num_actv =0; cont=0;
    fixbl_arcs.assign(data->narcs,-1);
    bool put=true;
    bool infeas;
    for(int i=0;i<sz;++i){
        vi->n_zerom =0;
        vi->n_nviol = 0;
        put=true;
        infeas=false;
        //if(vi->type==0)std::cout<<"try: type "<<vi->type<<" "<<vi->prgbl<<std::endl;
        if(recheck_collct) put = vi->check_updt_Viol(topo, infeas);
        if(infeas) return -1;
		put=false;
        
        if(put && !vi->prgbl){
            actvS[vi->id_vi] = fsize+csize;
            //if(vi->type==0)std::cout<<"in type: "<<vi->type<<" srial: "<<vi->serial_nmbr<<std::endl;
            ++csize;
            ++num_actv;
            vi = vi->next;
        }else{
        	//std::cout<<"out: "<<vi->serial_nmbr<<" id: "<<vi->id_vi<<std::endl;
        	//vi->print();
        	++cont;
        	vi = locals.move_to_end(vi);
        }
    }
    if(sz)locals.begin->prev = locals.end->next = 0;

    return cont;
}

//-------------------------------------------------------------------------------------------

void
LocalCutManager::clean_collection(){
	LocalCut * loc = locals.begin;
	for(int i=locals.sizeOfCollection;i--;){
		locals.begin = locals.begin->next;
		if(loc->toadd) delete loc;
		loc = locals.begin;
	}
	locals.begin = locals.end = 0;
	locals.sizeOfCollection = locals.discarted = 0;
	locals.empty = true;
	fixbl_arcs.assign(data->narcs,-1);
	 num_actv = gend = 0;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  main methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int 
LocalCutManager::localc0_generation_main( const double * ystar,  const int * closed, int sz, int curr_id){
	int opnd=0;
    int * vars_;
	 
	vars_ = new int[sz];
	for(int a=sz;a--;){
		 vars_[a]=closed[a];
		 //std::cout<<"closed: "<<closed[a]<<std::endl;
	}
	
    int added=0;
	added += make_localcut0( ystar, sz, vars_, curr_id);	
 
 	return added;
}

//-------------------------------------------------------------------------------------------

int
LocalCutManager::make_localcut0(const double * ystar, int sz,  int* vars_, int curr_id){
  	double suml=0;
	int id_arc;
	LocalCut *loc=0;
	if(ystar){
		for(int i=sz;i--;){
  			suml += ystar[vars_[i]];
 		}
 	}
	//std::cout<<std::endl;
 	loc = locals.createNewLocalCut(sz, vars_, curr_id, ttgend);
 
    if(loc!=0){
    	//std::cout<<"GlobalCutManager::make_globalcut try add: "<<std::endl;
        int added = locals.addLocalCut(loc);
         if(added){
        	loc->hs = (1.0 - suml);
        	//std::cout<<"GlobalCutManager::make_globalcut add global: "<<std::endl;
        	//gloc->print();
        	++ttgend;
        	++gend;
        	++num_actv;
        	//++ttgend_feas; //stat can be removed
        	return 1;
        }
    }
    return 0;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
LocalCutManager::localc1_generation_main(double lb, double ub, const double * ystar, const int * y, const double * rc, int curr_id, int max){
    int added=0;
    int id_arc, arc;
    double rc_val;
    double rc_ttp=0;
    double rc_ttm=0;
	std::list<Pair2> rc_p;
	std::list<Pair2> rc_m;
	std::vector<int>Tp;
	std::vector<int>Tm;
	
    for(int a=data->narcs;a--;){
     	if(y[a]>=0) continue;
    	rc_val=rc[a];
    	if(rc_val>0)
    		rc_p.push_back(Pair2(a, rc_val));
    	else if(rc_val<0) rc_m.push_back(Pair2(a, -rc_val));
    }
    
    rc_p.sort(compPair2());//decreasing
    rc_m.sort(compPair2());//decreasing
	//std::cout<<"positive"<<std::endl;
	rc_ttp = check_fixable(rc_p, Tp, lb, ub ,0);
	form_t(rc_p, Tp, lb, ub, rc_ttp);
	added += make_localcut1( Tp, ystar, curr_id, 0);
	if(added+curr_id>=max)return added;
	
	//std::cout<<"negative"<<std::endl;
	rc_ttm = check_fixable(rc_m, Tm, lb, ub ,1);
	form_t(rc_m, Tm, lb, ub, rc_ttm);
	added += make_localcut1( Tm, ystar, curr_id+added, 1);

 	return added;
}
 
//-------------------------------------------------------------------------------------------

int 
LocalCutManager::make_localcut1(std::vector<int>& T, const double * ystar, int curr_id, int oneor0){
	double rhs;
 	double suml=0;
	int id_arc;
	LocalCut *loc=0;
	if(T.size()<=1){T.clear(); return 0;}
	//std::cout<<" make localcut: ";
	for(int i=T.size();i--;){
		id_arc = arc_map[T[i]];
 		suml += ystar[id_arc];
		//std::cout<<T[i]<<" ";
	}
	//std::cout<<std::endl;
	if(oneor0==0){
		rhs=1;
		loc = locals.createNewLocalCut(T, curr_id, ttgend, -1, rhs);
	}else{
		rhs = (T.size()-1);
		loc = locals.createNewLocalCut(T, curr_id, ttgend, 1, rhs);
	}
	
	T.clear();
    if(loc!=0){
        int added = locals.addLocalCut(loc);
        if(added){
        	loc->hs = loc->sense*(rhs - suml);
        	++ttgend;
        	++gend;
        	++num_actv;
        	//++ttgend_sellm; //stat can be removed
        	return 1;
        }
    }
    return 0;
}

//-------------------------------------------------------------------------------------------

void
LocalCutManager::form_t(std::list<Pair2>& rc_, std::vector<int>& T, double lb, double ub, double rc_tt){
	int arc;
	double rc_val;

	while(!rc_.empty()){
		arc = rc_.front().fst;
		rc_val = rc_.front().snd;
		rc_.pop_front();
 		rc_tt+=rc_val;
 		//std::cout<<"arc "<<arc<<" rc: "<<rc_val<<" sum: "<<lb+rc_tt<<" ub: "<<ub<<std::endl;
		if(lb+rc_tt>=ub){
			T.push_back(arc);
			rc_tt=rc_val;
		}else break;
	}
	rc_.clear();
}

//-------------------------------------------------------------------------------------------

int
LocalCutManager::check_fixable(std::list<Pair2>& rc_, std::vector<int>& T, double lb, double ub, int oneor0){
 	int arc, fixto;
	double rc_val;
		
	if(oneor0==0){fixto=0;}
	else{ fixto=1;}
	
	while(!rc_.empty()){
		arc = rc_.front().fst;
		rc_val = rc_.front().snd;
		rc_.pop_front();
 		if(lb+rc_val>=ub){
 			std::cout<<"fixable "<<arc<<" to "<<fixto<<std::endl; abort();
			fixbl_arcs[arc]= fixto;
		}else{
		 	T.push_back(arc); 
		 	return rc_val;
		 }
 	}
	return 0;
}
 
 
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------

int 
LocalCutManager::localc2_generation_main( const double * ystar,  const double * topo, int sz, int curr_id){
    int * vars_;
    double* coef_;
    double rhs=1.0;
  		
	vars_ = new int[data->narcs];
	coef_ = new double[data->narcs];
	for(int a=data->narcs;a--;){
		vars_[a]=a;
		if(topo[a]<0.5){
			coef_[a] = 1.0;
		}else if(topo[a]>0.5){
			coef_[a] = -1.0;
			rhs-=1.0;
		}else{ std::cout<<"LocalCutManager::localc2_generation_main"<<std::endl;  abort();}
	}
	
	int added=0;
	added += make_localcut2( ystar, data->narcs, vars_, coef_, rhs, curr_id);	
 
 	return added;
}

//-----------------------------------------------------------------------------------

int 
LocalCutManager::make_localcut2(const double * ystar, int sz,  int* vars_, double* coef_, double rhs_, int curr_id){
	double suml=0;
	int id_arc;
	LocalCut *loc=0;
	if(ystar){
		for(int i=sz;i--;){
  			suml += coef_[i]*ystar[vars_[i]];
 		}
 	}
	//std::cout<<std::endl;
 	loc = locals.createNewLocalCut(sz, vars_, coef_, rhs_, curr_id, ttgend);
 
    if(loc!=0){
    	//std::cout<<"GlobalCutManager::make_globalcut try add: "<<std::endl;
         int added = locals.addLocalCut(loc);
         if(added){
        	loc->hs = (rhs_ - suml);
        	//std::cout<<"GlobalCutManager::make_globalcut add global: "<<std::endl;
        	//gloc->print();
        	++ttgend;
        	++gend;
        	++num_actv;
        	//++ttgend_opt; //stat can be removed
        	return 1;
        }
    }
    return 0;
}
	
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------


void
LocalCutManager::reposition_locals(int added){	
	if(locals.sizeOfCollection == num_actv) return;
	
	int num_adv = num_actv - added;
	int sz = locals.sizeOfCollection;
	LocalCut * last_actv = locals.begin ;
	locals.advance(last_actv, num_adv-1);
	LocalCut * trgt = locals.end ;
	
	if(num_adv == 0){ 
		for(int i = added; i-- ; ){
			locals.begin = trgt;
			locals.end = trgt->prev;
			locals.end->next = 0;
			trgt->prev = 0;
			trgt->next = last_actv;
			last_actv->prev = trgt;
			last_actv = trgt;
			trgt = locals.end ;
		}
		return;
	}

	for(int i = added; i--  ; ){
		locals.end = trgt->prev;
		locals.end->next = 0;
		trgt->next = last_actv->next;
		last_actv->next->prev = trgt;
		last_actv->next = trgt;
		trgt->prev = last_actv;
		trgt = locals.end ;
	}
	
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  Volume Integration methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


void
LocalCutManager::add_local_vi(int added, int * actvS, int & actvSSz, double * dualsol, double *lhsol,
    				 double * h, double * dstar, double * dual_lb, double * dual_ub ){
    //std::cout<<"add_local_vi: "<<actvSSz<<std::endl;
    int idx;
    LocalCut * vi = locals.end;
    for(int cont = added; cont--;){
        idx = actvSSz+cont; //if(actvS[vi->id_vi]>=0){ std::cout<<actvS[vi->id_vi]<<" already taken !!!!!!! for: "<<vi->id_vi<<std::endl;abort();}
        //std::cout<<"add vi: "<< vi->id_vi<<" idx: "<<idx<<std::endl;
        actvS[vi->id_vi]=idx;
        dstar[idx] =0;
        dual_lb[idx] = 0;
        dual_ub[idx] = 1e31;
        h[idx] = vi->hs;
        
        lhsol[vi->id_vi]=0;
        dualsol[vi->id_vi] =0;

        
        //vi->print();
        vi = vi->prev;
    }
    actvSSz += added;
}

//-------------------------------------------------------------------------------

int
LocalCutManager::compute_localc_sg( const double * x, const int * actvS, int actvSSz,  double * v){
    //std::cout<<"compute_flowpack_sg"<<std::endl;
    int index, id_arc;
    int sz = num_actv;
    LocalCut *vi = locals.begin;
    for(int n=0;n<sz;++n){
        index = actvS[vi->id_vi];
        v[index] = vi->get_total_rhs();
        for(int a=vi->size;a--;){
            id_arc = arc_map[vi->vars[a]];
            if(id_arc<0) continue;
            v[index] -=  vi->coef_at(a)*x[id_arc];
        }
		
		if(index>=actvSSz){ std::cout<<"localindex: "<<index<<"/"<<actvSSz<<std::endl; abort(); }
		v[index] *= vi->sense;
		
		//if(vi->type!=1){vi->print(); std::cout<<"viol "<<v[index]<<" rhs: "<<vi->get_total_rhs()<<" type "<<vi->type<<std::endl;}

        if(v[index]<=0){
            ++vi->n_nviol;
            if(vi->n_nviol>=lim_to_remv && vi->n_zerom>0) v[index]=0;
        }else vi->n_nviol=0;   
        vi = vi->next;
    }    
    return 0;
}

//-------------------------------------------------------------------------------

int
LocalCutManager::compute_localc_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0){

    int index, id_arc;
    int sz = num_actv;
    LocalCut *vi = locals.begin;
    for(;sz--;){   
        index = actvS[vi->id_vi];
        //std::cout<<" idx: "<<index<<std::endl;
        if(dual[index]==0){
            ++vi->n_zerom;
            vi = vi->next;
            continue;        
        }else vi->n_zerom = 0; 

        B0 +=   vi->sense*dual[index]*vi->get_total_rhs();
        for(int a=vi->size;a--;){
            id_arc = arc_map[vi->vars[a]];
            if(id_arc<0) continue;
            rc[id_arc]-= vi->sense*vi->coef_at(a)*dual[index];   
        }   
        vi = vi->next;
    }
    return 0;
}

//-------------------------------------------------------------------------------

double
LocalCutManager::arc_dg_imp(int arc, const double * xy, const double * h, const int * actvS, int actvSSz){
    int index;
    int sz = num_actv;
    LocalCut *vi = locals.begin;
    double gam;
    double dg=0;
    for(;sz--;){
        index = actvS[vi->id_vi];
        gam = locals.LocalCut_hasArc(vi, arc);
        dg += -vi->sense*gam*h[index];
        vi = vi->next;
    }
    return dg;
}







