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
    num_actv = gend = ttgend= 0;
}

//-------------------------------------------------------------------------------------------

int
LocalCutManager::reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize){
    int cont;
    int idxf = data->nnodes*data->ndemands;
    int sz = locals.sizeOfCollection;
    LocalCut* vi = locals.begin;
    num_actv =0; cont=0;
    fixbl_arcs.assign(data->narcs,-1);
    //std::cout<<"sz: "<<sz<<std::endl;
    for(int i=0;i<sz;++i){
        vi->n_zerom =0;
        vi->n_nviol = 0;
        //std::cout<<"in: id_vi "<<vi->id_vi<<std::endl;
        if(vi->check_updt_Viol(topo) && !vi->prgbl){
            actvS[vi->id_vi] = fsize+csize;
            //std::cout<<"in: "<<vi->serial_nmbr<<" id: "<<vi->id_vi<<std::endl;
            //vi->print();
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
    return cont;
}

//-------------------------------------------------------------------------------------------

void
LocalCutManager::clean_collection(){
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
LocalCutManager::localc_generation_main(double lb, double ub, const double * ystar, const double * y, const double * rc, int curr_id, int max){
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
    	id_arc = arc_map[a];
    	if(id_arc<0) continue;
    	rc_val=rc[id_arc];
    	if(rc_val>0)
    		rc_p.push_back(Pair2(a, rc_val));
    	else if(rc_val<0) rc_m.push_back(Pair2(a, -rc_val));
    }
    
    rc_p.sort(compPair2());//decreasing
    rc_m.sort(compPair2());//decreasing
	std::cout<<"positive"<<std::endl;
	rc_ttp = check_fixable(rc_p, Tp, lb, ub ,0);
	form_t(rc_p, Tp, lb, ub, rc_ttp);
	added += make_localcut( Tp, ystar, y, curr_id+added, 0);
	if(added+curr_id>=max)return added;
	
	std::cout<<"negative"<<std::endl;
	rc_ttm = check_fixable(rc_m, Tm, lb, ub ,1);
	form_t(rc_m, Tm, lb, ub, rc_ttm);
	added += make_localcut( Tm, ystar, y, curr_id, 1);

 	return added;
}
 
//-------------------------------------------------------------------------------------------

int 
LocalCutManager::make_localcut(std::vector<int>& T, const double * ystar, const double * y, int curr_id, int oneor0){
	double rhs;
	double sum=0;
	double suml=0;
	int id_arc;
	LocalCut *loc=0;
	std::cout<<" make localcut: ";
	for(int i=T.size();i--;){
		id_arc = arc_map[T[i]];
		sum += y[id_arc];
		suml += ystar[id_arc];
		std::cout<<T[i]<<" ";
	}
	std::cout<<std::endl;
	if(oneor0==0){
		std::cout<<"sum: "<<sum<<" suml: "<<suml<<" rhs: "<<1<<std::endl;
		if(sum > 1 && suml>1){
			loc = locals.createNewLocalCut(T, curr_id, ttgend, 2, 1.0);
			
		}
	}else{
		std::cout<<"sum: "<<sum<<" suml: "<<suml<<" rhs: "<<(T.size()-1)<<std::endl;
		rhs = (T.size()-1);
		if(sum < rhs && suml<rhs){
			loc = locals.createNewLocalCut(T, curr_id, ttgend, 1, rhs);
		}
	}
	
    if(loc!=0){
        int added = locals.addLocalCut(loc);
        if(added){
        	loc->hs = suml;
        	++ttgend;
        	++gend;
        	++num_actv;
        	return 1;
        }else delete loc;
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
 		std::cout<<"arc "<<arc<<" rc: "<<rc_val<<" sum: "<<lb+rc_tt<<" ub: "<<ub<<std::endl;
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
 			//std::cout<<"fixable "<<arc<<" to "<<fixto<<std::endl;
			fixbl_arcs[arc]= fixto;
		}else{
		 	T.push_back(arc); 
		 	return rc_val;
		 }
 	}
	return 0;
}
 
