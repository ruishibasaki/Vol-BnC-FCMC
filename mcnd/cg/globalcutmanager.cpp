//
//  globalcutmanager.cpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#include "globalcutmanager.hpp"

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  initializing methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void
GlobalCutManager::initialize(const Data * d, int lim) {
    data=d;
    globals.initialize(data->narcs);
 
    lim_to_remv = lim;
    num_actv  = 0;
}

//-------------------------------------------------------------------------------------------

int
GlobalCutManager::reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize, bool recheck_collct){
    int cont;
    int sz = globals.sizeOfCollection;
    GlobalCut* vi = globals.begin;
    num_actv =0; cont=0;
    bool put=true;
    bool infeas;
    for(int i=0;i<sz;++i){
        vi->n_zerom =0;
        vi->n_nviol = 0;
        put=true;
        infeas=false;

        if(recheck_collct) put = vi->check_updt_Viol(topo, infeas);
    	if(infeas) return -1;
		put=false;
        if(put && !vi->purgbl){
            actvS[vi->id_vi] = fsize+csize;
            //std::cout<<"in: "<<vi->serial_nmbr<<" id: "<<vi->id_vi<<" "<<recheck_collct<<std::endl;
            //vi->print();
            ++csize;
            ++num_actv;
            vi = vi->next;
        }else{
        	//std::cout<<"out: "<<vi->serial_nmbr<<" id: "<<vi->id_vi<<std::endl;
        	//vi->print();
        	++cont;
        	vi = globals.move_to_end(vi);
        }
    }
    if(sz)globals.begin->prev = globals.end->next = 0;

    return cont;
}

//-------------------------------------------------------------------------------------------

void
GlobalCutManager::clean_collection(){
	globals.begin = globals.end = 0;
	globals.sizeOfCollection = 0;
	globals.empty = true;
 	 num_actv = 0;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  main methods
//-------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

int
GlobalCutManager::globalc_generation_main2( const double * colub, const int * topo){
    int clsd=0;
	
	int * closed_arcs = new int[data->narcs];
	if(colub !=0 ){
		for(int a=data->narcs;a--;){
			if(colub[a]<0.5) 
				closed_arcs[clsd++]=a;
		}
	}else if(topo !=0 ){
		for(int a=data->narcs;a--;){
			if(topo[a]==0) 
				closed_arcs[clsd++]=a;
		}
	}else{
	 	delete [] closed_arcs;
		return 0;
	} 
	
	int * vars_ =  new int[clsd];
	for(int a=clsd;a--;){
 		 vars_[a] = closed_arcs[a];
	}
 	delete [] closed_arcs;
    int added=0;
	added += make_globalcut( 0, clsd, vars_, 0);	
 	//if(added)std::cout<<"GlobalCutManager::globalc_generation_main2 added global serial_num: "<<ttgend-1<<std::endl;
 	return added;
}
 
//-------------------------------------------------------------------------------------------

int
GlobalCutManager::globalc_generation_main( const double * ystar, const int * closed, int cont0, int curr_id){
    int * vars_;
	 
	vars_ = new int[cont0];
	for(int a=cont0;a--;){
		 vars_[a]=closed[a];
	}
	
    int added=0;
	added += make_globalcut( ystar, cont0, vars_, curr_id);	
 
 	return added;
}
 
//-------------------------------------------------------------------------------------------

int 
GlobalCutManager::make_globalcut(const double * ystar, int sz,  int* vars_, int curr_id){
  	double suml=0;
	int id_arc;
	GlobalCut *gloc=0;
	if(ystar){
		for(int i=sz;i--;){
  			suml += ystar[vars_[i]];
 		}
 	}
	//std::cout<<std::endl;
 	gloc = globals.createNewGlobalCut(sz, vars_, curr_id, ttgend);
 
    if(gloc!=0){
    	//std::cout<<"GlobalCutManager::make_globalcut try add: "<<std::endl;
        int added = globals.addGlobalCut(gloc);
         if(added){
        	gloc->hs = (1.0 - suml);
        	//std::cout<<"GlobalCutManager::make_globalcut add global: "<<std::endl;
        	//gloc->print();
        	++ttgend;
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
GlobalCutManager::reposition_globals(int added){	
	if(globals.sizeOfCollection == num_actv) return;
	
	int num_adv = num_actv - added;
	int sz = globals.sizeOfCollection;
	GlobalCut * last_actv = globals.begin ;
	globals.advance(last_actv, num_adv-1);
	GlobalCut * trgt = globals.end ;
	
	if(num_adv == 0){ 
		for(int i = added; i-- ; ){
			globals.begin = trgt;
			globals.end = trgt->prev;
			globals.end->next = 0;
			trgt->prev = 0;
			trgt->next = last_actv;
			last_actv->prev = trgt;
			last_actv = trgt;
			trgt = globals.end ;
		}
		return;
	}

	for(int i = added; i--  ; ){
		globals.end = trgt->prev;
		globals.end->next = 0;
		trgt->next = last_actv->next;
		last_actv->next->prev = trgt;
		last_actv->next = trgt;
		trgt->prev = last_actv;
		trgt = globals.end ;
	}
	
}

//-------------------------------------------------------------------------------------------

void 
GlobalCutManager::collect_globals(const BCP_vec<BCP_var*>& vbd, int & currnum){
	int sz = globals.track.size();
  	std::deque<GlobalCut *>& trackgloc = globals.track;
  	GlobalCut *gloc;
  	for(int i=sz; i--;){
  		gloc = trackgloc[i];
  		gloc->purgbl=false;
  		if(gloc->check_viol(vbd)){
  			//std::cout<<"GlobalCutManager::collect_globals:: gloc "<<gloc->serial_nmbr<<std::endl;
  			gloc->id_vi = currnum++;
  			globals.insert_end(gloc);
  		}
	}
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  Volume Integration methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


void
GlobalCutManager::add_global_vi(int added, int * actvS, int & actvSSz,  double * dualsol, double *lhsol,
    				 			double * h, double * dstar, double * dual_lb, double * dual_ub){
    //std::cout<<"add_local_vi: "<<actvSSz<<std::endl;
    int idx;
    GlobalCut * vi = globals.end;
    for(int cont = added; cont--;){
        idx = actvSSz+cont; //if(actvS[vi->id_vi]>=0){ std::cout<<actvS[vi->id_vi]<<" already taken !!!!!!! for: "<<vi->id_vi<<std::endl;abort();}
        //std::cout<<"add vi: "<< vi->id_vi<<" idx: "<<idx<<std::endl;
        actvS[vi->id_vi]=idx;
        dstar[idx] =0;
        dual_lb[idx] = 0;
        dual_ub[idx] = 1e31;
        h[idx] = vi->hs;
        
        dualsol[vi->id_vi] =0;
        lhsol[vi->id_vi] = 0;

        //vi->print();
        vi = vi->prev;
    }
    actvSSz += added;
}

//-------------------------------------------------------------------------------

int
GlobalCutManager::compute_cover_sg( const double * x, const int * actvS, int actvSSz,  double * v){
    //std::cout<<"compute_flowpack_sg"<<std::endl;
    int index, id_arc;
    int sz = num_actv;
    GlobalCut *vi = globals.begin;
    for(int n=0;n<sz;++n){
        index = actvS[vi->id_vi];
        
        v[index] = 1.0;
        for(int a=vi->size;a--;){
            id_arc = arc_map[vi->vars[a]];
            if(id_arc<0) continue;
            v[index] -=  x[id_arc];
        } 
		if(index>=actvSSz){ std::cout<<"globalindex: "<<index<<"/"<<actvSSz<<std::endl; abort(); }
        if(v[index]<=0){
            ++vi->n_nviol;
            if(vi->n_nviol>=lim_to_remv && vi->n_zerom>0) v[index]=0;
        }else vi->n_nviol=0;
    	
    	//if(v[index]>0){ std::cout<<"global: "<<vi->serial_nmbr<<" = "<<v[index]<<std::endl;}
         vi = vi->next;
    }    
    return 0;
}

//-------------------------------------------------------------------------------

int
GlobalCutManager::compute_cover_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0){

    int index, id_arc;
    int sz = num_actv;
    GlobalCut *vi = globals.begin;
    for(;sz--;){   
        index = actvS[vi->id_vi];
        if(dual[index]==0){
            ++vi->n_zerom;
            vi = vi->next;
            continue;        
        }else vi->n_zerom = 0; 

        B0 += dual[index];
        for(int a=vi->size;a--;){
            id_arc = arc_map[vi->vars[a]];
            if(id_arc<0) continue;
            rc[id_arc]-= dual[index];   
        }   
        vi = vi->next;
    }
    return 0;
}

//-------------------------------------------------------------------------------

double
GlobalCutManager::arc_dg_imp(int arc, const double * xy, const double * h, const int * actvS, int actvSSz){
    int index;
    int sz = num_actv;
    GlobalCut *vi = globals.begin;
    double gam;
    double dg=0;
    for(;sz--;){
        index = actvS[vi->id_vi];
        gam = globals.GlobalCut_hasArc(vi, arc);
        dg += gam*h[index];
        vi = vi->next;
    }
    return dg;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  manager methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void 
GlobalCutManager::collect_global(double* collb, double* colub, int & curr_numrows){
	const GlobalCut* vi=0;
	for(int i=globals.track.size();i--;){
		vi = globals.track[i];
		
	}

}

















