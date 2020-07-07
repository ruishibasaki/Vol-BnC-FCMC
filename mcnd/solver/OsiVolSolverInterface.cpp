// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdlib>
#include <numeric>
#include <cassert>
#include <cmath>
#include <iomanip>

#include "CoinHelperFunctions.hpp"

#include "OsiVolSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// warm start methods
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

CoinWarmStart *
OsiVolSolverInterface::getEmptyWarmStart () const{
	//std::cout<<"getEmptyWarmStart"<<std::endl;
    return new WarmStartDual();
}

//---------------------------------------------------------------------------

CoinWarmStart* 
OsiVolSolverInterface::getWarmStart() const{
    //std::cout<<"getWarmStart() "<<getNumRows()<<std::endl;
    return new WarmStartDual(narcs, solution, getNumRows(), dual, &cover_manager->covers, &localc_manager->locals, &globalc_manager->globals);
    
} 

//---------------------------------------------------------------------------

bool 
OsiVolSolverInterface::setWarmStart(const CoinWarmStart* warmstart){
    const WarmStartDual* ws = dynamic_cast<const WarmStartDual*>(warmstart);
	
    if (! ws){
        HotStartSet = false;
    	return false;
    }
    if(HotStart_) delete HotStart_;
    HotStart_ = ws->clone_ws();
    //std::cout<<"setWarmStart "<<HotStart_->size()<<" "<<getNumRows()<<std::endl;
    HotStartSet = true;
    return true;
};

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::markHotStart(){
    //std::cout<<"markHotStart() "<<getNumRows()<<std::endl;
    if(HotStart_) delete HotStart_;
    HotStart_ = new WarmStartDual(narcs, solution, getNumRows(), dual, &cover_manager->covers, &localc_manager->locals, &globalc_manager->globals); 
    //std::cout<<"markHotStart() now "<<std::endl;
    HotStartSet = true;
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::unmarkHotStart(){
    //std::cout<<"unmarkHotStart() "<<std::endl;
    if(HotStart_) delete HotStart_;
    HotStart_ = 0;
    HotStartSet = false;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// 				SETTERS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

void 
OsiVolSolverInterface::setRowSetBounds(const int* indexFirst,
                                       const int* indexLast,
                                       const double* boundList){
    //std::cout<<"rowsetbounds "<<std::endl;
    while (indexFirst != indexLast) {
        rowlb[*indexFirst] =  boundList[0];
        rowub[*indexFirst] = boundList[1];
        ++indexFirst;
        boundList += 2;
    }
    //std::cout<<"ok"<<std::endl;
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::setColSetBounds(const int* indexFirst,
                                       const int* indexLast,
                                       const double* boundList){
    //std::cout<<"colsetbounds "<<std::endl;
    while (indexFirst != indexLast) {
        collb[*indexFirst] =  boundList[0];
        colub[*indexFirst] = boundList[1];
         //std::cout<<"col: "<<*indexFirst<<" = "<<boundList[0]<<" "<<boundList[1]<<std::endl;
        ++indexFirst;
        boundList += 2;
    }
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::map_topology(){
	nz_arcs.clear();
    szopnd=0;
    szunfxd=0;
    for(int arc= narcs;arc--;){
        if(collb[arc] == 1){
            nz_arcs.push_back(arc);
            arc_map[arc] = -2;
            VItopo[arc] =  yhit[arc] = 1.0;
            ++szopnd;
             //std::cout<<"opened arc: "<<arc<<std::endl;
        }else if(colub[arc] == 1){
            ++szunfxd;
            nz_arcs.push_front(arc);
            VItopo[arc] = 0.0;
            yhit[arc] = -1.0;
            //std::cout<<"unfix: "<<arc<<" idx: "<<arc_map[arc]<<std::endl;
        }else{
        	VItopo[arc] = yhit[arc] = 0.0;
             //std::cout<<"closed arc: "<<arc<<std::endl;
            arc_map[arc] = -1;
        }
    }

    for(int i = 0; i<szunfxd;++i){
    	arc_map[nz_arcs[i]] = i;
    }
    
    sznz = nz_arcs.size();
    cover_manager->set_arc_map(arc_map);
	localc_manager->set_arc_map(arc_map);
	globalc_manager->set_arc_map(arc_map);
}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::map_duals(){
    fsize= csize = 0;
    bool flag=false;
    const Arc* item; const Demand* itemd;
    CoinFillN(actv, numrows_, -1);
    
    for(int k=0; k<ndemands; ++k ){
        for(int i=0; i<nnodes; ++i ){
            flag = false;
            
            for(int a=sznz; a--;){
                int arc = nz_arcs[a];
                item = &data->arcs[arc];
                if((i+1) == item->i && (colub[narcs+k*narcs+arc]>0)){
                    flag = true;
                    break;
                }else if((i+1) == item->j && (colub[narcs+k*narcs+arc]>0)){
                    flag = true;
                    break;
                }
            }
            itemd = &data->d_k[k];
            if( (i+1) == itemd->D){
                flag = true;
            }else if( (i+1) == itemd->O){
                flag = true;
            }
            if(flag){
                actv[k*nnodes + i] = fsize++;
            }/*else{
            	HotStartSet =false;
             }*/
        }
    }
	//std::cout<<"dual size: "<<fsize<<" < "<<nnodes*ndemands<<std::endl;
    int ret = cover_manager->reset_and_map_collection(fsize, yhit, dual, actv, csize, recheck_collct);
    int ret2 = localc_manager->reset_and_map_collection(fsize, yhit, dual, actv, csize, recheck_collct);
    int ret3 = globalc_manager->reset_and_map_collection(fsize, yhit, dual, actv, csize, recheck_collct);
	if(ret<0 || ret2<0 || ret3<0 ){ mode=-2; }
	//std::cout<<"OsiVolSolverInterface::map_duals()::INFEASIBLE DETECTED "<<ret<<" "<<ret2<<" "<<ret3<<std::endl;}
	//std::cout<<"fsize: "<<fsize<<" "<<csize<<std::endl;
    /*if(ret>num_purgbl && in_strong_branch){
    	HotStartSet =true;	
    }else*/
    /*if(HotStart_ !=0 ){
    	HotStartSet =true;	
    }else if(!in_strong_branch){
    	num_purgbl=ret;
    }*/
    
    //std::cout<<"OsiVolSolverInterface::map_duals "<<in_strong_branch<<" "<<ret<<" "<<HotStartSet<<std::endl;
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::translate_hotstart(){
    int fidx = ndemands*nnodes;
    int idx, sz;
 	for(int i=fidx; i--;){
		idx = actv[i];
		if(idx>=0){
			volprob_->dsol[idx] = HotStart_->dual[i];
		}
	}
	sz = cover_manager->covers.sizeOfCollection;
	Cover* vi = cover_manager->covers.begin;
	for(int i=sz; i--;){
		idx = actv[vi->id_vi];
		if(idx>=0){
			volprob_->dsol[idx] = HotStart_->get_mapped(vi->serial_nmbr);
			//std::cout<<" wsvi idx "<<idx<<" serial "<<vi->serial_nmbr<<" : "<<volprob_->dsol[idx]<<std::endl;
		}//else std::cout<<idx<<"wsvi id "<<vi->id_vi<<" : "<<HotStart_->dual_[idx]<<std::endl;
		vi = vi->next;
	}

	sz = localc_manager->locals.sizeOfCollection;
	LocalCut* vilc = localc_manager->locals.begin;
	for(int i=sz; i--;){
		idx = actv[vilc->id_vi];
		if(idx>=0){
			volprob_->dsol[idx] = HotStart_->get_mapped(vilc->serial_nmbr);
			//std::cout<<" wsvi idx "<<idx<<" serial "<<vilc->serial_nmbr<<" : "<<volprob_->dsol[idx]<<std::endl;
		}//else std::cout<<idx<<"wsvi id "<<vi->id_vi<<" : "<<HotStart_->dual_[idx]<<std::endl;
		vilc = vilc->next;
	}

	sz = globalc_manager->globals.sizeOfCollection;
	GlobalCut* vigc = globalc_manager->globals.begin;
	for(int i=sz; i--;){
		idx = actv[vigc->id_vi];
		if(idx>=0){
			volprob_->dsol[idx] = HotStart_->get_mapped(vigc->serial_nmbr);
			//std::cout<<" wsvi idx "<<idx<<" serial "<<vigc->serial_nmbr<<" : "<<volprob_->dsol[idx]<<std::endl;
		}//else std::cout<<idx<<"wsvi id "<<vi->id_vi<<" : "<<HotStart_->dual_[idx]<<std::endl;
		vigc = vigc->next;
	}

}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::set_start(){
    //std::cout<<"OsiVolSolverInterface::set_start "<<std::endl;
    int idx, sz;
    int fidx = ndemands*nnodes;
   
    VItt = VIub=-1e31;
    cover_manager->colub = colub;
    CoinFillN(volprob_->dual_ub.v, getNumRows(), 0.0);
    CoinFillN(volprob_->dual_lb.v, getNumRows(), 0.0);
    CoinFillN(volprob_->dsol.v, getNumRows(), 0.0);
    CoinFillN(rc_, narcs, 0.0);

    for(int i=fidx; i--;){
        //std::cout<<old_index[i]<<" "<<Iu[i]<<" value: ";
        //std::cout<<old_dual[old_index[i]]<<std::endl;
        idx = actv[i];
        if(idx>=0){
            volprob_->dual_lb[idx] = -1.0e31;
            volprob_->dual_ub[idx] = 1.0e31;
        } 
    }
    sz = cover_manager->num_actv;
    Cover* vi = cover_manager->covers.begin;
    for(int i=sz; i--;){
        idx = actv[vi->id_vi];
        if(idx>=0){
            volprob_->dual_ub[idx] = 1.0e31;
        }
        vi = vi->next;
    }
    sz = localc_manager->num_actv;
    LocalCut* vilc = localc_manager->locals.begin;
    for(int i=sz; i--;){
        idx = actv[vilc->id_vi];
        if(idx>=0){
            volprob_->dual_ub[idx] = 1.0e31;
        }
        vilc = vilc->next;
    }
    sz = globalc_manager->num_actv;
    GlobalCut* vigc = globalc_manager->globals.begin;
    for(int i=sz; i--;){
        idx = actv[vigc->id_vi];
        if(idx>=0){
            volprob_->dual_ub[idx] = 1.0e31;
        }
        vigc = vigc->next;
    }
    
    
    if(HotStartSet){
    	if(volprob_->parm.maxsgriters>250)
    		volprob_->parm.maxsgriters=250;
		translate_hotstart();
    }else volprob_->parm.maxsgriters=1000;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Solve triggers
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


void 
OsiVolSolverInterface::solveFromHotStart(){
    //std::cout<<"solveFromHotStart() maxiter: "<<volprob_->parm.maxsgriters<<" dsize: "<<numrows_<<std::endl;
    //CoinDisjointCopyN(HotStart_, numrows_, dual);
    //OsiVolAuxInfo * auxinfo  = static_cast<OsiVolAuxInfo*>(OsiSolverInterface::getApplicationData());
    resolve();
    
}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::initialSolve(){
    // set every entry to 0.0 in the dual solution
    //std::cout<<"INITIAL SOLVE"<<std::endl;
    //CoinFillN(u, getNumRows(), 0.0);
    //resolve();
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::resolve(){
    //std::cout<<"mode: "<<mode<<" numrows: "<<numrows_<<std::endl;
	if(mode<0){if(mode==-2) retval=-2; return;}
	else if(mode ==3){
		markHotStart();
	}
    
    map_topology();
    map_duals();
    if(mode==-2){retval=-2; return;}
    volprob_->value=0;
    volprob_->active_size = fsize + csize;
    volprob_->psize = szunfxd + data->ndemands*sznz;
    volprob_->active_psize = szunfxd;
    set_start();
    retval = volprob_->solve(*this, true);
	//if(!in_strong_branch)std::cout<<"volsol: "<<volprob_->value<<" "<<min_lower_bound<<" "<<volprob_->iter()<<std::endl;

    if(volprob_->value< min_lower_bound && in_strong_branch){
     	volprob_->value = min_lower_bound;
    }

    translate_sol();
    if(mode ==3){
		unmarkHotStart();
	}
	 
 	
    // extract the solutio
    // the lower bound on the objective value
}

//---------------------------------------------------------------------------

void OsiVolSolverInterface::test_opposites(BCP_vec<int>& changed_pos, BCP_vec<double>& new_bd, int * yfix,
										 const BCP_vec<BCP_var*>& vars,  double bestsolv){
    markHotStart();
    double currlb = volprob_->value;
    mode=0;
    std::deque<int> totest;
    for(int a=0;a<narcs;++a){
		collb[a] = yfix[a]==-1 ? vars[a]->lb() : yfix[a];
		colub[a] = yfix[a]==-1 ? vars[a]->ub() : yfix[a];
		/*for(int k=0;k<ndemands;++k){
			collb[narcs+k*narcs+a] = vars[narcs+k*narcs+a]->lb();
			colub[narcs+k*narcs+a] = vars[narcs+k*narcs+a]->ub();
		}*/
		if(yfix[a]==-1 && (solution[a]<0.1 || solution[a]>0.9))totest.push_back(a);
	}
	
    volprob_->parm.ascent_check_invl = 50;
	volprob_->parm.ascent_first_check = 50;

    int arc, ret;
    while(!totest.empty()){
		arc = totest.back();
		totest.pop_back();
 		if(solution[arc]>0.9){
 			collb[arc]=0;
 			colub[arc]=0;
		}else{
			collb[arc]=1;
 			colub[arc]=1;
		}
		mode=0;
  		recheck_collct=true;
  		map_topology();
    	map_duals();
		ret=1;
 		if(mode==-2){ ret=-1;}
 		else{
 			volprob_->active_size = fsize + csize;
    		volprob_->psize = szunfxd + data->ndemands*sznz;
    		volprob_->active_psize = szunfxd;
    		set_start();
    		ret = volprob_->solve(*this, true);
 		}
    	
 		
		if(ret<0 ||  volprob_->value<=0 || volprob_->value>bestsolv){
 			 ret = -1;
		}
		
		//std::cout<<arc<<" OsiVolSolverInterface::test_opposites "<<volprob_->value<<std::endl;
 		if( ret<0 ){
			if(solution[arc]>0.9){
				//std::cout<<arc<<" OsiVolSolverInterface::test_opposites FIX 1 "<<solution[arc]<<std::endl;
				yfix[arc]=1;
				collb[arc]=1;
 				colub[arc]=1;
 				changed_pos.push_back(arc);
				new_bd.push_back(1.0);
				new_bd.push_back(1.0);
			}else{
				//std::cout<<arc<<" OsiVolSolverInterface::test_opposites FIX 0 "<<solution[arc]<<std::endl;
				yfix[arc]=0;
				collb[arc]=0;
 				colub[arc]=0;
 				changed_pos.push_back(arc);
				new_bd.push_back(0.0);
				new_bd.push_back(0.0);
			}
		}else{
			collb[arc]= vars[arc]->lb();
 			colub[arc]= vars[arc]->ub();
		}
		 
	}
	has_checked=true;
	volprob_->value = currlb;
    unmarkHotStart();
    
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// VOLUME USER HOOKS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


int
OsiVolSolverInterface::addVI(int iter,double lcost, const VOL_dvector& xstar,
          const VOL_dvector& x, VOL_dvector& dstar,  VOL_dvector& dualu, VOL_dvector& dual_lb, VOL_dvector& dual_ub,
          VOL_dvector& rc, VOL_dvector& h, int & actvSSz){
    //return 0;
    bool letgen=false;
    if(mode == 0) return 0;
    
    if(lcost >= VItt ){
        VItt = lcost;
        letgen = true;
        int arc;
        for(int a=szunfxd; a--; ){ 
        	arc = nz_arcs[a];
        	rc_[arc] = rc[a];
        	for(int k=ndemands; k--; ){
				solution[narcs+k*narcs+arc] = x[szunfxd + k*sznz + a];
			}	
        }
        for(int a=szunfxd; a<sznz; ++a){ 
        	arc = nz_arcs[a];
        	for(int k=ndemands; k--; ){
				solution[narcs+k*narcs+arc] = x[szunfxd + k*sznz + a];
				//if(arc==52 || arc==8)std::cout<<"better "<<narcs+k*narcs+arc<<" "<<x[szunfxd + k*sznz + a]<<" "
				//<<rc[szunfxd + k*sznz + a]<<" ub: "<<colub[narcs+k*narcs+arc]/data->d_k[k].quantity<<std::endl;
 			}
        }
        
    }

    if( (mode == 1) && iter>0 && iter%intvlVI==0){
        translate_primal(xstar);
        int num_new_sets=ss_manager->cutset_generation_main( yhit, VItopo, false);
        VIub=-1e31;
    }
    
    if(numrows_>=maxNumrows_) return 0;

    if(letgen && iter>=minIterVI){

        int num_covers = cover_manager->cover_generation_main(xstar.v, x.v, &ss_manager->sets, numrows_, maxNumrows_);
        num_covers += cover_manager->mincard_generation_main(xstar.v, x.v, &ss_manager->sets, numrows_+num_covers, maxNumrows_);
        cover_manager->add_cover_vi(num_covers, actv, actvSSz, h.v, dstar.v, dualu.v, dual_lb.v,  dual_ub.v );
        if(num_covers>0){
        	numrows_ +=  num_covers;
        	cover_manager->reposition_covers(num_covers);
            //std::cout<<std::setprecision(10)<<"iter: "<<iter<<" added cuts "<<num_covers<<" L: "<<VItt<<std::endl;
        	/*Cover * c = cover_manager->covers.begin ;
        	for(int i =cover_manager->covers.sizeOfCollection; i--; ){
        		std::cout<<"vi "<<c->id_vi<<std::endl;
        		c = c->next;
			}*/
        } 

    }
    return 0;
}


//-----------------------------------------------------------------------

int
OsiVolSolverInterface::removeVI( int & actvSSz, VOL_dvector& pstarv, VOL_dvector& dstaru,  VOL_dvector& dualu){

	int oldsz  = actvSSz;
    cover_manager->covers.desactvCover(cover_manager->lim_to_remv, actv, actvSSz, cover_manager->num_actv, pstarv.v, dstaru.v, dualu.v);
    localc_manager->locals.desactvLocalc(localc_manager->lim_to_remv, actv, actvSSz, localc_manager->num_actv, pstarv.v, dstaru.v, dualu.v);
    globalc_manager->globals.desactvGlobalc(globalc_manager->lim_to_remv, actv, actvSSz, globalc_manager->num_actv, pstarv.v, dstaru.v, dualu.v);

    if(oldsz == actvSSz) return 0;
    	
 	int id;
	int sz = cover_manager->num_actv;
	int sz2 = sz+localc_manager->num_actv;
	int sz3 = sz2+globalc_manager->num_actv;
 	std::vector<TrioF> collect(sz3);
	
	Cover* vic = cover_manager->covers.begin;
 	for(int i=0;i<sz;++i){
 		//std::cout<<"globalc: "<<i<<" id: "<<vigc->id_vi<<" srnb: "<<vigc->serial_nmbr<<std::endl;		
 		id = actv[vic->id_vi];
		collect[i].fst = pstarv[id];	
		collect[i].snd = dstaru[id];	
		collect[i].trd = dualu[id];	
		vic = vic->next;
	}
	
	LocalCut* vilc = localc_manager->locals.begin;
 	for(int i=sz;i<sz2;++i){
 		//std::cout<<"globalc: "<<i<<" id: "<<vigc->id_vi<<" srnb: "<<vigc->serial_nmbr<<std::endl;		
		id = actv[vilc->id_vi];
		collect[i].fst = pstarv[id];	
		collect[i].snd = dstaru[id];	
		collect[i].trd = dualu[id];		
		vilc = vilc->next;
	}
	
	GlobalCut* vigc = globalc_manager->globals.begin;
 	for(int i=sz2;i<sz3;++i){
 		//std::cout<<"globalc: "<<i<<" id: "<<vigc->id_vi<<" srnb: "<<vigc->serial_nmbr<<std::endl;		
		id = actv[vigc->id_vi];
		collect[i].fst = pstarv[id];	
		collect[i].snd = dstaru[id];	
		collect[i].trd = dualu[id];				
		vigc = vigc->next;
	}
	//std::cout<<"first: "<<fidx+sz3<<" sz: "<<actvSSz<<std::endl;
	//redo identifications
	vic = cover_manager->covers.begin;
 	for(int i=0;i<sz;++i){
 		id = fsize+i;
		actv[vic->id_vi] = id;
 		pstarv[id] = collect[i].fst;	
		dstaru[id] = collect[i].snd;	
		dualu[id] = collect[i].trd;	
		//std::cout<<"cover: "<<i<<" id: "<< actv[vic->id_vi]<<" -> "<<pstarv[actv[vic->id_vi]]<<" "<<dstaru[actv[vic->id_vi]]<<" "<<dualu[actv[vic->id_vi]]<<std::endl;			
 		vic = vic->next;
	}
 	vilc = localc_manager->locals.begin;
	for(int i=sz;i<sz2;++i){
		id = fsize+i;
		actv[vilc->id_vi] = id;
 		pstarv[id] = collect[i].fst;	
		dstaru[id] = collect[i].snd;	
		dualu[id] = collect[i].trd;
		//std::cout<<"local: "<<i<<" id: "<< actv[vilc->id_vi]<<" -> "<<pstarv[actv[vilc->id_vi]]<<" "<<dstaru[actv[vilc->id_vi]]<<" "<<dualu[actv[vilc->id_vi]]<<std::endl;				
		vilc = vilc->next;
	}
	vigc = globalc_manager->globals.begin;
	for(int i=sz2;i<sz3;++i){
		id = fsize+i;
		actv[vigc->id_vi] = id;
 		pstarv[id] = collect[i].fst;	
		dstaru[id] = collect[i].snd;	
		dualu[id] = collect[i].trd;	
		//std::cout<<"global: "<<i<<" id: "<< actv[vigc->id_vi]<<" -> "<<pstarv[actv[vigc->id_vi]]<<" "<<dstaru[actv[vigc->id_vi]]<<" "<<dualu[actv[vigc->id_vi]]<<std::endl;			 				
		vigc = vigc->next;
	}
	collect.clear();
	if(fsize+sz3 != actvSSz){ std::cout<<"cont: "<<fsize+sz3<<" sz: "<<actvSSz<<std::endl; abort();}
	//abort();
    return 0;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::compute_rc(const VOL_dvector& dualu, VOL_dvector& rc, int actvSSz){
    const Arc* item;
    int arc, id;
    double dvali, dvalj;
    for(int a=szunfxd; a--;){
        arc = nz_arcs[a];
        item = &data->arcs[arc];
        rc[a] = item->f;
        //addrc[a] =0;
        for(int k=ndemands; k--; ){
        	id = actv[k*nnodes + item->j-1];
        	dvalj = id>=0 ? dualu[id]: 0.0;
        	id = actv[k*nnodes + item->i-1];
        	dvali = id>=0 ? dualu[id]: 0.0;
            rc[szunfxd + k*sznz + a] = item->c[k]*data->d_k[k].quantity -dvalj + dvali;
        }
    }
    for(int a=szunfxd; a<sznz;++a){
        arc = nz_arcs[a];
        item = &data->arcs[arc];
        for(int k=ndemands; k-- ;){
        	id = actv[k*nnodes + item->j-1];
        	dvalj = id>=0 ? dualu[id]: 0.0;
        	id = actv[k*nnodes + item->i-1];
        	dvali = id>=0 ? dualu[id]: 0.0;
            rc[szunfxd + k*sznz + a] = item->c[k]*data->d_k[k].quantity -dvalj + dvali;
        }
    }

    B0=0;
    for(int k=0; k<data->ndemands; ++k){
        B0 += ( dualu[actv[k*nnodes + data->d_k[k].D-1]] - dualu[actv[k*nnodes + data->d_k[k].O-1]]);
        
    }
    cover_manager->compute_cover_rc( dualu.v, actv,  actvSSz, rc.v,   B0);
    localc_manager->compute_localc_rc( dualu.v, actv,  actvSSz, rc.v,   B0);
    globalc_manager->compute_cover_rc( dualu.v, actv,  actvSSz, rc.v,   B0);
    //double b = B0;
    //cover_manager->compute_cover_rc( dual.v,  actv,  actvSSz, addrc, b);
    return 0;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::compute_sg(const VOL_dvector& x, int actvSSz, VOL_dvector& v){
    //std::cout<<"compute_sg"<<std::endl;
    for(int n=actvSSz; n--; )
    	v[n] = 0;
    	
    int arc, basek, basef;
    int id;
    const Arc* item; const Demand* itemd;
    for(int k=ndemands; k--; ){
        itemd = &data->d_k[k];
        basek = k*nnodes;
        id = actv[basek+ itemd->O-1];
        if(id>=0)v[id] -= 1;
        id=actv[basek + itemd->D-1];
        if(id>=0)v[id] += 1;
        for(int a=0; a<sznz; ++a){
            arc = nz_arcs[a];
            item = &data->arcs[arc];
            basef = szunfxd+k*sznz;
            id = actv[basek + item->i-1];
            if(id>=0)v[id] += x[basef+ a];
            id = actv[basek + item->j-1];
            if(id>=0)v[id] -= x[basef+ a];
        }
    }
    cover_manager->compute_cover_sg( x.v, actv,  actvSSz, v.v);
    localc_manager->compute_localc_sg( x.v, actv,  actvSSz, v.v);
    globalc_manager->compute_cover_sg( x.v, actv,  actvSSz, v.v);

    return 0;
}


//-----------------------------------------------------------------------

int
OsiVolSolverInterface::solve_subproblem(const VOL_dvector& xstar,
                                        const VOL_dvector& dualu,  VOL_dvector& rc,
                                        double& lcost, VOL_dvector& x,
                                        double& pcost){
    int arc;
    double cost_a, ratio;
    pcost= 0;
    lcost =0;
    for(int a=szunfxd; a--; ){
    	//if(xstar[a]>1)std::cout<<"WHAT? "<<xstar[a]<<std::endl;
        arc = nz_arcs[a];
        cost_a = knapsack(a, rc.v, x.v);
        rc[a] += cost_a;
        if(cost_a > volprob_->parm.dual_limit) return -1; 
        //addrc[a]+= rc[a];
        //std::cout<<a<<", "<<arc<<" addrc: "<<addrc[a]<<std::endl;
        //if(addrc[a]<1e-10 && addrc[a]>-1e-10) addrc[a] = 0;
        /*if( addrc[a] < 0.0){
            x[a] =1.0;
        }else x[a] = 0.0;*/
    }
    for(int a=szunfxd; a<sznz;++a){
        arc = nz_arcs[a];
        cost_a = knapsack(a, rc.v, x.v);
        lcost += cost_a + data->arcs[arc].f;
        if(cost_a > volprob_->parm.dual_limit) return -1; 
        //pcost += data->arcs[arc].f;
        //for(int k=ndemands; k--; )
         //   pcost += data->arcs[arc].c[k] * x[szunfxd + k*sznz + a];
    }
    
    return 0;
}

//---------------------------------------------------------------------------

int
OsiVolSolverInterface::resolve_subproblem(const VOL_dvector& dualu, VOL_dvector& rc,
                       double& lcost,
                       VOL_dvector& x,double& pcost){
    
    int arc;
    lcost += B0;
    for(int a=szunfxd; a--; ){
        arc = nz_arcs[a];
        //if(rc[a]<1e-10 && rc[a]>-1e-10) rc[a] = 0;
        if( rc[a]<0 /*|| (rc[a]==0 && x[a]==1.0)*/){
            lcost += rc[a];
            x[a]=1.0;
            
        }else{
            //if(x[a]==1.0){ std::cout<<"rarc: "<<a<<" "<<arc<<" rc: "<<rc[a]<<" y: "<<x[a]<<std::endl; abort();}
            x[a]=0.0;
            for(int k=0; k<ndemands; ++k)
                x[szunfxd + k*sznz + a]=0.0;
            
        }
    }
    //std::cout<<"lcost: "<<lcost<<std::endl;
    mark_topo( x,lcost);
    return 0;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::additional_settings(int iter, double& lcost, VOL_dvector& dualu, VOL_dvector& rc, VOL_dvector& h,
                          VOL_dvector& x, const VOL_dvector& xhist, int actvSSz){
   return 0;
   /* if(cover_manager->covers.sizeOfCollection==0 || iter==0){
        cover_manager->compute_cover_rc( dualu.v, actv,  actvSSz, rc.v,   B0);
        return 0;
    }*/
    //double ret=0;
	
	
    //cover_manager->recompute_mult_neg( dualu.v, addrc, rc.v, x.v, actv, actvSSz);
    //cover_manager->recompute_mult_pos( dualu.v, h.v, addrc, x.v, actv);
    /*int sz = cover_manager->num_actv;
    int index;
    Cover *vi = cover_manager->covers.begin;
    for(;sz--;){
        index = actv[vi->id_vi];
        B0 +=  dualu[index]*vi->get_total_rhs();
        vi = vi->next;
    }*/
    /*for(int a=szunfxd; a--; ){
    	//rc[a] = addrc[a];
        if(rc[a]<1e-5 && rc[a]>-1e-5) rc[a] = 0;
        if(rc[a]==0 ){
            ret = arc_dg_imp(a, x.v, h.v , actvSSz);
            if(ret>0 && xhist[a]>=0.0){
                //std::cout<<"yh: "<<xhist[a]<<std::endl;
                //std::cout<<"opaa zero to one: "<<a<<" rc: "<<addrc[a]<<" rc+: "<<rc[a]<<" y: "<<x[a]<<std::endl;
                x[a] = 1.0 ;
            }else if(ret<0 && xhist[a]<=0.0)  x[a] = 0.0 ;
        }
    }*/
    //return 0;
}

//---------------------------------------------------------------------------
//  solving methods
//---------------------------------------------------------------------------

double
OsiVolSolverInterface::knapsack(int a, const double * rc, double* x){
    double kpsack =0;
    double fillUp =0;
    double rcost, flow, xval;
    int arc = nz_arcs[a];
    int basex;
    std::list<HeapCell> heap;
    //get reduced cost for each commodity in arc e
    for(int k=ndemands; k--; ){
    	rcost = rc[szunfxd + k*sznz + a];
        if(rcost<0.0 && colub[narcs+k*narcs+arc]>0.0){
            heap.push_back(HeapCell(k, double(rcost/data->d_k[k].quantity)));
       }
        xval= collb[narcs+k*narcs+arc];
        if(xval>0){
         	fillUp += xval;
        	xval /= data->d_k[k].quantity;
        	if(xval>1){   xval = 1.0; }
        	kpsack += rcost * xval;
        	x[szunfxd + k*sznz + a] = xval;
        }else x[szunfxd + k*sznz + a] = 0;
    }
    
    double capa = data->arcs[arc].capa;
    if(fillUp >= capa || heap.empty() ){ 
    	if(fillUp > capa) kpsack = volprob_->parm.dual_limit+1;
    	heap.clear(); return kpsack;
    }
    heap.sort(comp());
    //std::stable_sort(heap.begin(), heap.end(), comp());
     
    int comm;
    while(heap.size()>0){ 
        comm = heap.back().k;
        basex = szunfxd + comm*sznz;
        if(fillUp < capa){
        	//if(colub[narcs+comm*narcs+arc]<data->arcs[arc].b[comm])
        	// std::cout<<narcs+comm*narcs+arc<<" confirm: "<<arc+1<<" comm: "<<comm+1<<" "<<colub[narcs+comm*narcs+arc]<<std::endl;
        	flow = std::min((capa - fillUp), (colub[narcs+comm*narcs+arc]-collb[narcs+comm*narcs+arc]));
        	if(flow<0)flow=0;
        	xval = flow/data->d_k[comm].quantity;
        	if(xval>1){  xval =1.0; }

            x[basex + a] += xval;
            fillUp += flow;
            kpsack += rc[basex + a] * xval;
            
        } 
        heap.pop_back();
    }
    return kpsack;
}

//---------------------------------------------------------------------------

int
OsiVolSolverInterface::heuristics(const VOL_problem& p,
                                  const VOL_dvector& x, double& heur_val){ return 0;}
                                  

//---------------------------------------------------------------------------

double
OsiVolSolverInterface::arc_dg_imp(int a, const double * xy, const double * h, int actvSSz){
    double dg=0;
    int arc = nz_arcs[a];
    for(int k=0; k<ndemands; ++k){
        //std::cout<<narcs+k*narcs+ a<<" "<<k*nnodes + data->arcs[a].i-1<<" "<< k*nnodes + data->arcs[a].j-1<<std::endl;
        dg += h[actv[k*nnodes + data->arcs[arc].i-1]]*xy[szunfxd + k*sznz + a];
        dg += -h[actv[k*nnodes + data->arcs[arc].j-1]]*xy[szunfxd + k*sznz + a];
    }
    
    dg += cover_manager->arc_dg_imp(a, xy, h, actv,  actvSSz);
	dg += localc_manager->arc_dg_imp(a, xy, h, actv,  actvSSz);
	dg += globalc_manager->arc_dg_imp(a, xy, h, actv,  actvSSz);
    //if(dg>0) std::cout<<"arc: "<<arc<<" dg: "<<dg<<std::endl;
    return dg;
}

//---------------------------------------------------------------------------
//  auxiliary methods
//---------------------------------------------------------------------------


int
OsiVolSolverInterface::mark_topo( VOL_dvector& x, double lcost){
    
    //std::cout<<"Lt: "<<lcost<<std::endl;
    if(lcost<=VIub || mode > 1 || mode==0) return 0;
    int arc;
    VIub =lcost;
    for(int a=szunfxd; a--; ){
        arc = nz_arcs[a];
        VItopo[arc] = x[a];
    }
    return 0;
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::translate_primal( const VOL_dvector& xhist){
    int arc;
    for(int a=szunfxd; a--;){
        arc = nz_arcs[a];
        yhit[arc] = xhist[a];
    }
}

//-----------------------------------------------------------------------

void 
OsiVolSolverInterface::add_external_cover(const std::deque<Pair2>& c){
    if(numrows_>=maxNumrows_) return ;

	int ret = cover_manager->add_external_cover(c, numrows_, maxNumrows_);
	if(ret){
	
		cover_manager->add_cover_vi(1, actv, volprob_->active_size, volprob_->viol.v, volprob_->dsol.v,
								volprob_->dsol.v, volprob_->dual_lb.v,  volprob_->dual_ub.v );       
		cover_manager->reposition_covers(1);	
		++numrows_;
	}
	
}

//-----------------------------------------------------------------------

int 
OsiVolSolverInterface::add_external_localc( const int * y, const double * sol, int sz, int type){ 
    if(numrows_>=maxNumrows_) return 0;

	int ret=0;
	switch(type){
		case 0:{
			ret = localc_manager->localc0_generation_main(solution,  y, sz, numrows_);
			break;
		}
		case 1:{
			ret = localc_manager->localc1_generation_main(volprob_->value, upper_bound, solution, 
														y, rc_, numrows_, maxNumrows_);
			break;
		}
		case 2:{
			ret = localc_manager->localc2_generation_main(solution, sol, sz, numrows_);
			break;
		}
		default: break;
	}
		
		
	if(ret){
		//std::cout<<"add local: "<<std::endl;
		localc_manager->add_local_vi(ret, actv, volprob_->active_size, dual, lhs_, volprob_->viol.v,
								volprob_->dsol.v, volprob_->dual_lb.v,  volprob_->dual_ub.v );   
		localc_manager->reposition_locals(ret);	
		numrows_ += ret;
	}
	return ret;
}


//-----------------------------------------------------------------------

int 
OsiVolSolverInterface::add_external_globalc( const int * y, int cont0){
	
	if(cont0==0) return 0;     
 	int ret = globalc_manager->globalc_generation_main(solution , y, cont0, numrows_);
 
	if(numrows_>=maxNumrows_) return 0;
	
	if(ret){
		GlobalCut* gloc = globalc_manager->globals.track.back();
		globalc_manager->globals.insert_end(gloc);
		globalc_manager->num_actv += ret;
		globalc_manager->add_global_vi(ret, actv, volprob_->active_size, dual, lhs_, volprob_->viol.v, 
								volprob_->dsol.v, volprob_->dual_lb.v,  volprob_->dual_ub.v );   
		globalc_manager->reposition_globals(ret);	
		numrows_ += ret;
		 
	}
	return ret;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//Solution
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

void 
OsiVolSolverInterface::translate_sol(){
    
    CoinFillN(solution, narcs, 0.0);

    int arc;
    double addvalue;
    double iters = volprob_->parm.maxsgriters; 
    double coeff = iters>1? volprob_->iter()/iters: 0.0;
    if(coeff>0.7) coeff=0.7;
	//if(!in_strong_branch)std::cout<<std::setprecision(10)<<"OsiVolSolverInterface::translate_sol: "<<volprob_->value<<" numrows: "<<numrows_<<" iters: "
	//<<volprob_->iter()<<"/"<<volprob_->parm.maxsgriters<<" coef: "<<coeff<<std::endl;

	for(int a=szunfxd; a--;){
		arc = nz_arcs[a];
		if(HotStartSet) solution[arc] = (coeff*volprob_->psol[a] + (1.0-coeff)*HotStart_->primal[arc]);
		else solution[arc] = volprob_->psol[a];
		
		//std::cout<<arc<<" "<<volprob_->psol[a]<<std::endl;// /double(volprob_->iter())<<std::endl;
		
	}
	for(int a=szunfxd; a<sznz;++a){
		arc = nz_arcs[a];
		solution[arc] = 1;
		//for(int k=ndemands; k--; )
			//solution[narcs+k*narcs+arc] = volprob_->psol[szunfxd + k*sznz + a];
	}
	translate_dualsol();
    
    const Arc* item;
    for(int a=sznz; a--;){
		arc = nz_arcs[a];
		item = &data->arcs[arc];
		for(int k=ndemands; k--; )
			rc_[narcs+k*narcs+arc] = item->c[k]*data->d_k[k].quantity - dual[k*nnodes + item->j-1] + dual[k*nnodes + item->i-1];
	}
}

//-----------------------------------------------------------------------

void 
OsiVolSolverInterface::translate_dualsol(){
	int idx;
	for(int i=ndemands*nnodes; i-- ;){
        idx = actv[i];
        if(idx>=0){
            //if(i>=ndemands*nnodes) std::cout<<i<<": "<<volprob_->dsol[idx]<<std::endl;
            dual[i] = volprob_->dsol[idx];
            lhs_[i] = volprob_->viol[idx];
        }else{
            dual[i] =0;
            lhs_[i] =0;
        }
    }
    int sz = cover_manager->covers.sizeOfCollection;
    Cover* vi = cover_manager->covers.end;    
	for(;sz--;){
		idx = actv[vi->id_vi];
		if(idx>=0){
            dual[vi->id_vi] = volprob_->dsol[idx];
            lhs_[vi->id_vi] = volprob_->viol[idx];
        }else{
            dual[vi->id_vi] =0;
            lhs_[vi->id_vi] =0;
        }
        //std::cout<<idx<<" sol "<<vi->serial_nmbr<<" "<<dual[vi->id_vi]<<std::endl;
		vi = vi->prev;
	}
	sz = localc_manager->locals.sizeOfCollection;
    LocalCut* vilc = localc_manager->locals.end;    
	for(;sz--;){
		//std::cout<<"ok "<<sz<<" "<<vi->id_vi<<std::endl;
		idx = actv[vilc->id_vi];
		if(idx>=0){
            dual[vilc->id_vi] = volprob_->dsol[idx];
            lhs_[vilc->id_vi] = volprob_->viol[idx];
        }else{
            dual[vilc->id_vi] =0;
            lhs_[vilc->id_vi] =0;
        }
		//std::cout<<idx<<" sol "<<vilc->serial_nmbr<<" "<<dual[vilc->id_vi]<<std::endl;
		vilc = vilc->prev;
		//std::cout<<"vi next: "<<vi<<std::endl;
	}
	sz = globalc_manager->globals.sizeOfCollection;
    GlobalCut* vigc = globalc_manager->globals.end;    
	for(;sz--;){
		//std::cout<<"ok "<<sz<<" "<<vi->id_vi<<std::endl;
		idx = actv[vigc->id_vi];
		if(idx>=0){
            dual[vigc->id_vi] = volprob_->dsol[idx];
            lhs_[vigc->id_vi] = volprob_->viol[idx];
        }else{
            dual[vigc->id_vi] =0;
            lhs_[vigc->id_vi] =0;
        }
		//std::cout<<idx<<" sol "<<vigc->serial_nmbr<<" "<<dual[vigc->id_vi]<<std::endl;
		vigc = vigc->prev;
		//std::cout<<"vi next: "<<vi<<std::endl;
	}
}

//-----------------------------------------------------------------------

void 
OsiVolSolverInterface::reset_dualsol(const std::map<int, double>& dual_map){
 	int idx;
 	int fx =  ndemands*nnodes;
 	std::map<int, double>::const_iterator it;
	for(int i=fx; i-- ;){
        idx = actv[i];
        if(idx>=0){
        	it = dual_map.find(i);
            dual[i] =  it != dual_map.end() ? it->second : 0.0 ;
        }else{
            dual[i] =0;
        }
    }

   int sz = cover_manager->covers.sizeOfCollection;
    Cover* vi = cover_manager->covers.end;    
	for(; sz--;){
		idx = actv[vi->id_vi];
		if(idx>=0){
            it = dual_map.find(fx+vi->serial_nmbr);
            //if(it != dual_map.end()) std::cout<<"restdual: "<<vi->serial_nmbr<<" "<<it->second<<std::endl;
            dual[vi->id_vi] =  it != dual_map.end() ? it->second : 0.0 ;
        }else{
            dual[vi->id_vi] =0;
        }
		vi = vi->prev;
	}
	sz = localc_manager->locals.sizeOfCollection;
    LocalCut* vilc = localc_manager->locals.end;    
	for(; sz--;){
		idx = actv[vilc->id_vi];
		if(idx>=0){
            it = dual_map.find(fx+vilc->serial_nmbr);
            //if(it != dual_map.end()) std::cout<<"restdual: "<<vilc->serial_nmbr<<" "<<it->second<<std::endl;
            dual[vilc->id_vi] =  it != dual_map.end() ? it->second : 0.0 ;
        }else{
            dual[vilc->id_vi] = 0;
        }
		vilc = vilc->prev;
	}
	sz = globalc_manager->globals.sizeOfCollection;
    GlobalCut* vigc = globalc_manager->globals.end;    
	for(; sz--;){
		idx = actv[vigc->id_vi];
		if(idx>=0){
            it = dual_map.find(fx+vigc->serial_nmbr);
            dual[vigc->id_vi] =  it != dual_map.end() ? it->second : 0.0 ;
        }else{
            dual[vigc->id_vi] = 0;
        }
		vigc = vigc->prev;
	}

}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

bool 
OsiVolSolverInterface::isProvenPrimalInfeasible()const{
	if(mode==-2) return true;
	if(retval==-1 || volprob_->value<0 || volprob_->value>(volprob_->parm.dual_limit + 0.0001)
		|| (isPrimalObjectiveLimitReached() && !has_sol)){
		return true;
	}else return false;
}
    
//-----------------------------------------------------------------------

bool 
OsiVolSolverInterface::isPrimalObjectiveLimitReached() const{
	if((retval==0)){
		if(volprob_->value>(volprob_->parm.dual_limit + 0.0001)) return false;
		return (volprob_->value>(upper_bound + 0.0001))? true : false;
	}else return false;
}
 
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//loaders
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

void
OsiVolSolverInterface::loadProblem(const int numcols, const int numrows,
                                   const int* start, const int* index,
                                   const double* value,
                                   const double* collb_, const double* colub_,
                                   const double* obj,
                                   const double* rowlb_, const double* rowub_){
    
    //std::cout<<"OsiVolSolverInterface::loadProblem nrows: "<<numrows<<std::endl;
    unmarkHotStart();
    OsiVolAuxInfo * auxinfo  = static_cast<OsiVolAuxInfo*>(OsiSolverInterface::getApplicationData());
    
    maxNumVI =auxinfo->maxNumVI;
    intvlVI = auxinfo->intvlVI;
	minIterVI = auxinfo->minIterVI;
	
    nnodes =  data->nnodes;
    ndemands = data->ndemands;
    narcs = data->narcs;
    
    volprob_->value =0;
    retval = 0;

    maxNumcols_ = numcols_ = numcols;
    numrows_ = numrows + cover_manager->covers.sizeOfCollection;
    numrows_ += localc_manager->locals.sizeOfCollection;
    numrows_ += globalc_manager->globals.sizeOfCollection;
    maxNumrows_ = numrows_+maxNumVI;
    if(numrows_> auxinfo->maxPos){
    	std::cout<<"OsiVolSolverInterface::loadProblem PROBLEM "<<numrows_<<" > "<<auxinfo->maxPos<<std::endl;
    	abort();
    }
    if(maxNumrows_ > auxinfo->maxPos)
    	maxNumrows_ =  auxinfo->maxPos;
    
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::deleteRows(const int num, const int * rowIndices){ 
	//std::cout<<"deleteRows "<<numrows_<<std::endl;
	numrows_ -= num;
 	Cover* vi;
	while(!cover_manager->purgbl.empty()){
		vi = cover_manager->purgbl.back();
		cover_manager->covers.remove_nodel(vi);
		cover_manager->purgbl.pop_back();
	}

	LocalCut* vilc;
	while(!localc_manager->purgbl.empty()){
		vilc = localc_manager->purgbl.back();
		//std::cout<<"outlocal "<<vilc->serial_nmbr<<std::endl;
		localc_manager->locals.remove_nodel(vilc);
		localc_manager->purgbl.pop_back();
	}
	
	rebuild_collections();
	volprob_->active_size = fsize+localc_manager->num_actv + cover_manager->num_actv + globalc_manager->num_actv;
	//abort();
}

//-----------------------------------------------------------------------------

bool 
OsiVolSolverInterface::deleteGlobalRows(){ 
 
	std::deque<GlobalCut *> purgbl;
 	GlobalCut * gloc = globalc_manager->globals.begin;
	for(int i = globalc_manager->globals.sizeOfCollection; i--;){
		if(gloc->purgbl){
 			purgbl.push_back(gloc);
			gloc->purgbl=false;
		} 
		gloc=gloc->next;
	}
	bool ret = purgbl.empty();
	while(!purgbl.empty()){
		gloc =purgbl.back();
		//std::cout<<"outglobal "<<gloc->serial_nmbr<<std::endl;
		--numrows_;
		globalc_manager->globals.remove_nodel(gloc);
		purgbl.pop_back();
	}
	return !ret;
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::rebuild_collections(){
	int fidx = nnodes*ndemands;

	Cover* vi = cover_manager->covers.end;
	int sz = cover_manager->covers.sizeOfCollection;
	int sz2 = sz+ localc_manager->locals.sizeOfCollection;
	int sz3 = sz2+globalc_manager->globals.sizeOfCollection;
 	std::vector<Trio1> collect(sz3);
	for(int i=sz;i--;){
		collect[i].fst = actv[vi->id_vi];
		collect[i].snd = dual[vi->id_vi];
		collect[i].trd = lhs_[vi->id_vi];
		vi = vi->prev;
	}
	LocalCut* vilc = localc_manager->locals.end;
 	for(int i=sz2-1;i>=sz;--i){
		collect[i].fst = actv[vilc->id_vi];
		collect[i].snd = dual[vilc->id_vi];
		collect[i].trd = lhs_[vilc->id_vi];
		vilc = vilc->prev;
	}
	GlobalCut* vigc = globalc_manager->globals.end;
 	for(int i=sz3-1;i>=sz2;--i){
 		//std::cout<<"globalc: "<<i<<" id: "<<vigc->id_vi<<" srnb: "<<vigc->serial_nmbr<<std::endl;		
		collect[i].fst = actv[vigc->id_vi];
		collect[i].snd = dual[vigc->id_vi];
		collect[i].trd = lhs_[vigc->id_vi];
		vigc = vigc->prev;
	}

	vi = cover_manager->covers.end;
	for(int i=sz;i--;){
		//std::cout<<" new indice: "<<i+fidx<<" id: "<<vi->id_vi<<" srnb: "<<vi->serial_nmbr<<std::endl;	
		if(vi->id_vi>=numrows_)	 actv[vi->id_vi] = -1;
		vi->id_vi = i+fidx;
		actv[vi->id_vi] = collect[i].fst;
		dual[vi->id_vi] = collect[i].snd;
		lhs_[vi->id_vi] = collect[i].trd;	 
		vi = vi->prev;
	}
 	vilc = localc_manager->locals.end;
	for(int i=sz2-1;i>=sz;--i){
		//std::cout<<" new indicelc: "<<i+fidx<<" id: "<<vilc->id_vi<<" srnb: "<<vilc->serial_nmbr<<std::endl;		
		if(vilc->id_vi>=numrows_) actv[vilc->id_vi] = -1;
		vilc->id_vi = i+fidx;
		actv[vilc->id_vi] = collect[i].fst;
		dual[vilc->id_vi] = collect[i].snd;
		lhs_[vilc->id_vi] = collect[i].trd;	
		vilc = vilc->prev;
	}
	vigc = globalc_manager->globals.end;
	for(int i=sz3-1;i>=sz2;--i){
		//std::cout<<" new indicegc: "<<i+fidx<<" id: "<<vigc->id_vi<<" srnb: "<<vigc->serial_nmbr<<std::endl;		
		if(vigc->id_vi>=numrows_) actv[vigc->id_vi] = -1;
		vigc->id_vi = i+fidx;
		actv[vigc->id_vi] = collect[i].fst;
		dual[vigc->id_vi] = collect[i].snd;
		lhs_[vigc->id_vi] = collect[i].trd;	
		vigc = vigc->prev;
	}
	collect.clear();
	cover_manager->num_actv =  cover_manager->covers.sizeOfCollection;
	localc_manager->num_actv =  localc_manager->locals.sizeOfCollection;
	globalc_manager->num_actv =  globalc_manager->globals.sizeOfCollection;
	/*vi = cover_manager->covers.end;
 	for(;sz--;){
		std::cout<<"vi serial: "<<vi->serial_nmbr<<" id: "<<vi->id_vi<<" sol: "<<dual[vi->id_vi]<<std::endl;
		//std::cout<<"sol "<<fidx+sz<<" "<<vi->id_vi<<" "<<dual[fidx+sz]<<std::endl;
		vi = vi->prev;
	}
	vilc = localc_manager->locals.end;
	for(;sz2--;){
		std::cout<<"vi serial: "<<vilc->serial_nmbr<<" id: "<<vilc->id_vi<<" sol: "<<dual[vilc->id_vi]<<std::endl;
		//std::cout<<"sol "<<fidx+sz<<" "<<vi->id_vi<<" "<<dual[fidx+sz]<<std::endl;
		vilc = vilc->prev;
	}*/
	
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::rowRimAllocator_()
{
    /*rowub = new double[maxNumrows_];
    rowlb = new double[maxNumrows_];
    //rowsense_ = new char[maxNumrows_];
    //rhs_      = new double[maxNumrows_];
    //rowrange_ = new double[maxNumrows_];
    dual = new double[maxNumrows_];
    lhs_      = new double[maxNumrows_];
    actv =  new int [maxNumrows_];
    std::cout<<"OsiVolSolverInterface::rowRimAllocator_"<<std::endl;*/

}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::colRimAllocator_()
{
    /*colub  = new double[maxNumcols_];
    collb  = new double[maxNumcols_];
    //continuous_ = new bool[maxNumcols_];
    //objcoeffs_ = new double[maxNumcols_];
    solution    = new double[maxNumcols_];
    rc_        = new double[maxNumcols_];
    if(!yhit)yhit =  new double [narcs];
    if(!VItopo)VItopo = new double [narcs];
    if(!addrc)addrc = new double [narcs];
    if(!arc_map)arc_map =  new int [narcs];
    std::cout<<"OsiVolSolverInterface::colRimAllocator_"<<std::endl;*/


}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::gutsOfDestructor_()
{
	/*std::cout<<"OsiVolSolverInterface::gutsOfDestructor_()"<<std::endl;
    //rowMatrix_.clear();
    //colMatrix_.clear();
    //rowMatrixCurrent_ = true;
    //colMatrixCurrent_ = true;
    //delete[] continuous_; continuous_ = 0;
    //delete[] rowsense_;	rowsense_ = 0;
    //delete[] rhs_;	rhs_ = 0;
    //delete[] rowrange_;	rowrange_ = 0;
    //delete[] objcoeffs_;	objcoeffs_ = 0;
    if(rowub){delete[] rowub;    rowub = 0;}
    if(rowlb){delete[] rowlb;    rowlb = 0;}
    if(colub){delete[] colub; colub = 0;}
    if(collb){delete[] collb;    collb = 0;}
    if(solution){delete[] solution;	        solution = 0;}
    if(dual){delete[] dual;	        dual = 0;}
    if(rc_){delete[] rc_;     rc_ = 0;}
    if(lhs_){delete[] lhs_;    lhs_ = 0;}
    if(actv){delete[] actv; actv=0;}
    numrows_ = maxNumrows_ = 0;
    numcols_ = maxNumcols_ = 0;*/
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Constructors, destructors clone and assignment
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

OsiVolSolverInterface::OsiVolSolverInterface () :
colub(0),
collb(0),
rowub(0),
rowlb(0),
solution(0),
dual(0),
rc_(0),
lhs_(0),
yhit(0),
arc_map(0),
VItopo(0),
addrc(0),
actv(0),
HotStart_(0),
volprob_(0)
{

    maxNumrows_ = maxNumcols_ = 0;
    mode =1;
    num_purgbl=0;
    min_lower_bound=0;
    has_checked= recheck_collct=in_strong_branch = HotStartSet = false;
}



//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::initialize(OsiVolAuxInfo & osidata){
	data = osidata.data;
	nnodes =  data->nnodes;
    ndemands = data->ndemands;
    narcs = data->narcs;
    
    cover_manager = osidata.cover_manager;
    ss_manager = osidata.ss_manager;
    localc_manager = osidata.localc_manager;
    globalc_manager =  osidata.globalc_manager;
    lpchecker = osidata.lpchecker;

    osidata.init_solver(nnodes*ndemands);
    volprob_ = &osidata.volprob;

    maxNumVI =osidata.maxNumVI;
    intvlVI = osidata.intvlVI;
    
	int maxvi = osidata.maxPos;
	
	
    rowub = new double[maxvi];
    rowlb = new double[maxvi];
    dual = new double[maxvi];
    lhs_      = new double[maxvi];
    actv =  new int [maxvi];
    
   	int ncol = narcs+ndemands*narcs;
   	colub  = new double[ncol];
    collb  = new double[ncol];
    solution    = new double[ncol];
    rc_        = new double[ncol];
    yhit =  new double [narcs];
    VItopo = new double [narcs];
    arc_map =  new int [narcs];
    
    CoinFillN(rowub, maxvi, 0.0);
    CoinFillN(rowlb, maxvi, 0.0);
    CoinFillN(dual, maxvi, 0.0);
    CoinFillN(lhs_, maxvi, 0.0);
    CoinFillN(actv, maxvi, -1);
    
	CoinFillN(colub, ncol, OsiVolInfinity);
	CoinFillN(collb, ncol, 0.0);
    CoinFillN(rc_, ncol, 0.0);
    
    osidata.rowub = rowub;
    osidata.rowlb = rowlb;
    osidata.colub = colub;
    osidata.collb = collb;
    osidata.dual = dual;
    osidata.lhs_ = lhs_;
    osidata.actv = actv;
	osidata.rc_ = rc_;
	osidata.VItopo = VItopo;
	osidata.arc_map = arc_map;
	osidata.yhit = yhit;
	osidata.solution = solution;
    //if(!addrc)addrc = new double [data->narcs];
}

//---------------------------------------------------------------------------

OsiSolverInterface *
OsiVolSolverInterface::clone(bool copyData) const {
    OsiVolSolverInterface * c;
    if(copyData)
        c= new OsiVolSolverInterface(*this);
    else
        c =new OsiVolSolverInterface();
    
    return c;
}

//-----------------------------------------------------------------------

OsiVolSolverInterface::OsiVolSolverInterface(const OsiVolSolverInterface& x) :
OsiSolverInterface(x),
volprob_()
{
	has_checked = x.has_checked;
	HotStart_ = x.HotStart_;
    maxNumrows_ = x.maxNumrows_;
    maxNumcols_ = x.maxNumcols_;
    mode = x.mode;
    recheck_collct=HotStartSet = false;
    num_purgbl=x.num_purgbl;
    min_lower_bound=x.min_lower_bound;
    in_strong_branch = x.in_strong_branch;
    colub = x.colub;
	collb = x.collb;
	rowub = x.rowub;
	rowlb = x.rowlb;
	solution = x.solution;
	dual = x.dual;
	rc_ = x.rc_;
	lhs_ = x.lhs_;
	yhit = x.yhit;
	arc_map = x.arc_map;
	VItopo = x.VItopo;
	addrc = x.addrc;
	actv = x.actv;
	
	data = x.data;   
	cover_manager = x.cover_manager;
    ss_manager = x.ss_manager;
    localc_manager = x.localc_manager;
    globalc_manager = x.globalc_manager;
    volprob_ = x.volprob_;
    lpchecker = x.lpchecker;
	
}

//-----------------------------------------------------------------------

OsiVolSolverInterface&
OsiVolSolverInterface::operator=(const OsiVolSolverInterface& x){
    if (&x == this)
        return *this;
    
    colub = x.colub;
	collb = x.collb;
	rowub = x.rowub;
	rowlb = x.rowlb;
	solution = x.solution;
	dual = x.dual;
	rc_ = x.rc_;
	lhs_ = x.lhs_;
	yhit = x.yhit;
	arc_map = x.arc_map;
	VItopo = x.VItopo;
	addrc = x.addrc;
	actv = x.actv;
	HotStart_ = x.HotStart_;
	has_checked = x.has_checked;
	
    mode = x.mode;
    data = x.data;   
    
	cover_manager = x.cover_manager;
    ss_manager = x.ss_manager;
    localc_manager = x.localc_manager;
    globalc_manager = x.globalc_manager;
    volprob_ = x.volprob_;
    lpchecker = x.lpchecker;
    return *this;
}

//-----------------------------------------------------------------------

OsiVolSolverInterface::~OsiVolSolverInterface (){
	unmarkHotStart();
}







