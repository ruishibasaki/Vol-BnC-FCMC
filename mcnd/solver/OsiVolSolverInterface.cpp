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


//---------------------------------------------------------------------------
// WarmStart related methods */

CoinWarmStart *
OsiVolSolverInterface::getEmptyWarmStart () const{
    return (dynamic_cast<CoinWarmStart *>(new CoinWarmStartDual())) ;
}

//---------------------------------------------------------------------------

CoinWarmStart* 
OsiVolSolverInterface::getWarmStart() const{
    std::cout<<"getWarmStart() "<<getNumRows()<<std::endl;
    return new WarmStartDual(getNumRows(), dual, actv);
    
}

//---------------------------------------------------------------------------

bool 
OsiVolSolverInterface::setWarmStart(const CoinWarmStart* warmstart){
	
    HotStart_ = dynamic_cast<const WarmStartDual*>(warmstart);
	
    if (! HotStart_){
        const CoinWarmStartDual* hs = dynamic_cast<const CoinWarmStartDual*>(warmstart);
        if(!hs){
        	HotStartSet = false;
        	return false;
        }
    	HotStart_ = new WarmStartDual(hs);
    }
    std::cout<<"setWarmStart "<<HotStart_->size()<<" "<<getNumRows()<<std::endl;

    HotStartSet = true;
    return true;
};

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::markHotStart(){
    std::cout<<"markHotStart() "<<getNumRows()<<std::endl;
    if(HotStart_) delete HotStart_;
    HotStart_ = new WarmStartDual(getNumRows(), dual, actv); 
    //std::cout<<"markHotStart() now "<<std::endl;
    HotStartSet = true;
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::unmarkHotStart(){
    std::cout<<"unmarkHotStart() "<<std::endl;
    if(HotStart_) delete HotStart_;
    HotStart_ = 0;
    HotStartSet = false;
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::solveFromHotStart(){
    //std::cout<<"solveFromHotStart() maxiter: "<<volprob_.parm.maxsgriters<<" dsize: "<<numrows_<<std::endl;
    //CoinDisjointCopyN(HotStart_, numrows_, dual);
    //OsiVolAuxInfo * auxinfo  = static_cast<OsiVolAuxInfo*>(OsiSolverInterface::getApplicationData());
    resolve();
    
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::translate_sol(){
    
    CoinFillN(solution, getNumCols(), 0.0);

    int arc, idx;
    for(int a=szopnd; a<sznz;++a){
        arc = nz_arcs[a];
        
        solution[arc] = volprob_.psol[a];
        //std::cout<<arc<<" "<<volprob_.psol[a]<<" "<<h1[arc]<<std::endl;// /double(volprob_.iter())<<std::endl;
        for(int k=ndemands; k--; )
            solution[narcs+k*narcs+arc] = volprob_.psol[szunfxd + k*sznz + a];
    }
    for(int a=szopnd; a--;){
        arc = nz_arcs[a];
        solution[arc] = 1;
        for(int k=ndemands; k--; )
            solution[narcs+k*narcs+arc] = volprob_.psol[szunfxd + k*sznz + a];
    }
    
    int cidx = ndemands*nnodes + cover_manager->covers.sizeOfCollection;
    for(int i=cidx; i-- ;){
        idx = actv[i];
        if(idx>=0){
            //if(i>=ndemands*nnodes) std::cout<<i<<": "<<volprob_.dsol[idx]<<std::endl;
            dual[i] = volprob_.dsol[idx];
            lhs_[i] = volprob_.viol[idx];
        }else{
            dual[i] =0;
            lhs_[i] =0;
        }
    }
}


//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::setRowSetBounds(const int* indexFirst,
                                       const int* indexLast,
                                       const double* boundList){
    std::cout<<"rowsetbounds "<<std::endl;
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
    std::cout<<"colsetbounds "<<std::endl;
    while (indexFirst != indexLast) {
        collb[*indexFirst] =  boundList[0];
        colub[*indexFirst] = boundList[1];
        ++indexFirst;
        boundList += 2;
        //std::cout<<"col: "<<*indexFirst<<std::endl;
    }
    nz_arcs.clear();
    szopnd=0;
    szunfxd=0;
    for(int arc=data->narcs;arc--;){
        if(collb[arc] == 1){
            nz_arcs.push_back(arc);
            arc_map[arc] = -2;
            VItopo[arc] =  yhit[arc] = 1.0;
            ++szopnd;
            std::cout<<"opened arc: "<<arc<<std::endl;
            continue;
        }else if(colub[arc] == 1){
            arc_map[arc] = szunfxd++;
            nz_arcs.push_front(arc);
            std::cout<<"unfix: "<<arc<<" idx: "<<arc_map[arc]<<std::endl;
        }else{
            std::cout<<"closed arc: "<<arc<<std::endl;
            arc_map[arc] = -1;
        }
        VItopo[arc] = yhit[arc] = 0.0;
    }
    sznz = nz_arcs.size();
    cover_manager->set_arc_map(arc_map);
    
}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::map_duals(){

    fsize= csize = 0;
    bool flag=false;
    const Arc* item; const Demand* itemd;
    CoinFillN(actv, maxNumcols_, -1);
    for(int k=0; k<ndemands; ++k ){
        for(int i=0; i<nnodes; ++i ){
            flag = false;
            
            for(int a=sznz; a--;){
                int arc = nz_arcs[a];
                item = &data->arcs[arc];
                if((i+1) == item->i){
                    flag = true;
                    break;
                }else if((i+1) == item->j){
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
            }else mode =0;
        }
    }
    if(cover_manager->reset_and_map_collection(fsize, VItopo, dual, actv, csize))
    	HotStartSet =false;
    
    //std::cout<<"ok"<<std::endl;
}

//---------------------------------------------------------------------------

void 
OsiVolSolverInterface::translate_hotstart(){
    int fidx = ndemands*nnodes;
    int idx, sz;
	for(int i=fidx; i--;){
		idx = actv[i];
		if(idx>=0){
			volprob_.dsol[idx] = HotStart_->at(i);
		}
	}
	sz = cover_manager->num_actv;
	Cover* vi = cover_manager->covers.begin;
	for(int i=sz; i--;){
		idx = actv[vi->id_vi];
		if(idx>=0){
			volprob_.dsol[idx] = HotStart_->at(vi->id_vi);
			std::cout<<"vi: "<<vi->serial_nmbr<<" id "<<vi->id_vi<<" : "<<volprob_.dsol[idx]<<std::endl;
		}
		vi = vi->next;
	}
}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::set_start(){
    std::cout<<"OsiVolSolverInterface::set_start "<<std::endl;
    int idx, sz;
    int fidx = ndemands*nnodes;
    volprob_.dsol = 0.0;
    volprob_.dual_lb = 0.0; volprob_.dual_ub = 0.0;
    
    for(int i=fidx; i--;){
        //std::cout<<old_index[i]<<" "<<Iu[i]<<" value: ";
        //std::cout<<old_dual[old_index[i]]<<std::endl;
        idx = actv[i];
        if(idx>=0){
            volprob_.dual_lb[idx] = -1.0e31;
            volprob_.dual_ub[idx] = 1.0e31;
        }
    }
    sz = cover_manager->num_actv;
    Cover* vi = cover_manager->covers.begin;
    for(int i=sz; i--;){
        idx = actv[vi->id_vi];
        //std::cout<<"vi "<<vi->id_vi<<": "<<idx<<std::endl;
        if(idx>=0){
            volprob_.dual_ub[idx] = 1.0e31;
        	//std::cout<<"vi: "<<vi->serial_nmbr<<" id "<<vi->id_vi<<" : "<<volprob_.dsol[idx]<<std::endl;
        }
        vi = vi->next;
    }
    if(HotStartSet){
		translate_hotstart();
    }
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Solve triggers
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

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
	if(mode==-1) return;
	std::cout<<"mode: "<<mode<<std::endl;
    int i;
    map_duals();
    volprob_.active_size = fsize + csize;
    
    if(maxNumrows_ < fsize + csize+1000)
        maxNumrows_ = fsize + csize+1000;
    
    volprob_.dsize = maxNumrows_;
    volprob_.psize = szunfxd + data->ndemands*sznz;
    volprob_.dsol.allocate(maxNumrows_);
    volprob_.dual_lb.allocate(maxNumrows_);
    volprob_.dual_ub.allocate(maxNumrows_);
    set_start();

    std::cout<<"re solve "<<szunfxd<<" dsize: "<<numrows_<<" maxdsz: "<<maxNumrows_<<std::endl;
    
    
    // Set the dual starting point
    retval=volprob_.solve(*this, true);
    numrows_ +=  cover_manager->gend;
    std::cout<<std::setprecision(10)<<"solve: "<<volprob_.value<<" retval: "<<retval<<" iters: "<<volprob_.iter()<<std::endl;
    translate_sol();
    if(!isProvenPrimalInfeasible()){
        lagrangeanCost_ = volprob_.value;
    }else lagrangeanCost_=0;
    
    // extract the solutio
    // the lower bound on the objective value
}

//---------------------------------------------------------------------------

void OsiVolSolverInterface::direct_solve(const std::deque<int>& topo, const CoinWarmStart* warmstart){
    nz_arcs = topo;
    sznz = szopnd = nz_arcs.size();
    szunfxd = 0;
    markHotStart();
    setWarmStart(warmstart);
    resolve();
    
    CoinDisjointCopyN(HotStart_->dual(), getNumRows(), dual);
    unmarkHotStart();
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// VOLUME USER HOOKS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


int
OsiVolSolverInterface::addVI(int iter,double lcost, const VOL_dvector& xstar,
          const VOL_dvector& x, VOL_dvector& dual_lb, VOL_dvector& dual_ub,
          VOL_dvector& v, VOL_dvector& h, int & actvSSz){
    if(mode == 0) return 0;
    if(maxNumVI<=(cover_manager->covers.sizeOfCollection)) return 0;
    if(iter%intvlVI==0 ){
        translate_primal(xstar);
        int num_new_sets=ss_manager->cutset_generation_main( yhit, VItopo, false);
        VIub=-1e31;
    }
    
    if(lcost >= VItt ){
        VItt = lcost;
        if(iter<100) return 0;
        int num_covers = cover_manager->cover_generation_main(xstar.v, x.v, &ss_manager->sets, actvSSz);
        cover_manager->add_cover_vi(num_covers, actv, actvSSz, h.v, dual_lb.v,  dual_ub.v );
        
        if(num_covers>0){
            std::cout<<std::setprecision(10)<<"iter: "<<iter<<" added cuts "<<num_covers<<" L: "<<VItt<<std::endl;
        }
        
        
    }
    return 0;
}


//-----------------------------------------------------------------------

int
OsiVolSolverInterface::removeVI( int & actvSSz,VOL_dvector& pstarv, VOL_dvector& dstaru,  VOL_dvector& dualu){
    //cover_manager->covers.removeCover(lim_to_remv, actv, actvSSz, pstarv.v, dstaru.v, dualu.v);
    return 0;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::compute_rc(const VOL_dvector& dual, VOL_dvector& rc, int actvSSz){
    const Arc* item;
    int arc;

    for(int a=szunfxd; a--;){
        arc = nz_arcs[a];
        item = &data->arcs[arc];
        rc[a] = item->f;
        addrc[a] =0;
        for(int k=ndemands; k--; ){
            rc[szunfxd + k*sznz + a] = item->c[k] - dual[actv[k*nnodes + item->j-1]] + dual[actv[k*nnodes + item->i-1]];
        }
    }
    for(int a=szunfxd; a<sznz;++a){
        arc = nz_arcs[a];
        item = &data->arcs[arc];
        for(int k=ndemands; k-- ;){
            rc[szunfxd + k*sznz + a] = item->c[k] - dual[actv[k*nnodes + item->j-1]] + dual[actv[k*nnodes + item->i-1]];
        }
    }

    B0=0;
    for(int k=0; k<data->ndemands; ++k){
        B0 += data->d_k[k].quantity * ( dual[actv[k*nnodes + data->d_k[k].D-1]] - dual[actv[k*nnodes + data->d_k[k].O-1]]);
        
    }
    double b = B0;
    cover_manager->compute_cover_rc( dual.v,  actv,  actvSSz, addrc, b);
    return 0;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::compute_sg(const VOL_dvector& x, int actvSSz, VOL_dvector& v){
    
    v =0;
    int arc;
    const Arc* item; const Demand* itemd;
    for(int k=ndemands; k--; ){
        itemd = &data->d_k[k];
        v[actv[k*nnodes + itemd->O-1]] -= itemd->quantity;
        v[actv[k*nnodes + itemd->D-1]] += itemd->quantity;
        for(int a=0; a<sznz; ++a){
            arc = nz_arcs[a];
            v[actv[k*nnodes + data->arcs[arc].i-1]] += x[szunfxd+k*sznz+ a];
            v[actv[k*nnodes + data->arcs[arc].j-1]] -= x[szunfxd+ k*sznz+ a];
        }
    }
    cover_manager->compute_cover_sg( x.v, actv,  actvSSz, v.v);
    return 0;
}


//-----------------------------------------------------------------------

int
OsiVolSolverInterface::solve_subproblem(const VOL_dvector& xstar,
                                        const VOL_dvector& dual,  VOL_dvector& rc,
                                        double& lcost, VOL_dvector& x,
                                        double& pcost){
    int arc;
    double cost_a, ratio;
    pcost= 0;
    lcost =0;
    for(int a=szunfxd; a--; ){
        arc = nz_arcs[a];
        rc[a] += knapsack(a, rc.v, x.v);
        addrc[a]+= rc[a];
        if(addrc[a]<1e-10 && addrc[a]>-1e-10) addrc[a] = 0;
        if( addrc[a] < 0.0){
            x[a] =1.0;
        }else x[a] = 0.0;
    }
    for(int a=szunfxd; a<sznz;++a){
        arc = nz_arcs[a];
        cost_a = knapsack(a, rc.v, x.v);
        lcost += cost_a + data->arcs[arc].f;
        pcost += data->arcs[arc].f;
        for(int k=ndemands; k--; )
            pcost += data->arcs[arc].c[k] * x[szunfxd + k*sznz + a];
    }
    
    return 0;
}

//---------------------------------------------------------------------------

int
OsiVolSolverInterface::resolve_subproblem(const VOL_dvector& dual, VOL_dvector& rc,
                       double& lcost,
                       VOL_dvector& x,double& pcost){
    
    int arc;
    lcost += B0;
    for(int a=szunfxd; a--; ){
        arc = nz_arcs[a];
        if( rc[a]<0 || (rc[a]==0 && x[a]==1.0)){
            lcost += rc[a];
            x[a]=1.0;
            pcost += data->arcs[arc].f;
            for(int k=0; k<ndemands; ++k)
                pcost += data->arcs[arc].c[k] * x[szunfxd + k*sznz + a];
        }else{
            if(x[a]==1.0) std::cout<<"rarc: "<<a<<" rc: "<<rc[a]<<" y: "<<x[a]<<std::endl;
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
OsiVolSolverInterface::additional_settings(int iter, double& lcost, VOL_dvector& dual, VOL_dvector& rc, VOL_dvector& h,
                          VOL_dvector& x, const VOL_dvector& xhist, const VOL_dvector& dstar, int actvSSz){
   
    if(cover_manager->covers.sizeOfCollection==0 || iter==0){
        cover_manager->compute_cover_rc( dual.v, actv,  actvSSz, rc.v,   B0);
        return 0;
    }
    double ret;
    //std::cout<<"rc[id_arc]: "<<rc[arc_map[59]]<<std::endl;
    //std::cout<<"opaa one: "<<33<<" rc: "<<addrc[33]<<" rc+: "<<rc[33]<<" x: "<<x[33]<<std::endl;
    cover_manager->recompute_mult_neg( dual.v, addrc, rc.v, x.v, actv, actvSSz);
    cover_manager->recompute_mult_pos( dual.v, h.v, addrc, dstar.v , x.v, actv);
    cover_manager->compute_cover_rc( dual.v, actv,  actvSSz, rc.v,   B0);
    for(int a=szunfxd; a--; ){
        if(rc[a]<1e-10 && rc[a]>-1e-10) rc[a] = 0;
        if(rc[a]==0 ){
            ret = arc_dg_imp(a, x.v, h.v , actvSSz);
            if(ret>0 && xhist[a]>=0.0){
                //std::cout<<"yh: "<<xhist[a]<<std::endl;
                //std::cout<<"opaa zero to one: "<<a<<" rc: "<<addrc[a]<<" rc+: "<<rc[a]<<" y: "<<x[a]<<std::endl;
                x[a] = 1.0 ;
            }else if(ret<0 && xhist[a]<=0.0)  x[a] = 0.0 ;
        }
    }
    return 0;
}

//---------------------------------------------------------------------------
//  solving methods
//---------------------------------------------------------------------------

double
OsiVolSolverInterface::knapsack(int a, const double * rc, double* x){
    double kpsack =0;
    double fillUp =0;
    int arc = nz_arcs[a];
    std::list<HeapCell> heap;
    //get reduced cost for each commodity in arc e
    for(int k=ndemands; k--; ){
        if(rc[szunfxd + k*sznz + a]<0.0){
            heap.push_back(HeapCell(k,rc[szunfxd + k*sznz + a]));
        }else x[szunfxd + k*sznz + a]=0.0;
    }
    
    heap.sort(comp());
    //std::stable_sort(heap.begin(), heap.end(), comp());
    
    double capa;
    int comm;
    while(heap.size()>0){
        capa = data->arcs[arc].capa;
        comm = heap.back().k;
        if(fillUp < capa){
            x[szunfxd + comm*sznz + a] = std::min((capa - fillUp),  data->d_k[comm].quantity);
            fillUp += x[szunfxd + comm*sznz + a];
            kpsack += heap.back().rc_ * x[szunfxd + comm*sznz + a];
            heap.pop_back();
            //std::cout<<"in"<<std::endl;
        }else{
            x[szunfxd + comm*sznz + a] = 0.0;
            heap.pop_back();
        }
    }
    return kpsack;
}

//---------------------------------------------------------------------------

int
OsiVolSolverInterface::heuristics(const VOL_problem& p,
                                  const VOL_dvector& x, double& heur_val){ return 0;}

//---------------------------------------------------------------------------
//  auxiliary methods
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
    //if(dg>0) std::cout<<"arc: "<<arc<<" dg: "<<dg<<std::endl;
    return dg;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::mark_topo( VOL_dvector& x, double lcost){
    
    //std::cout<<"Lt: "<<lcost<<std::endl;
    if(lcost<=VIub) return 0;
    int arc;
    VIub =lcost;
    for(int a=szopnd; a<sznz;++a){
        arc = nz_arcs[a];
        VItopo[arc] = x[a];
    }
    return 0;
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::translate_primal( const VOL_dvector& xhist){
    int arc;
    for(int a=szopnd; a<sznz;++a){
        arc = nz_arcs[a];
        yhit[arc] = xhist[a];
    }
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
    
    std::cout<<"OsiVolSolverInterface::loadProblem nrows: "<<numrows<<std::endl;
    OsiVolAuxInfo * auxinfo  = static_cast<OsiVolAuxInfo*>(OsiSolverInterface::getApplicationData());
    data = auxinfo->data;
    cover_manager = auxinfo->cover_manager;
    ss_manager = auxinfo->ss_manager;
    
    nnodes =  data->nnodes;
    ndemands = data->ndemands;
    narcs = data->narcs;
    
    VItt = VIub=-1e31;
    lim_to_remv = cover_manager->lim_to_remv;
    maxNumVI =auxinfo->maxNumVI;
    intvlVI = auxinfo->intvlVI;
    
    
    gutsOfDestructor_();
    maxNumcols_ = numcols_ = numcols;
    maxNumrows_ = numrows+maxNumVI;
    numrows_ = numrows+cover_manager->covers.sizeOfCollection;
    //cover_manager->reset_collection();
    
    if (maxNumrows_ > 0) {
        rowRimAllocator_();
        CoinFillN(rowub, maxNumrows_, 0.0);
        CoinFillN(rowlb, maxNumrows_, 0.0);
        // Set the initial dual solution
        CoinFillN(dual, maxNumrows_, 0.0);
        CoinFillN(lhs_, maxNumrows_, 0.0);
        
        if (rowub) {  CoinDisjointCopyN(rowub_, numrows, rowub); }
        if (rowlb) { CoinDisjointCopyN(rowlb_, numrows, rowlb); }
        
    }
    
    if (maxNumcols_ > 0) {
        colRimAllocator_();
        if (colub) {
            CoinDisjointCopyN(colub_, numcols, colub);
        } else {
            CoinFillN(colub, numcols, OsiVolInfinity);
        }
        if (collb) {
            CoinDisjointCopyN(collb_, numcols, collb);
        } else {
            CoinFillN(collb, numcols, 0.0);
        }
        // Set the initial rc solution
        CoinFillN(rc_, numcols, 0.0);
        int c;
        for ( c=0; c<numcols; c++ ) {
            if ( fabs(collb[c]) < fabs(colub[c]) ) {
                solution[c] = collb[c];
            }
            else {
                solution[c] = colub[c];
            }
        }
    }
    
    
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::rowRimAllocator_()
{
    rowub = new double[maxNumrows_];
    rowlb = new double[maxNumrows_];
    //rowsense_ = new char[maxNumrows_];
    //rhs_      = new double[maxNumrows_];
    //rowrange_ = new double[maxNumrows_];
    dual = new double[maxNumrows_];
    lhs_      = new double[maxNumrows_];
    actv =  new int [maxNumrows_];
    //std::cout<<"OsiVolSolverInterface::rowRimAllocator_"<<std::endl;

}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::colRimAllocator_()
{
    colub  = new double[maxNumcols_];
    collb  = new double[maxNumcols_];
    //continuous_ = new bool[maxNumcols_];
    //objcoeffs_ = new double[maxNumcols_];
    solution    = new double[maxNumcols_];
    rc_        = new double[maxNumcols_];
    if(!yhit)yhit =  new double [narcs];
    if(!VItopo)VItopo = new double [narcs];
    if(!addrc)addrc = new double [narcs];
    if(!arc_map)arc_map =  new int [narcs];
    //std::cout<<"OsiVolSolverInterface::colRimAllocator_"<<std::endl;


}

//---------------------------------------------------------------------------

void
OsiVolSolverInterface::gutsOfDestructor_()
{
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
    lagrangeanCost_ = 0.0;
    numrows_ = maxNumrows_ = 0;
    numcols_ = maxNumcols_ = 0;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Constructors, destructors clone and assignment
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

OsiVolSolverInterface::OsiVolSolverInterface () :
//rowsense_(0),
//rhs_(0),
//rowrange_(0),
//continuous_(0),
//objcoeffs_(0),
//objsense_(1.0),
volprob_(),
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
HotStart_(0)
{
    lagrangeanCost_ = maxNumrows_ = maxNumcols_ = 0;
    mode =1;
    HotStartSet = false;
}

//---------------------------------------------------------------------------

OsiVolSolverInterface::OsiVolSolverInterface(const char * volparfile):
volprob_(volparfile),
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
HotStart_(0)
{
    lagrangeanCost_ = maxNumrows_ = maxNumcols_ = 0;
    mode =1;
    HotStartSet = false;
}


//---------------------------------------------------------------------------

OsiSolverInterface *
OsiVolSolverInterface::clone(bool copyData) const {
    OsiVolSolverInterface * c;
    if(copyData)
        c= new OsiVolSolverInterface(*this);
    else
        c =new OsiVolSolverInterface();
    
    //c->volprob_.parm.temp_dualfile =volprob_.parm.temp_dualfile ;
    c->volprob_.parm.ubinit = volprob_.parm.ubinit;
    c->volprob_.parm.printflag = volprob_.parm.printflag;
    c->volprob_.parm.printinvl = volprob_.parm.printinvl ;
    c->volprob_.parm.maxsgriters = volprob_.parm.maxsgriters;
    c->volprob_.parm.maxtime = volprob_.parm.maxtime;
    //c->volprob_.parm.heurinvl = volprob_.parm.heurinvl;
    c->volprob_.parm.greentestinvl =volprob_.parm.greentestinvl ;
    c->volprob_.parm.yellowtestinvl = volprob_.parm.yellowtestinvl;
    c->volprob_.parm.redtestinvl = volprob_.parm.redtestinvl;
    c->volprob_.parm.lambdainit = volprob_.parm.lambdainit ;
    c->volprob_.parm.alphainit = volprob_.parm.alphainit;
    c->volprob_.parm.alphamin = volprob_.parm.alphamin ;
    c->volprob_.parm.alphafactor =volprob_.parm.alphafactor ;
    c->volprob_.parm.alphaint = volprob_.parm.alphaint;
    c->volprob_.parm.primal_abs_precision = volprob_.parm.primal_abs_precision;
    c->volprob_.parm.gap_abs_precision = volprob_.parm.gap_abs_precision;
    c->volprob_.parm.gap_rel_precision =volprob_.parm.gap_rel_precision ;
    //c->volprob_.parm.ubfeas = volprob_.parm.ubfeas;
    //c->volprob_.parm.infeas_trigger = volprob_.parm.infeas_trigger;
    //c->volprob_.parm.granularity = volprob_.parm.granularity ;
    return c;
}

//-----------------------------------------------------------------------

OsiVolSolverInterface::OsiVolSolverInterface(const OsiVolSolverInterface& x) :
OsiSolverInterface(x),
volprob_(),
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
HotStart_(0)
{
    lagrangeanCost_ = maxNumrows_ = maxNumcols_ = 0;
    mode =1;
    HotStartSet = false;
    operator=(x);
}

//-----------------------------------------------------------------------

OsiVolSolverInterface&
OsiVolSolverInterface::operator=(const OsiVolSolverInterface& rhs){
    if (&rhs == this)
        return *this;
    
    OsiSolverInterface::operator=(rhs);
    gutsOfDestructor_();
    
    if (rhs.maxNumrows_) {
        maxNumrows_ = rhs.maxNumrows_;
        rowRimAllocator_();
        const int rownum = rhs.numrows_;
        CoinDisjointCopyN(rhs.rowub, rownum, rowub);
        CoinDisjointCopyN(rhs.rowlb, rownum, rowlb);
        //CoinDisjointCopyN(rhs.rowsense_, rownum, rowsense_);
        //CoinDisjointCopyN(rhs.rhs_, rownum, rhs_);
        //CoinDisjointCopyN(rhs.rowrange_, rownum, rowrange_);
        CoinDisjointCopyN(rhs.dual, rownum, dual);
        CoinDisjointCopyN(rhs.lhs_, rownum, lhs_);
    }
    if (rhs.maxNumcols_) {
        maxNumcols_ = rhs.maxNumcols_;
        colRimAllocator_();
        const int colnum = numcols_;
        CoinDisjointCopyN(rhs.colub, colnum, colub);
        CoinDisjointCopyN(rhs.collb, colnum, collb);
        // CoinDisjointCopyN(rhs.continuous_, colnum, continuous_);
        //CoinDisjointCopyN(rhs.objcoeffs_, colnum, objcoeffs_);
        CoinDisjointCopyN(rhs.solution, colnum, solution);
        CoinDisjointCopyN(rhs.rc_, colnum, rc_);
    }
    //volprob_.parm.granularity = 0.0;
    return *this;
}

//-----------------------------------------------------------------------

OsiVolSolverInterface::~OsiVolSolverInterface (){
    if(yhit) delete [] yhit;
    if(arc_map) delete [] arc_map;
    if(VItopo) delete [] VItopo;
    if(addrc) delete [] addrc;
    gutsOfDestructor_();
}









