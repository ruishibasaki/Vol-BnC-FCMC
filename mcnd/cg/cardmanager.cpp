//
//  covermanager.cpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#include "covermanager.hpp"



//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//   min cardinality methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
CoverManager::mincard_generation_main(const double * ystar, const double * y, const CutSetCollection * sets, int curr_id, int max){
    CutSet * cutset;
    int added=0;
    int sz = sets->sizeOfCollection;
     //sets->print();
    cutset = sets->begin;
    for(int i=0;i<sz;++i){

        added += mincard_generation(cutset->ss_size, cutset->SS_arcs,cutset->ss_ksize, cutset->SS_comm, cutset->uss, cutset->dss, ystar, y, curr_id+added);
        if(added+curr_id>=max) break;
        added += mincard_generation(cutset->s_ssize, cutset->S_Sarcs, cutset->s_sksize, cutset->S_Scomm, cutset->us_s, cutset->ds_s, ystar, y, curr_id+added);
        if(added+curr_id>=max) break;
        //std::cout<<"i: "<<i<<std::endl;
        cutset = cutset->next;
    }
    gend+=added;
   // ttgend_card+=added; //stat can be removed
    return added;
}

//-------------------------------------------------------------------------------------------

int
CoverManager::mincard_generation(int ss_size, const int * SS_arcs, int ss_ksize, const int * SS_comm, double uss, double dss,
								 const double * ystar, const double * y, int id_vi){
    
    std::deque<Pair2> ss_;
    double  u1, delta, rhs;
    int card;
    int added=0;
    double viol=0;
    
    //std::cout<<"dss:  "<<dss<<" uss: "<<uss<<std::endl;
    u1 =  card_cutset_preprocess(ss_size, SS_arcs,  ss_ksize, SS_comm,  ss_, y, ystar, dss);
    delta = dss - u1;
    if(delta<=0){ss_.clear(); return 0;}
    
    //std::cout<<"delta: "<<delta<<std::endl;
    card = make_cardcs( ss_,  delta );
     
    Cover* cover = covers.createNewCover(ss_, id_vi, ttgend);
    //std::cout<<"done add cove"<<std::endl;
    cover->rhs = card;
    cover->owner = SS_arcs;
    cover->maxsize = ss_size;
    cover->hs = checkViol(cover, ystar);
    if(cover->hs>0){
        added = covers.addCover(cover, ystar);
        if(added){
        	++ttgend;
        	++num_actv;
        	//cover->print();
        }
         
    }else delete cover;
    
    ss_.clear();
    return added;
}


//-------------------------------------------------------------------------------------------

int
CoverManager::make_cardcs( std::deque<Pair2> & ss_, double delta ){
    int sz =(int)ss_.size();
    int card =0;
    int arc;
    double bl;
    double ttu=0;
    std::stable_sort(ss_.begin(), ss_.end(), compPair2());
    //std::cout<<"new: "<<std::endl;
    for(int n=0;n<sz;n++){
        arc = ss_[n].fst;
        bl = ss_[n].snd;
        ttu+= bl;
        if(ttu>=delta){
            break;
        }
        card = n+1;
    }
    
    return card+1;
}

//----------------------------------------------------------------------------------

double
CoverManager::card_cutset_preprocess(int sz, const int * ss_, int szk, const int * ss_k,  std::deque<Pair2>& ss_deque,
                                const double *y, const double *ystar, double dss){
    int arc, id_arc, k;
    double U1=0;
    double bl=0;
    double maxf=0;
    bool skip=true;
    //std::cout<<"cut: "<<dss<<std::endl;
    int narcs = data->narcs;
    for(int i=sz;i--;){
        arc = ss_[i];
        id_arc = arc_map[arc];
        maxf=0;
        for(int ik=szk; ik--;){
        	k = ss_k[ik];
        	maxf+=colub[narcs+k*narcs+arc];
        }
        bl = fmin(maxf,data->arcs[arc].capa);
        //std::cout<<"arc "<<arc <<" bl: "<<bl<<" "<<fmin(dss,data->arcs[arc].capa)<<std::endl;
        if(bl < fmin(dss,data->arcs[arc].capa)) skip=false;
    
        
        if(id_arc>=0){
            if(y[id_arc]>0.5 || ystar[id_arc]>0.9){
                U1+=bl;
            }else{
                ss_deque.push_back(Pair2(arc, bl));
            }
        }else if(id_arc == -2){ U1+=bl;}
    }
    if(skip) return 1e30;
    else return U1;
}

