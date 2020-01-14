//
//  covermanager.cpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#include "covermanager.hpp"



//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  initializing methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void
CoverManager::initialize(const Data * d, int lim) {
    data=d;
    coverlift.set_data(d);
    covers.initialize(data->narcs);
    lim_to_remv = lim;
    num_actv = gend = ttgend= 0;
}

//-------------------------------------------------------------------------------------------

int
CoverManager::reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize){
    int cont;
    int idxf = data->nnodes*data->ndemands;
    int sz = covers.sizeOfCollection;
    Cover* vi = covers.begin;
    num_actv =0; cont=0;
    for(int i=0;i<sz;++i){
        vi->n_zerom =0;
        vi->n_nviol = 0;
        //std::cout<<"in: id_vi "<<vi->id_vi<<std::endl;
        if(check_updt_Viol(vi, topo)){
            actvS[vi->id_vi] = fsize+csize;
            //std::cout<<"in: "<<vi->id_vi<<std::endl;
            ++csize;
            ++num_actv;
            vi = vi->next;
        }else{
        	std::cout<<"out: "<<vi->id_vi<<std::endl;
        	++cont;
        	vi = covers.move_to_end(vi);
        }
    }
    gend=0;
    return cont;
}

//-------------------------------------------------------------------------------------------

void
CoverManager::clean_collection(){
	covers.begin = covers.end = 0;
	covers.sizeOfCollection = covers.discarted = 0;
	covers.empty = true;
	 num_actv = gend = 0;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  main methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
CoverManager::cover_generation_main(const double * ystar, const double * y,const CutSetCollection * sets, int actvSSz){
    CutSet * cutset;
    int added=0;
    int sz = sets->sizeOfCollection;
    //sets->print();
    cutset = sets->begin;
    for(int i=0;i<sz;++i){

        added+=cover_generation(cutset->ss_size, cutset->SS_arcs, cutset->uss, cutset->dss, ystar, y, actvSSz+added, cutset->id );
        //if(added>0) break;
        added+=cover_generation(cutset->s_ssize, cutset->S_Sarcs, cutset->us_s, cutset->ds_s, ystar, y, actvSSz+added, -cutset->id );
        //if(added>0) break;
        //std::cout<<"i: "<<i<<std::endl;
        cutset = cutset->next;
    }
    gend+=added;
    return added;
}

//------------------------------------------------------------------------------------------

int
CoverManager::cover_generation(int ss_size, const int * SS_arcs, double uss, double dss, const double * ystar, const double * y, int actvSSz,  int id_owner_){
    
    std::deque<Trio1> ss_;
    std::deque<Pair2> lift_down;
    std::deque<Pair2> lift_up;
    std::deque<Pair2> luc;
    double delta, u1;
    int added=0;
    CoverL* vi;
    Cover * cover;
    //std::cout<<"dss:  "<<dss<<" uss: "<<uss<<std::endl;
    u1 =  cutset_preprocess( ss_size, dss, SS_arcs,  ss_, lift_down, y, ystar);
    delta = dss - u1;
    //std::cout<<"delta in cover: "<<delta<<std::endl;

    if(delta<=0){ss_.clear();lift_down.clear();lift_up.clear(); return 0;}

    restrict_cutset(lift_down, lift_up, ss_, ystar, delta, dss, uss);
    //std::cout<<"delta in cover1: "<<delta<<std::endl;
    int id = data->nnodes*data->ndemands + covers.sizeOfCollection;
    cover = make_cover(delta, ss_, ystar, lift_down, luc, id, id_owner_ );
    //std::cout<<&luc<<std::endl;
    
    vi = coverlift.lift_cover(luc, lift_down , lift_up , int( delta), int(cover->rhs), int(dss));
    //std::cout<<"vi "<<std::endl;
    cover->addLftd(vi);
    
    //std::cout<<"done add cove"<<std::endl;
    cover->owner = SS_arcs;
    cover->maxsize = ss_size;
    cover->hs = checkViol(cover, ystar);
    if(cover->hs>0){
        added = covers.addCover(cover, ystar);
        if(added){
        	++ttgend;
        	++num_actv;
        }
        //if(added) std::cout<<"cover: "<<cover->get_total_rhs()<<std::endl;
        /*if(added){
        //std::cout<<"add cut"<<std::endl;
        //for(int a=ss_size; a--;){
         //   std::cout<<SS_arcs[a]<<":("<<data->arcs[SS_arcs[a]].i<<"-"<<data->arcs[SS_arcs[a]].j<<") y: "<<y[arc_map[SS_arcs[a]]]<<" y*: "<<ystar[arc_map[SS_arcs[a]]]<<std::endl;
        //}
        std::cout<<"cover: "<<cover->get_total_rhs()<<std::endl;
        for(int a=cover->size; a--;){
            std::cout<<cover->C[a]<<" y: "<<y[arc_map[cover->C[a]]]<<" y*: "<<ystar[arc_map[cover->C[a]]]<<std::endl;
        }
        if(cover->hasLftd){
            vi = cover->Lftd;
            for(int a=0;a<vi->size;++a){
                std::cout<<vi->SSmC[a]<<" ry: "<<y[arc_map[vi->SSmC[a]]]<<" y*: "<<ystar[arc_map[vi->SSmC[a]]]<<std::endl;
            }
        }}*/
    }else delete cover;
    
    luc.clear();
    ss_.clear();
    lift_down.clear();
    lift_up.clear();
    
    return added;
}


//-------------------------------------------------------------------------------------------

Cover *
CoverManager::make_cover(double& delta, const std::deque<Trio1> & ss_, const double * ystar, std::deque<Pair2>& lift_down, std::deque<Pair2>& cover, int id_vi, int id_owner_ ){
    int sz =(int)ss_.size();
    int arc;
    double bl;
    std::deque<Trio1> candidates;
    //std::cout<<sz<<std::endl;
    for(int n=0;n<sz;n++){
        arc = ss_[n].fst;
        bl = ss_[n].trd;
        //std::cout<<"R: "<<arc<<" ("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa/bl: "<<data->arcs[arc].capa<<"/"<<bl<<" y*: "<<ystar[arc]<<std::endl;
        if(delta<=bl){
            cover.push_back(Pair2(arc,1.0));
            //std::cout<<"c: "<<arc<<" cap: "<<bl<<" y*: "<<ystar[arc]<<std::endl;
            //std::cout<<"c: "<<arc<<" ("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa: "<<data->arcs[arc].capa<<" y*: "<<ystar[arc]<<std::endl;
        }else{
            candidates.push_back(Trio1(arc, double(ystar[arc_map[arc]]), bl));
        }
    }
    if(!candidates.empty()){
        minimize_cover(delta, candidates, cover, ystar, lift_down);
    }
    candidates.clear();
    return covers.createNewCover(cover, delta,  id_vi, id_owner_, ttgend);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  minimal Cover heuristic methods
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------


void
CoverManager::minimize_cover(double& delta, std::deque<Trio1>& candidates, std::deque<Pair2>& cover, const double * ystar,
                                std::deque<Pair2>& lift_down){
    int arc;
    double bl;
    std::stable_sort(candidates.begin(),candidates.end(),compTrio12());//decreasing order
    //std::cout<<"minimize cover delta:"<<delta<<std::endl;
    while(!candidates.empty()){
        arc = candidates.front().fst;
        bl =  candidates.front().trd;
        //std::cout<<"+c1: "<<arc<<" ("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa: "<<data->arcs[arc].capa<<" y*: "<<ystar[arc]<<std::endl;
        if(delta<=bl){
            cover.push_back(Pair2(arc, 1.0));
            //std::cout<<"c: "<<arc<<" ("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa: "<<bl<<" y*: "<<ystar[arc]<<std::endl;
        }else{delta -= bl; lift_down.push_back(Pair2(arc, double(ystar[arc_map[arc]])));
            //std::cout<<"n/c: "<<arc<<" ("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa: "<<bl<<" y*: "<<ystar[arc]<<std::endl;
        }
        candidates.pop_front();
    }
}

//-------------------------------------------------------------------------------------------

double
CoverManager::cutset_preprocess(int sz, double dss, const int * ss_,  std::deque<Trio1>& ss_deque, std::deque<Pair2> & lift_down,
                            const double *y, const double *ystar){
    int arc, id_arc;
    double U1=0;
    double bl;
    //std::cout<<"cut: "<<dss<<std::endl;
    
    for(int i=sz;i--;){
        arc = ss_[i];
        id_arc = arc_map[arc];
        bl = fmin(dss,data->arcs[arc].capa);
        //std::cout<<"arc "<<arc <<<<" bl: "<<bl<<std::endl;
        
        if(id_arc>=0){
            if(y[id_arc]>0.9){
                lift_down.push_back(Pair2(arc, double(ystar[id_arc])));
                U1+=bl;
                //std::cout<<"n/c: "<<arc<<" "<<y[id_arc]<<" cap: "<<bl<<" y*: "<<ystar[id_arc]<<std::endl;
            }else{
                ss_deque.push_back(Trio1(arc, double(ystar[id_arc]), bl));
            }
        }else if(id_arc == -2){ U1+=bl;}
    }
    return U1;
}

//-------------------------------------------------------------------------------------------

void
CoverManager::restrict_cutset(std::deque<Pair2> & lift_down, std::deque<Pair2> & lift_up, std::deque<Trio1> & ss_,
                              const double *ystar,double & delta, double dss, double uss){
    
    //for(int i=ss_.size();i--;) std::cout<<"ss_: "<<ss_[i].fst<<std::endl;
    
    form_c1(lift_down, ss_, ystar, delta, dss, uss);
    //form_c0(lift_up, ss_, ystar, dss, uss);
    
    //std::cout<<std::endl;
}

//-------------------------------------------------------------------------------------------

void
CoverManager::form_c1(std::deque<Pair2> & lift_down, std::deque<Trio1> & ss_,
                      const double *ystar,double & delta, double dss, double uss){
    int arc, id_arc;
    double capa;
    std::deque<Trio1> aux;
    std::stable_sort(ss_.begin(),ss_.end(),compTrio1());//decreasing order
    //std::cout<<"form_c1: "<<std::endl;
    arc = ss_.front().fst;
    id_arc = arc_map[arc];
    while(!ss_.empty() && ystar[id_arc]>=0.9){
        capa = ss_.front().trd;
        if(delta-capa<=0){
            aux.push_back(Trio1(arc, double(ystar[id_arc]), capa));
            ss_.pop_front();
            arc = ss_.front().fst;
            continue;
            //break;
        }
        lift_down.push_back(Pair2(arc, double(ystar[id_arc])));
        delta-= capa;
        ss_.pop_front();
        //std::cout<<"n/c "<<arc<<" :("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa: "<<capa<<" y*: "<<ystar[arc]<<std::endl;
        arc = ss_.front().fst;
    }
    //std::cout<<"stop : delta =  "<<delta<<std::endl;
    ss_.insert(ss_.end(), aux.begin(), aux.end());
    aux.clear();
}

//-------------------------------------------------------------------------------------------

void
CoverManager::form_c0(std::deque<Pair2> & lift_up, std::deque<Trio1> & ss_,
                      const double *ystar, double dss, double uss){
    int arc, id_arc;
    double capa;
    std::deque<Trio1> aux;
    std::stable_sort(ss_.begin(),ss_.end(),compTrio2());//decreasing order
    //std::cout<<"form_c0: "<<std::endl;
    while(!ss_.empty()){
        arc = ss_.front().fst;
        capa  = ss_.front().trd;
        id_arc = arc_map[arc];

        if( ystar[id_arc]>=1.0){
            aux.push_back(Trio1(arc, double(ystar[id_arc]), capa));
            ss_.pop_front();
            continue;
            //break;
        }
        if(uss-capa <dss){
            //std::cout<<"stop :"<<uss<<" <=> "<<dss<<std::endl;
            break;
        }
        lift_up.push_back(Pair2(arc,double(ystar[id_arc])));
        uss-=capa;
        //std::cout<<arc<<" :("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<") capa: "<<data->arcs[arc].capa<<" y*: "<<ystar[arc]<<" **: "<<ss_.front().trd<<std::endl;
        ss_.pop_front();
    }
    ss_.insert(ss_.end(), aux.begin(), aux.end());
    aux.clear();
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

double
CoverManager::checkViol(const Cover * c, const double *y){
    
    double sum=0;
    double comp=c->get_total_rhs();
    int sz = c->get_total_sz();
    int id_arc;
    for(int a=0;a<sz;++a){
        id_arc = arc_map[c->at(a)];
        if(id_arc>=0){
            sum+= c->gamma_at(a)*y[id_arc];
        }
        if(sum>=comp){return 0.0;}
    }
    return comp-sum;
}

//-------------------------------------------------------------------------------------------

bool
CoverManager::check_updt_Viol(Cover * c, const double *y){
    
    double sum=0;
    double comp=c->get_total_rhs();
    int sz = c->get_total_sz();
    int id_arc;
    c->rhs_dimsh=0;
    for(int a=0;a<sz;++a){
    	id_arc = c->at(a);
        //std::cout<<"arc: "<<id_arc<<": "<<y[id_arc]<<std::endl;
        sum+= c->gamma_at(a)*y[id_arc];
        if(sum>=comp){return false;}
    }
    c->rhs_dimsh = sum;
    return true;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  Multipliers methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

double
CoverManager::update_dual_pos( const double * dstar, std::vector<Trio1>& ws, double * dual, double * h){
    int idx;
    int size = ws.size();
    for(int n=size; n--;){
        idx = ws[n].fst;
        if(idx<0){continue; }
        //std::cout<<idx<<" oldalph: "<<dual[idx]<<" newalph:  "<<ws[n].snd<<std::endl;
        dual[idx] += ws[n].snd;
    }
    return 0;
}

//-------------------------------------------------------------------------------------------

double
CoverManager::set_new_mult_pos(double *rc, std::vector<Trio1>& ws, const double * dual,
                            const std::deque<int>& con_arcs, const std::vector<int> & con_arcs_map,
                            const std::vector<Pair2>&  con_arcs_wnid){
    int arc, idc, index, id_arc;
    int size = num_actv;
    int tt;
    double sum, ret, diff, div;
    double alphsum;
    Pair2 item;
    for(int a=con_arcs.size(); a--;){
        arc = con_arcs[a];
        id_arc = arc_map[arc];
        tt = con_arcs_map[arc];
        if(rc[id_arc]<0){ continue;}
        sum=0; alphsum=0;
        for(int c = 0;c<tt;c++){
            item = con_arcs_wnid[arc*size+c];
            sum += ws[item.fst].snd * item.snd;
            index = ws[item.fst].fst;
            alphsum += dual[index];
        }
        diff = sum - rc[id_arc];
        //std::cout<<"arc: "<<arc<<" sum: "<<sum<<" dif: "<<diff<<std::endl;
        if(diff>1e-10){
            div = rc[id_arc]/double(tt);
            for(int c = 0;c<tt;c++){
                item = con_arcs_wnid[arc*size+c];
                index = ws[item.fst].fst;
                if(alphsum > 1e-10) div = rc[id_arc]*(dual[index]/(double)(alphsum*item.snd));
                else div = div/(double)(item.snd);
                //std::cout<<std::setprecision(15)<<"recomp: "<<index<<" : "<<rc[arc]<<" * "<<dual[index]<<" / "<<alphsum<<" w: "<<ws[item.fst].snd<<" div: "<<div<<std::endl;
                ws[item.fst].snd = fmin(div,ws[item.fst].snd);
            }
        }
    }
    for(int a=con_arcs.size(); a--;){
        arc = con_arcs[a];
        tt = con_arcs_map[arc];
        id_arc = arc_map[arc];

        sum=0;
        for(int c = 0;c<tt;c++){
            item = con_arcs_wnid[arc*size+c];
            sum += ws[item.fst].snd * item.snd;
        }
        //std::cout<<arc<<" new rc: "<<rc[arc]<<std::endl;
        if(rc[id_arc] - sum < -1e-3  && rc[id_arc] >0)std::cout<<"PROBLEM arc "<<arc<<" : "<<rc[arc]<<" < "<<sum<<std::endl;
    }
    return 0;
}

//-------------------------------------------------------------------------------------------

double
CoverManager::update_rc_neg(double dimsh, int nvi, const Cover * vi, std::vector<Trio1>& ws, std::vector<Pair2> & con_arcs_map, double *rc){
    int arc, id_arc;
    //std::cout<<"confirm: "<<dimsh<<std::endl;
    for(int a=vi->get_total_sz(); a--;){
        arc = vi->at(a);
        id_arc = arc_map[arc];
		if(id_arc < 0) continue;
        rc[id_arc] -= vi->gamma_at(a)*dimsh;
        con_arcs_map[arc].snd -= ws[nvi].snd;
    }
    ws[nvi].snd += dimsh;
    return 0;
}

//-------------------------------------------------------------------------------------------

double
CoverManager::update_con_arcs(std::deque<Pair2>& con_arcs, const double* fk, const double *rc){
    int arc, id_arc;
    for(int a = con_arcs.size();a--;){
        arc = con_arcs[a].fst;
        id_arc = arc_map[arc];

        if(fk[id_arc]<=0){
            con_arcs[a].snd = rc[id_arc] - fk[id_arc];
        }else con_arcs[a].snd  =  rc[id_arc];
    }
    return 0;
}


//-------------------------------------------------------------------------------------------

double
CoverManager::set_new_mult_neg(double *rc, std::vector<Trio1>& ws,  double * dual, std::vector<Pair2> & con_arcs_map,
                                std::deque<Pair2>& con_arcs, const std::vector<Pair2>& con_arcs_wnid,
                                const std::vector<const Cover *>& addrs,const double *fk){
    
    int arc, index, tt, id_arc;
    int size = num_actv;
    int szarcs = con_arcs.size();
    double diff, mult, alphsum;
    Pair2 item;
    
    std::stable_sort(con_arcs.begin(), con_arcs.end(), compPair2());
    while(!con_arcs.empty()){
        arc = con_arcs.front().fst;
        id_arc = arc_map[arc];
        mult = con_arcs.front().snd;
        con_arcs.pop_front();
        tt = con_arcs_map[arc].fst;
        alphsum = con_arcs_map[arc].snd;

        //std::cout<<"arc: "<<arc<<" fk: "<<fk[id_arc]<<" rc: "<<rc[id_arc]<<" mult: "<<mult<<" tt: "<<tt<<" al: "<<alphsum<<std::endl;
        if(mult > -1e-10 || alphsum  < 1e-10 ) continue;

        for(int c=0; c<tt; c++){
            item = con_arcs_wnid[arc*size+c];
            if(ws[item.fst].trd == 1.0){  continue; }
            index = ws[item.fst].fst;
            diff = mult*(dual[index]/(double)(alphsum*item.snd));
            
            if(alphsum < 1e-10) std::cout<<"PROBLEMA set_new_mult_neg alphsum < 1e-10"<<std::endl;
            if(dual[index]+diff<0){ diff = -dual[index]; }//std::cout<<"PROBLEMA set_new_mult_neg dual[index]<0 "<<dual[index]<<std::endl;}
            ws[item.fst].trd = 1.0;
            dual[index] += diff;
            update_rc_neg(diff, item.fst, addrs[item.fst], ws, con_arcs_map, rc);

        }
        
		
        if(rc[id_arc]<1e-10 && rc[id_arc]>-1e-10) rc[id_arc] = 0;
        
        update_con_arcs( con_arcs, fk, rc);
        std::stable_sort(con_arcs.begin(), con_arcs.end(), compPair2());
        //std::cout<<"in arc: "<<arc<<" new rc: "<<rc[arc]<<std::endl;
    }
    return 0;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  Volume Integration methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

double
CoverManager::recompute_mult_neg(double * dual, double * rc, const double *fk,
                                 const double *xy, const int * actvS, int actvSSz){
    
    if(num_actv==0) return 0;
    int index, arc, idc, id_arc;
    int narcs = data->narcs;
    int size = num_actv;
    double alph, viol, min;
    Trio1 item;
    //std::cout<<"recompute_mult_neg"<<std::endl;
    std::vector<Pair2> con_arcs_map(narcs,Pair2(0,0));
    std::vector<Pair2> con_arcs_wnid(narcs*size,Pair2(-1,-1));
    std::vector<Trio1> ws(size,Trio1(-1,-1,-1));
    std::vector<const Cover *> addrs(size,0);
    std::deque<Pair2> con_arcs;
    std::deque<Pair2> con_arcs_aux;
    
    Cover * vi = covers.begin;
    for(int n=size; n--;){
        addrs[n] = vi;
        index = actvS[vi->id_vi];
        if(checkViol(vi, xy)){ ws[n].fst = -index; vi = vi->next; continue; }

        for(int a=vi->get_total_sz(); a--;){
            arc = vi->at(a);
            id_arc = arc_map[arc];
            if(id_arc < 0) continue;
            alph = rc[id_arc];
            if(alph>=0) continue;
            if(con_arcs_map[arc].fst==0){
                if(fk[id_arc]<=0) con_arcs_aux.push_back(Pair2(arc, alph - fk[id_arc]));
                else con_arcs.push_back(Pair2(arc, alph));
            }
            con_arcs_map[arc].snd += dual[index];
            idc = con_arcs_map[arc].fst++;
            con_arcs_wnid[arc*size+idc].fst = n;
            con_arcs_wnid[arc*size+idc].snd = vi->gamma_at(a);
        }

        item.fst = index;
        item.snd= dual[index];
        item.trd = 0.0;
        ws[n] = item;
        vi = vi->next;
    }
    
    set_new_mult_neg( rc,  ws,  dual, con_arcs_map, con_arcs, con_arcs_wnid, addrs, fk);
    update_con_arcs( con_arcs_aux, fk, rc);
    set_new_mult_neg( rc,  ws,  dual, con_arcs_map, con_arcs_aux, con_arcs_wnid, addrs, fk);

    con_arcs_map.clear();
    con_arcs_wnid.clear();
    ws.clear();
    addrs.clear();
    
    return 0;
}

//-------------------------------------------------------------------------------------------


double
CoverManager::recompute_mult_pos(double * dual, double * h,  double *rc,
                             const double * dstar, const double *xy, const int * actvS){
    if(num_actv==0) return 0;
    int index, arc, idc, id_arc;
    int narcs = data->narcs;
    int size = num_actv;
    double alph, viol, min;
    
    std::vector<int> con_arcs_map(narcs,0);
    std::vector<Pair2> con_arcs_wnid(narcs*size,Pair2());
    std::deque<int> con_arcs;
    std::vector<Trio1>ws(size,Trio1(-1,-1,-1));
    Cover * vi = covers.begin;
    //std::cout<<"recompute_mult "<<size<<std::endl;
    for(int n=size; n--;){
        index = actvS[vi->id_vi];
        if(!checkViol(vi, xy)){ ws[n].fst = -index; vi = vi->next; continue; }
        //std::cout<<"pos c: "<<n<<" rhs: "<<vi->get_total_rhs()<<" dual: "<<dual[index]<<" h: "<<h[index]<<std::endl;
        min = 1e30;
        for(int a=vi->get_total_sz(); a--;){
            //std::cout<<"qa "<<vi->id_vi<<std::endl;
            arc = vi->at(a);
            id_arc = arc_map[arc];
            if(id_arc < 0) continue;
            //std::cout<<"arc : "<<arc<<" id: "<<id_arc<<std::endl;
            alph = rc[id_arc]; //std::cout<<"alph : "<<alph<<std::endl;
            if(alph<0) continue;

            if(con_arcs_map[arc]==0)con_arcs.push_back(arc);
            if(min > alph ) min = alph;
            idc = con_arcs_map[arc]++;
            con_arcs_wnid[arc*size+idc].fst = n;
            con_arcs_wnid[arc*size+idc].snd = vi->gamma_at(a);
        }

        //else std::cout<<"new idx: "<<idx<<" viol: "<<h[idx]<<std::endl;
        ws[n].fst = index;
        ws[n].snd = min;
        ws[n].trd = vi->get_total_rhs();
        //std::cout<<"min: "<<min<<std::endl;
        vi = vi->next;
    }
    //std::cout<<"ok"<<std::endl;
    set_new_mult_pos(rc,  ws, dual, con_arcs, con_arcs_map, con_arcs_wnid);
    double ret = update_dual_pos(  dstar,  ws, dual, h);
    con_arcs_map.clear();
    con_arcs_wnid.clear();
    con_arcs.clear();
    ws.clear();
    
    return ret;
}


//-------------------------------------------------------------------------------

void
CoverManager::add_cover_vi(int added, int * actvS, int & actvSSz, double * h, double * dual, double * dual_lb, double * dual_ub ){
    //std::cout<<"add_flowpack_vi: "<<actvSSz<<std::endl;
    
    int idx;
    Cover * vi = covers.end;
    for(int cont = added; cont--;){
        idx = actvSSz+cont; if(actvS[vi->id_vi]>=0) std::cout<<actvS[vi->id_vi]<<" already taken !!!!!!! for: "<<vi->id_vi<<std::endl;
        actvS[vi->id_vi]=idx;
        dual[idx] =0;
        dual_lb[idx] = 0;
        dual_ub[idx] = 1e31;
        h[idx] = vi->hs;
        //std::cout<<"add vi: "<< vi->id_vi<<" idx: "<<idx<<" cont: "<<cont<<" h: "<<h[idx]<<std::endl;
        //vi->print();
        vi = vi->prev;
    }
    actvSSz += added;
    //std::cout<<" size: "<<actvSSz<<std::endl;
    //std::cout<<"pass add_flowpack_vi"<<std::endl;
}

//-------------------------------------------------------------------------------

int
CoverManager::compute_cover_sg( const double * x, const int * actvS, int actvSSz,  double * v){
    //std::cout<<"compute_flowpack_sg"<<std::endl;
    int index, id_arc;
    int sz = num_actv;
    Cover *vi = covers.begin;
    CoverL *lifted;
    for(int n=0;n<sz;++n){
        index = actvS[vi->id_vi];
        
        v[index] = vi->get_total_rhs();
        for(int a=vi->get_total_sz();a--;){
            id_arc = arc_map[vi->at(a)];
            if(id_arc<0) continue;
            v[index] -=  vi->gamma_at(a)*x[id_arc];

            //if(c->C[a]==74)std::cout<<"in c: "<<c->C[a]<<" : "<<y[c->C[a]]<<std::endl;
            //std::cout<<i<<" inc: "<<cover->C[a]<<" : "<<x[cover->C[a]]<<std::endl;
        }

        if(v[index]<=0){
            ++vi->n_nviol;
            if(vi->n_nviol>=lim_to_remv && vi->n_zerom>0) v[index]=0;
        }else vi->n_nviol=0;
        
        vi = vi->next;
    }
    //std::cout<<"passou compute_flowpack_sg"<<std::endl;
    
    return 0;
}

//-------------------------------------------------------------------------------

int
CoverManager::compute_cover_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0){
    //std::cout<<" compute_flowpack_rc"<<std::endl;
    int index, id_arc;
    int sz = num_actv;
    CoverL *lifted;
    Cover *vi = covers.begin;
    //std::cout<<"bB0: "<<B0<<std::endl;
    for(;sz--;){
        
        index = actvS[vi->id_vi];
        if(dual[index]==0){
            ++vi->n_zerom;
            vi = vi->next;
            continue;
            
        }else vi->n_zerom = 0; 
        //std::cout<<"mu: "<<vi->mu<<std::endl;

        //std::cout<<"vi_id: "<<vi->id_vi<<" srnbr: "<<vi->serial_nmbr<<" rhs: "<<vi->get_total_rhs()<<" idx: "<<index<<std::endl;
        B0 +=  dual[index]*vi->get_total_rhs();
        for(int a=vi->get_total_sz();a--;){
            id_arc = arc_map[vi->at(a)];
            //std::cout<<"arc"<<vi->at(a)<<" id_arc: "<<id_arc<<" gam: "<<vi->gamma_at(a)<<std::endl;
            if(id_arc<0) continue;
            rc[id_arc]-=  vi->gamma_at(a)*dual[index];
            //if(c->C[a]==74)std::cout<<"in c: "<<c->C[a]<<" : "<<y[c->C[a]]<<std::endl;
            //std::cout<<i<<" inc: "<<cover->C[a]<<" : "<<x[cover->C[a]]<<std::endl;
        }
        
        
        //std::cout<<"idx: "<<index<<" "<<(vi->ymult)*dual[index]<<" "<<dual[index]<<" dual: "<<dual[index]<<std::endl;
        //std::cout<<"idx: "<<index<<" dual: "<<dual[index]<<" mu: "<<vi->mu<<" ws:"<<vi->ws<<std::endl;
        //std::cout<<"idx: "<<index<<" dual: "<<dual[index]<<" lbd_c: "<<(vi->ymult)<<" mu: "<<vi->mu<<" cid: "<<vi->owner->id<<std::endl;
        vi = vi->next;
    }
    //std::cout<<"B0: "<<B0<<std::endl;

    //std::cout<<"passou compute_flowpack_sg"<<std::endl;
    return 0;
}

//-------------------------------------------------------------------------------

void
CoverManager::make_inactive(int index, const int* actvS, double * v){
    Cover* c  = covers[index];
    v[actvS[c->id_vi]] = 0;
    ++c->n_nviol;
    //std::cout<<"make_inactive index: "<<actvS[c->id_vi]<<std::endl;

}

//-------------------------------------------------------------------------------

double
CoverManager::arc_dg_imp(int arc, const double * xy, const double * h, const int * actvS, int actvSSz){
    int index;
    int sz = num_actv;
    CoverL *lifted;
    Cover *vi = covers.begin;
    double gam;
    double dg=0;
    for(;sz--;){
        index = actvS[vi->id_vi];
        gam = covers.cover_hasArc(vi, arc);
        dg += -gam*h[index];
        vi = vi->next;
    }
    return dg;
}
