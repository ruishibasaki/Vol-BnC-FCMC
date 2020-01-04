#include "lift.hpp"

#define MAX_INT  1000000000;

//=========================================================================================
//=========================================================================================
// Lift Class Methods
//=========================================================================================
//=========================================================================================

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// setters and initializing methods / Auxiliars
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void 
Lift::print(const std::deque<Pair2>& luc, int dbar){
    for(int i=(int)luc.size();i--;) std::cout<<luc[i].fst<<" : "<<luc[i].snd<<" capa: "<<min(dss, data->arcs[luc[i].fst].capa)<<std::endl;
    std::cout<<"dbar: "<<dbar<<std::endl;

}

//-----------------------------------------------------------------------------------------

void 
Lift::set_data(const Data * d){
    data = d;
    nnodes = data->nnodes;
    ndemands = data->ndemands;
    narcs = data->narcs;
    p = new int [narcs];
    w = new int [narcs];
    sz_states = 0;
    max_dbar = luc_tt_u = max_n = max_u = dim =0;
    min_dbar = MAX_INT;
    //knapslvr.set_data(d);

}

//-----------------------------------------------------------------------------------------

int
Lift::reset( const std::deque<Pair2> & lift_down, int dbar){
    clear();
    max_dbar = dbar;
    for(int i=(int)lift_down.size();i--;)
        max_dbar += min(dss , int(data->arcs[lift_down[i].fst].capa));
    
    if(sz_states>0){
        delete [] states;
        sz_states=0;
    }
    sz_states = max_n*max_dbar;
    states = new int [sz_states];
    luc_tt_u=0;
    max_up_dbar = max_dbar;
    
    return 0;
}


//-----------------------------------------------------------------------------------------

int
Lift::build(const std::deque<Pair2>& luc,const std::deque<Pair2>& lift_up,const std::deque<Pair2>& lift_down){
    Pair2 item;
    bool down_first=true;
    int g, capa, min_up_u, diff;
    int downsize = (int)lift_down.size();
    
    dim= (int)luc.size();
    min_up_u = min_dbar = MAX_INT;
    max_u =0;
    for(int i=dim;i--;){
        item = luc[i];
        g = int(item.snd);
        capa= min(dss, int(data->arcs[item.fst].capa));
        
        p[i] = g;
        w[i] = capa;
        if(capa<min_dbar) min_dbar = capa;
        if(capa>max_u) max_u = capa;
    }
    for(int i=downsize;i--;){
        capa = min(dss, int(data->arcs[lift_down[i].fst].capa));
        if(capa<min_dbar) min_dbar = capa;
        if(capa>max_u) max_u = capa;
    }
    if(downsize==0) down_first=false;
    for(int i=(int)lift_up.size();i--;){
        capa = min(dss, int(data->arcs[lift_up[i].fst].capa));
        diff = max_dbar - capa;
        if(diff>=0 && diff<min_dbar) min_dbar = diff;
        if(capa<min_up_u) min_up_u = capa;
        if(!down_first && capa>max_u) max_u = capa;

    }
    max_up_dbar = max_dbar - min_up_u;
    if(max_up_dbar<0) max_up_dbar =0;
    if(min_dbar>max_dbar) min_dbar = max_dbar;
    return 0;
}

//-----------------------------------------------------------------------------------------

void
Lift::add_lifted_var(bool down, CoverL* vi, std::deque<Pair2>& luc, int dbar, int u, int gam, int arc, int & ttgamcul, int & ttgamc1){
    
    gam = (gam<0)? 0 :gam;
    if(down){
        ttgamc1+=gam;
    }
    ttgamcul+=gam;
    luc.push_back(Pair2(arc,gam));
    p[dim] = gam;
    w[dim] = u;
    ++dim;
    if(down && gam>1){std::cout<<"var down: "<<arc<<" lifted gam: "<<gam<<" capa: "<<u<<std::endl;}
    //for(int i=(int)luc.size();i--;) std::cout<<luc[i].fst<<" capa: "<<min(dss,int(data->arcs[luc[i].fst].capa))<<std::endl;}
    else if(!down && gam>1) std::cout<<"var up: "<<arc<<" lifted gam: "<<gam<<" capa: "<<u<<std::endl;
    
    if(gam>0){vi->addvar(arc,gam, down);}
    
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// Down Lifting methods
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


void
Lift::lift_vars_down(CoverL* vi,  std::deque<Pair2> & lift_down ,std::deque<Pair2>& luc, int & dbar, int & ttgamcul, int & ttgamc1, int rhs, bool down_first){
    int arc;
    if(down_first){
        arc = lift_down.back().fst;
        down_lifting(0, dim, arc, vi, luc, dbar, ttgamcul, ttgamc1, rhs);
        lift_down.pop_back();
    }
    while(!lift_down.empty()){
        arc = lift_down.back().fst;
        down_lifting(dim-1, dim, arc, vi, luc, dbar, ttgamcul, ttgamc1, rhs);
        lift_down.pop_back();
    }
}

//-----------------------------------------------------------------------------------------


int
Lift::down_lifting(int strt, int end,int arc , CoverL* vi , std::deque<Pair2>& luc, int & dbar, int & ttgamcul, int & ttgamc1, int rhs){
    
    int z,  gam;
    int u;
    u = min(dss, int(data->arcs[arc].capa));
    z = solve( u, dbar, strt,  end, true);
    //double U=dbar+u; double Z=0; int zz = knapslvr.solve(luc,U,Z); if(Z!=z){print(luc,U);std::cout<<"PRRROBLEMMM "<<Z<< "  dwn "<<z<<std::endl;}
    //if(zz!=1){
    if(z<0){
        gam = ttgamcul - ttgamc1 -rhs +1.0;
    }else{
        gam = z /*Z*/ - rhs - ttgamc1;
    }
    //std::cout<<"solution "<<z<<std::endl;
    add_lifted_var(true, vi,luc, dbar, u, gam, arc, ttgamcul, ttgamc1);
    dbar+=u;
    return gam;
}




//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// Up Lifting methods
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void
Lift::lift_vars_up(CoverL* vi,  std::deque<Pair2> & lift_up ,std::deque<Pair2>& luc, int& dbar, int& ttgamcul, int& ttgamc1, int rhs, bool down_first){
    int arc;
    if(!down_first){
        arc = lift_up.front().fst;
        up_lifting(0, dim, arc, vi, luc, dbar, ttgamcul, ttgamc1, rhs);
        lift_up.pop_front();
    }
    while(!lift_up.empty()){
        arc = lift_up.front().fst;
        up_lifting(dim-1, dim, arc ,  vi , luc, dbar, ttgamcul,  ttgamc1,  rhs);
        lift_up.pop_front();
    }
}

//-----------------------------------------------------------------------------------------

int
Lift::up_lifting(int strt, int end,int arc , CoverL* vi , std::deque<Pair2>& luc, int & dbar, int & ttgamcul, int & ttgamc1, int rhs){
    
    int z,  gam;
    int u;
    
    u = min(dss, int(data->arcs[arc].capa));
    z = solve( -u, dbar, strt,  end, false);
    //double U=dbar-u; double Z=0; int zz = knapslvr.solve(luc,U,Z); if(Z!=z){print(luc,U);std::cout<<"PRRROBLEMMM "<<Z<< " up "<<z<<std::endl;}
    //if(z!=1){
    if(z<0){
        gam = ttgamc1 - ttgamcul + rhs - 1.0;
    }else{
        gam = rhs - z /*Z*/ + ttgamc1;
    }
    //std::cout<<"solution "<<z<<std::endl;
    add_lifted_var(false, vi,luc, dbar, u, gam, arc, ttgamcul, ttgamc1);
    return gam;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
// Main methods
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


CoverL*
Lift::lift_cover(std::deque<Pair2>& luc, std::deque<Pair2> & lift_down , std::deque<Pair2> & lift_up , int dbar, int rhs, int dss_){
    bool down_first = false;
    dss = dss_;
    //int pri=0;
    //luc.clear(); lift_down.clear(); lift_up.clear();
    //luc.push_back(Pair2(104,1)); luc.push_back(Pair2(108,1)); luc.push_back(Pair2(233,1));
    //for(int i=80;i<70;++i) luc.push_back(Pair2(i,1)); //luc[0].snd=1;
    //for(int i=74;i<78;i++) lift_down.push_front(Pair2(i,0));
    //for(int i=78;i<84;i++) lift_up.push_front(Pair2(i,0));
    //dbar = 2336;
    //std::cout<<"---------- START LIFTING ---------- "<<std::endl;
    //std::cout<<"delta: "<<dbar<<std::endl;
    
    //for(int i=luc.size();i--;) std::cout<<pri++<<" : "<<luc[i].snd<<" capa: "<<data->arcs[luc[i].fst].capa<<std::endl;
    int ttgamcul, ttgamc1;
    
    CoverL* vi = 0;
    max_n = (int) (lift_down.size()+lift_up.size());
    //std::cout<<"delta: "<<dbar<<std::endl;

    if(dbar<0) return vi;
    if(max_n==0)return vi;
    vi = new CoverL(max_n);
    max_n += luc.size();
    
    reset(lift_down,  dbar);
    build(luc, lift_up, lift_down);
    //std::cout<<"MAX_U "<<max_u<<std::endl;
    //std::cout<<"MIN_DBAR: "<<min_dbar<<std::endl;
    //std::cout<<"MAX_DBAR: "<<max_dbar<<std::endl;
    //std::cout<<"MAX_UP_DBAR: "<<max_up_dbar<<std::endl;
    
    ttgamcul= dim;
    ttgamc1 = 0;
    if(!lift_down.empty()){
        down_first = true;
        std::stable_sort(lift_down.begin(),lift_down.end(),compPair2());//decreasing order
        lift_vars_down(vi, lift_down, luc, dbar, ttgamcul, ttgamc1, rhs, down_first);
    }
    //std::cout<<"LIFT UP "<<std::endl;//pri=0;
    //for(int i=luc.size();i--;) std::cout<<pri++<<" : "<<luc[i].snd<<" capa: "<<data->arcs[luc[i].fst].capa<<std::endl;
    if(!lift_up.empty()){
        std::stable_sort(lift_up.begin(),lift_up.end(),compPair2());//decreasing order
        lift_vars_up(vi, lift_up, luc, dbar, ttgamcul, ttgamc1, rhs, down_first);
    }
    return vi;
}

//-----------------------------------------------------------------------------------------

int
Lift::solve( int u, int  dbar, int strt, int end, bool down){
    int z;
    dbar+=u;
    //std::cout<<"solve dbar: "<<dbar<< " dim: "<<dim<<std::endl;

    if(strt==end){
        if(dbar<=0) return 0;
        else return -1;
    }
    z = dpsolver(dbar, strt, end, down);
    
    return z;
}

//-----------------------------------------------------------------------------------------

int
Lift::dpsolver(int dbar, int strt, int end, bool down){
    
    int c, v;
    int cmax_up = min (max_up_dbar , max_dbar);
    int cmax, cmaxt;
    int n=0;
    if(strt==0){
        luc_tt_u = w[0];
        v = min(w[0],max_dbar);
        if(down){ cmax = min( luc_tt_u+ max_u ,max_dbar);}
        else{ cmax = cmax_up; }
        //std::cout<<"cmax: "<<cmax<<std::endl;
        
        for (c = 0; c < v; c++){
            states[c] = p[0];
        }for (c = v; c < cmax; c++){
            states[c] = MAX_INT;
        }
        ++strt;
    }

    for (int i = strt; i < end; i++){
        v = min(w[i],max_dbar); //std::cout<<i<<" id: "<<i<<" v: "<<v<<" p: "<<p[i]<<std::endl;
        luc_tt_u += w[i];
        if(down){ cmax = min(luc_tt_u + max_u ,max_dbar); cmaxt = min(luc_tt_u, max_dbar);}
        else{ cmax = cmax_up; cmaxt = min(luc_tt_u, cmax_up);}
        //std::cout<<"cmax: "<<cmax<<" cmaxt: "<<cmaxt<<std::endl;
        
        for (c = 0; c < v; c++){
            states[i*max_dbar+c] =  min(p[i], states[(i-1)*max_dbar+c]);
        }
        for (c = v; c < cmaxt; c++){
            states[i*max_dbar+c] = min(states[(i-1)*max_dbar+c], states[(i-1)*max_dbar+(c-w[i])] + p[i]);
        }
        for (c = cmaxt; c < cmax; c++){
            states[i*max_dbar+c] = MAX_INT;
        }
        n =i;
    }
    
    int val=0;
    if(dbar<=0) val = 0;
    else if(luc_tt_u < dbar) val = -1;
    else val = states[(n)*max_dbar+dbar-1];
    
    return val;
}

