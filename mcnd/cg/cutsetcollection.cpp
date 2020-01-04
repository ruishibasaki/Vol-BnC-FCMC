//
//  cutsetcollection.cpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#include "cutsetcollection.hpp"



//====================================================================================
//====================================================================================
// Class CutSet Collection methods
//====================================================================================
//====================================================================================

void
CutSetCollection::initialize(int M){
    int sizeOfInt=8*sizeof(unsigned int);
    int mirror_limit = (M%sizeOfInt);
    
    last_mirror = 0;
    sizeOfContainer = (M/sizeOfInt)+1;
    sizeOfMap = M;
    sizeOfCollection = 0;
    map =new Pair2[M];
    empty=true;
    for(int i=0;i<M;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
        //std::cout<<i<<" "<<map[i].fst<<" : "<<map[i].snd<<std::endl;
    }
    for(int i=0;i<mirror_limit;++i){
        last_mirror |= (1) << (i);
    }
    discarted=0;
    end = begin =0;
}

//----------------------------------------------------------------------------------

CutSetCollection::~CutSetCollection(){
    for(int i=0;i<sizeOfCollection;++i){
        CutSet* next = begin->next;
        delete begin;
        begin = next;
    }
    if(map) delete [] map;
}

//----------------------------------------------------------------------------------

int
CutSetCollection::addCutSet(CutSet * newSS_, const std::deque<int>& ss_arcs, double uss, double dss, const std::deque<int>& s_sarcs, double us_s, double ds_s){
    int ssz =(int) ss_arcs.size();
    int s_sz = (int)s_sarcs.size();
    newSS_->ss_size = ssz;
    newSS_->s_ssize = s_sz;
    newSS_->SS_arcs = new int [ssz];
    newSS_->S_Sarcs = new int [s_sz];
    
    newSS_->dss = dss;
    newSS_->uss = uss;
    newSS_->ds_s = ds_s;
    newSS_->us_s = us_s;
    for(int n=ssz;n--;){
        newSS_->SS_arcs[n] = ss_arcs[n];
    }
    for(int n=s_sz;n--;){
        newSS_->S_Sarcs[n] = s_sarcs[n];
    }
    
    ++sizeOfCollection;
    newSS_->id = sizeOfCollection;
    if(empty){
        end = begin = newSS_;
        empty=false;
    }else{
        end->next= newSS_;
        end = newSS_;
    }
    //for(int i=sizeOfMap;i--;)
    //    if(newSS_->nodeInS(map[i].fst,map[i].snd))std::cout<<i<<" ";
    //std::cout<<std::endl;
    return 1;
}
//----------------------------------------------------------------------------------

CutSet *
CutSetCollection::collect(const std::vector<int>& s_nodes){
    int sz = (int)s_nodes.size();
    CutSet * newSS_ = new CutSet(sizeOfContainer);
    
    for(int n=sz;n--;){
        if(s_nodes[n]) newSS_->addNode(map[n].fst, map[n].snd);
    }
    if(collected(newSS_)){++discarted;delete newSS_; return 0;}
    return newSS_;
}


//----------------------------------------------------------------------------------
const bool
CutSetCollection::collected(const CutSet * trySS_){
    bool pos;
    unsigned int tld;
    CutSet* SS_ = begin;
    for(int i=0;i<sizeOfCollection;++i){
        if((trySS_->s_size == SS_->s_size) || (trySS_->s_size == (sizeOfMap-SS_->s_size)) ){
            pos = true;
            
            for(int id=0;id<sizeOfContainer;++id){
                if(id==sizeOfContainer-1) tld = last_mirror ^ trySS_->SS_nodes[id];
                else tld = ~trySS_->SS_nodes[id];

                if((trySS_->SS_nodes[id]!=SS_->SS_nodes[id]) && (tld!=SS_->SS_nodes[id])){
                    pos=false;
                    break;
                }
            }
            if(pos){
                return true;
            }
        }
        SS_ = SS_->next;
    }
    return false;
}


//----------------------------------------------------------------------------------

void
CutSetCollection::advance(CutSet*& SS_, int n){
    if(n>sizeOfCollection) return;
    SS_ = begin;
    for(int i=1;i<=n;++i)
        SS_ = SS_->next;
}

//----------------------------------------------------------------------------------

void
CutSetCollection::print() const{
    CutSet* c = begin;
    for(int i=0;i<sizeOfCollection;++i){
        std::cout<<"set "<<i<<" "<<std::endl;
        for(int id=0;id<sizeOfContainer;++id){
            std::cout<<c->SS_nodes[id]<<" ";
        }
        std::cout<<std::endl;
        for(int id=0;id<c->ss_size;++id){
            std::cout<<c->SS_arcs[id]<<" ";
        }
        std::cout<<" reverse: "<<std::endl;
        for(int id=0;id<c->s_ssize;++id){
            std::cout<<c->S_Sarcs[id]<<" ";
        }
        std::cout<<std::endl;
        c = c->next;
    }
}




//====================================================================================
//====================================================================================
// Class CutSet methods
//====================================================================================
//====================================================================================

void CutSet::addNode(int iset, int node){
    ++s_size;
    setBit(SS_nodes,iset,node);
}

//----------------------------------------------------------------------------------

void CutSet::removeNode(int iset, int node){
    clearBit(SS_nodes,iset,node);
}

//----------------------------------------------------------------------------------

bool CutSet::nodeInS(int iset, int node) const{
    if(testBit(SS_nodes,iset,node)) return true;
    else return false;
}

//----------------------------------------------------------------------------------

void
CutSet::copyref(bool stosb, int& sb_ssz, int& ssb_sz, int *& sb_sarcs, int *& ssb_arcs, double &dss_ ) const{
    if(stosb){
        ssb_sz = ss_size;
        sb_ssz = s_ssize;
        ssb_arcs = SS_arcs;
        sb_sarcs = S_Sarcs;
        dss_ = dss;
    }else{
        ssb_sz = s_ssize;
        sb_ssz = ss_size;
        ssb_arcs = S_Sarcs;
        sb_sarcs = SS_arcs;
        dss_ = ds_s;
    }
}

//----------------------------------------------------------------------------------


void
CutSet::copyref(bool stosb, int& ssb_sz, int *& ssb_arcs, double &dss_ ) const{
    if(stosb){
        ssb_sz = ss_size;
        ssb_arcs = SS_arcs;
        dss_ = dss;
    }else{
        ssb_sz = s_ssize;
        ssb_arcs = S_Sarcs;
        dss_ = ds_s;
    }
}


//----------------------------------------------------------------------------------


void
CutSet::copyref(bool stosb, int& ssb_sz, int *& ssb_arcs, double &dss_, double &uss_ ) const{
    if(stosb){
        ssb_sz = ss_size;
        ssb_arcs = SS_arcs;
        dss_ = dss;
        uss_ = uss;
    }else{
        ssb_sz = s_ssize;
        ssb_arcs = S_Sarcs;
        dss_ = ds_s;
        uss_ = us_s;
    }
}

//====================================================================================
//====================================================================================
// MinCardCS methods
//====================================================================================
//====================================================================================


MinCardCS::MinCardCS(double uss_, double dss_, int id_, int rhs, int ss_size, const int * SS_arcs, int id_owner){
    next=0;
    prev =0;
    id = id_;
    card = rhs;
    csid = id_owner;
    arcs = SS_arcs;
    size = ss_size;
    uss = uss_;
    dss = dss_;
    bin=0;
    actv = 0;
}

//----------------------------------------------------------------------------------
//  card methods
//----------------------------------------------------------------------------------

void MinCardCS::mark(int binpos, int arcpos){
    setBit(bin,binpos,arcpos);
}

//====================================================================================
//====================================================================================
// MinCardCSCollection methods
//====================================================================================
//====================================================================================


void
MinCardCSCollection::initialize(int M){
    int sizeOfInt=8*sizeof(unsigned int);
    sizeOfContainer = (M/sizeOfInt)+1;
    sizeOfMap = M;
    sizeOfCollection = 0;
    map =new Pair2[M];
    empty=true;
    for(int i=0;i<M;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
        //std::cout<<i<<" "<<map[i].fst<<" : "<<map[i].snd<<std::endl;
    }
    sizeOfCollection = 0;
    empty=true;
    end = begin =0;
}

//----------------------------------------------------------------------------------


MinCardCSCollection::~MinCardCSCollection(){
    for(int i=0;i<sizeOfCollection;++i){
        MinCardCS* next = begin->next;
        delete begin;
        begin = next;
    }
    delete [] map;
}

//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------

void
MinCardCSCollection::clearbin(MinCardCS * card){
    for(int i=0;i<sizeOfContainer;++i){
        card->bin[i]=0;
    }
    card->actv =0;
}

//----------------------------------------------------------------------------------

void
MinCardCSCollection::mark(int arc, MinCardCS * card){
    card->mark(map[arc].fst, map[arc].snd);
}

//----------------------------------------------------------------------------------

int
MinCardCSCollection::equalbin(const MinCardCS * card1, const MinCardCS * card2){
    unsigned int eqlvl=0;
    for(int i=0;i<sizeOfContainer;++i){
        eqlvl = countSetBits(card1->bin[i] & card2->bin[i]);
    }
    return eqlvl;
}

//----------------------------------------------------------------------------------
//  VI methods
//----------------------------------------------------------------------------------

int
MinCardCSCollection::addMinCardCS(MinCardCS * newvi){
        
    if(collected(newvi)){ delete newvi; return 0;}
    else{
        newvi->bin = new unsigned int [sizeOfContainer];
        for(int i=0;i<sizeOfContainer;++i){ newvi->bin[i]=0;}

        ++sizeOfCollection;
        if(empty){
            end = begin = newvi;
            empty=false;
        }else{
            newvi->prev = end;
            end->next= newvi;
            end = newvi;
        }
        return 1;
    }
    return 0;
}


//------------------------------------------------------------------------------------

int
MinCardCSCollection::collected(MinCardCS * tryvi){
    MinCardCS* vi = begin;
    for(int i=0;i<sizeOfCollection;++i){
        if( vi->csid == tryvi->csid ){
            return 1;
        }
        vi = vi->next;
    }
    return 0;
}
