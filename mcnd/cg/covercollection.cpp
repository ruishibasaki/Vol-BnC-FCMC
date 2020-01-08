//
//  covervi.cpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#include "covercollection.hpp"


//====================================================================================
//====================================================================================
// Class CoverCollection methods
//====================================================================================
//====================================================================================

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  Construction/Destruction methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------


void
CoverCollection::initialize(int M){
    sizeOfMap = M;
    int sizeOfInt=8*sizeof(unsigned int);
    sizeOfIdSeq = (M/sizeOfInt)+1;
    
    sizeOfCollection = 0;
    map =new Pair2[M];
    empty=true;
    for(int i=0;i<M;++i){
        map[i].fst = i/sizeOfInt;
        map[i].snd = i%sizeOfInt;
        //std::cout<<i<<" "<<map[i].fst<<" : "<<map[i].snd<<std::endl;
        
    }
    newly_added=0;
    discarted=0;
    begin_actv = end_actv = end = begin =0;
}

//----------------------------------------------------------------------------------


CoverCollection::~CoverCollection(){
    for(int i=0;i<sizeOfCollection;++i){
        Cover* next = begin->next;
        delete begin;
        begin = next;
    }
    if(map) delete [] map;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  setter getter methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

Cover *
CoverCollection::createNewCover(const std::deque<Pair2>& c, double mu,  int id_vi,  int id_owner_, int serial_num_){
    int arc;
    int csize =(int) c.size();
    Cover * newC = new Cover(sizeOfIdSeq, csize,  id_vi, id_owner_, serial_num_);
    
    for(int n=csize;n--;){
        arc = c[n].fst;
        newC->addArc(map[arc].fst, map[arc].snd);
        newC->C[n]=arc;
        //std::cout<<"c: "<<arc<<std::endl;
    }
    newC->rhs = 1.0;
    return newC;
}

//----------------------------------------------------------------------------------

int 
CoverCollection::addCover(Cover * tryC, const double * xystar){
    int arc, ret;
    if(tryC->Lftd){
        CoverL* liftd = tryC->Lftd;
        for(int n=liftd->size;n--;){
            arc = liftd->SSmC[n];
            if(liftd->gamma[n]>0)
                tryC->addArc(map[arc].fst, map[arc].snd);
        }
    }
    ret = collected(tryC);
    if(ret==2){ ++discarted; delete tryC;return 0;}
    else if(ret==0){
        //ret = check_maximal(tryC , xystar);
        //if(ret) return 0;
        ++sizeOfCollection;
        if(empty){
            begin_actv = end_actv = end = begin = tryC;
            empty=false;
        }else{
            tryC->prev = end;
            end->next= tryC;
            end = tryC;
            end_actv->next= tryC;
            end_actv = tryC;
        }
        return 1;
    }else return 0;

    
}

//----------------------------------------------------------------------------------

void
CoverCollection::replace(Cover * out, Cover * in){
    in->next = out->next;
    in->prev = out->prev;
    in->id_vi = out->id_vi;
    if(out == begin)
        begin = in;
    else out->prev->next = in;
    if(out == end)
        end = in;
    else out->next->prev = in;
    
    delete out;
    
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  VI methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int
CoverCollection::collected(Cover * tryC){
    int ret=0;
    int vi1sz=0;
    int vi2sz=0;
    bool equal;
    Cover* C = begin;
    for(int i=0;i<sizeOfCollection;++i){
        vi1sz = tryC->get_total_sz();
        vi2sz = C->get_total_sz();            
        if( vi1sz == vi2sz ){
            equal = true;
            for(int id=0;id<sizeOfIdSeq;++id){
                if(tryC->id_seq[id]!=C->id_seq[id]){
                    equal=false;
                    break;
                }
            }
            
            if(equal){
                ret = compScalar(tryC, C);
                if(ret==0) return 0;
                else if(ret==1){
                    replace(C, tryC);
                    return 1;
                }
                return 2;
            }
        }
        C = C->next;
    }
    return 0;
}

//----------------------------------------------------------------------------------

int
CoverCollection::compScalar(Cover * c1, Cover * c2){
    int arc;
    int sz = c1->maxsize;
    double div;
    double factor;
    double cfactor;
    double gamma1, gamma2;
    bool equal = true;
    bool v1dominate = true;
    bool v2dominate = true;
    std::vector<PairF> mapset;
    
    gamma1 = c1->rhs;
    gamma2 = c2->rhs;
    if(c1->hasLftd) gamma1 += c1->Lftd->ttgamma1;
    if(c2->hasLftd) gamma2 += c2->Lftd->ttgamma1;
    factor = (gamma1) / double(gamma2);
    cfactor = 1/double(factor);
    
    mapCover(c1, c2, mapset);
    for(int n=0;n<sz;++n){
        arc = c1->owner[n];
        
        if( !getScalar(arc, mapset, gamma1, gamma2)) continue;
        
        div  = gamma1/(gamma2);
        if( (gamma2 != gamma1) && (div) != factor){
            equal = false;
        }
        
        if(div<1.0 && factor>1.0) v2dominate = false;
        if(div>1.0 && factor<1.0) v1dominate = false;
        //check dominance of vi2
        if(factor>=1.0){
            //need to drecrease rhs1 the minimum of "1/factor", if for some item vi1*(1/factor) < vi2 then not dominant
            if((gamma2) > (gamma1/factor)) v2dominate=false;
        }else{
            //have a max of "1/factor" to increase vi1 if vi1 < vi2. if the needed incrase is greater than max then not dominant
            if((div<1.0) && (1.0/div > cfactor) ) v2dominate=false;
        }
        //check dominance of vi1
        if(factor<=1.0){
            //need to drecrease rhs2 the minimum of "factor", if for some item factor*vi2 < vi1 then not dominant
            if((gamma1) > (gamma2*factor)) v1dominate=false;
        }else{
            //have a max of "factor" to increase vi2 if vi2 < vi1. if the needed increase is greater than max then not dominant
            if((div>1.0) && (div > factor)) v1dominate=false;
        }
    }
    
    mapset.clear();
    if(equal){
        //std::cout<<"equal"<<std::endl;
        return -1;
    }
    
    if(v2dominate){
        c1->print();
        c2->print();
        std::cout<<"v2dominate"<<std::endl;
        if(v1dominate) std::cout<<"!!!!! PROBLEM::compLiftPart domdom!!!!!"<<std::endl;
        return 2;
    }else if(v1dominate){
        c1->print();
        c2->print();
        std::cout<<"v1dominate"<<std::endl;
        return 1;
    }
    return 0;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//   modifying methods
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

Cover*
CoverCollection::move_to_end(Cover * trgt){
	Cover * ret;
    if(trgt == end ){
        return 0;
    }else if(trgt == begin){
        trgt->next->prev = 0;
        begin = trgt->next;
    }else{
     	trgt->next->prev = trgt->prev;
     	trgt->prev->next = trgt->next;
    }
    ret = trgt->next;
    trgt->next=0;
    end->next = trgt;
    trgt->prev=end;
    end = trgt;
    return ret; 
}

//------------------------------------------------------------------------------------

Cover*
CoverCollection::swap_to_end(Cover * trgt){
    if(trgt == end){return 0; }
    Cover * aux; Cover * ret;
    int id_aux;
    
    id_aux= end->id_vi;
    end->id_vi = trgt->id_vi;
    trgt->id_vi = id_aux;
    
    aux = end->prev;
    end->prev = trgt->prev;
    trgt->prev = aux;
    if(end->prev){ end->prev->next = end; }
    if(trgt->prev){ trgt->prev->next = trgt; }
    end->next = trgt->next;
    trgt->next = 0;
    if(end->next){ end->next->prev = end; }
    if(trgt == begin){ begin = end; }
    ret = end;
    end = trgt;
    return ret;
}

//----------------------------------------------------------------------------------

int
CoverCollection::swap_toend_destruct(Cover * out, Cover *& ret, bool destruct){
    Cover* new_end=0;
    int id;
    //std::cout<<"out: "<<out->id_vi<<" end: "<<end->id_vi<<std::endl;
    if(sizeOfCollection == 1){
        ret = 0;
        id =-1;
        begin =0;
        empty=true;
        //std::cout<<" EMPTY!!!"<<std::endl;
    }else{
        id = end->id_vi;
        ret = swap_to_end(out);
        if(ret==0) id=-1;
        new_end = end->prev;
        new_end->next =0;
        //std::cout<<" out_in_end: "<<end->id_vi<<" ret: "<<ret->id_vi<<" last_end_id: "<<id<<std::endl;
    }
    --sizeOfCollection;
    ++discarted;
    
    if(destruct) delete end;
    end = new_end;
    return id;
}


//----------------------------------------------------------------------------------


int
CoverCollection::removeCover(int lim, int * actvS, int & actvSSz, double * pstarv, double * dstaru, double * dualu){
    int idx_end, idx_out, id_end;
    Cover * ret=0;
    
    //std::cout<<"remove check "<<std::endl;
    //for(int h=2000;h<actvSSz;h++)std::cout<<"d: "<<h<<" pv: "<<pstarv[actvS[h]]<<" d*: "<<dstaru[actvS[h]]<<" d: "<<dualu[actvS[h]]<<std::endl;
    for(Cover *vi = begin; vi != 0;){
        idx_out = actvS[vi->id_vi];
        //std::cout<<"vi->id_vi: "<<vi->id_vi<<" "<<vi<<" "<<end<<" "<<vi->n_noviol<<" "<<idx_out<<std::endl;
        if(vi->n_nviol >= lim /*&& dualu[idx_out]<=0 */&& dstaru[idx_out]<=0){
            //std::cout<<"remove vi: "<<vi->id_vi<<" "<<idx_out<<std::endl;;
            id_end = swap_toend_destruct(vi, ret, true);
            if( id_end <0 ){
                actvS[actvSSz-1]=-1;
                dstaru[actvSSz-1] = 0;
                dualu[actvSSz-1] = 0;
                pstarv[actvSSz-1] = 0;
                //std::cout<<" id_end <0, idx: "<<idx_out<<" new end: "<<end->id_vi<<" sz: "<<sizeOfCollection<<" id_end: "<<id_end<<" now: "<<end->id_vi<<std::endl;
                
                --actvSSz;
                return 0;
            }
            
            idx_end = actvS[id_end];
            actvS[id_end]=-1;
            pstarv[idx_out] = pstarv[idx_end];
            dstaru[idx_out] = dstaru[idx_end];
            dualu[idx_out] = dualu[idx_end];
            dstaru[idx_end] = 0;
            dualu[idx_end] = 0;
            pstarv[idx_end] = 0;
            --actvSSz;
            //std::cout<<" size: "<<actvSSz<<std::endl;
            //std::cout<<" idx: "<<idx_out<<" new end: "<<end->id_vi<<" sz: "<<sizeOfCollection<<" idx_end: "<<idx_end<<" id_end: "<<id_end<<" now: "<<ret->id_vi<<std::endl;
            
            vi = ret;
        }else if( vi == end ){ vi=0; //std::cout<<"end"<<std::endl;
        }else vi = vi->next;
    }
    /*RFCI * t = begin;
     std::cout<<"idseq: ";
     for(int sz = sizeOfCollection;sz--;){
     std::cout<<t->id_vi<<" ";
     t = t->next;
     }
     std::cout<<std::endl;*/
    //std::cout<<"ok"<<std::endl;
    //for(int h=2000;h<actvSSz;h++)std::cout<<"d: "<<h<<" pv: "<<pstarv[actvS[h]]<<" d*: "<<dstaru[actvS[h]]<<" d: "<<dualu[actvS[h]]<<std::endl;
    return 0;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

Cover*
CoverCollection::operator[](int n){
    if(n>=sizeOfCollection){std::cout<<"PROBLEM::CoverCollection::operator[]   n>=sizeOfCollection "<<std::endl; return 0;}
    Cover* C = begin_actv;
    for(int i=1;i<=n;++i)
        C = C->next;
    return C;
}

//----------------------------------------------------------------------------------

const void
CoverCollection::advance(Cover*& C, int n){
    if(n>sizeOfCollection) return;
    C = begin_actv;
    for(int i=1;i<=n;++i)
        C = C->next;
}

//----------------------------------------------------------------------------------

const void
CoverCollection::print(){
    Cover* c = begin;
    for(int i=0;i<sizeOfCollection;++i){
        std::cout<<"cover "<<i<<" "<<std::endl;
        for(int id=0;id<sizeOfIdSeq;++id){
            std::cout<<c->id_seq[id]<<" ";
        }
        std::cout<<std::endl;
        for(int id=0;id<c->size;++id){
            std::cout<<c->C[id]<<" ";
        }
        std::cout<<std::endl;
        c = c->next;
    }
}

//----------------------------------------------------------------------------------

void
CoverCollection::mapCover(Cover * c1, Cover * c2, std::vector<PairF> & mapset){
    mapset.assign(sizeOfMap, PairF());
    for(int n=c1->size;n--;){
        mapset[c1->C[n]].fst = 1.0;
    }
    if(c1->hasLftd){
        for(int n=c1->Lftd->size;n--;){
            mapset[c1->Lftd->SSmC[n]].fst = c1->Lftd->gamma[n];
        }
    }
    for(int n=c2->size;n--;){
        mapset[c2->C[n]].snd = 1.0;
    }
    if(c2->hasLftd){
        for(int n=c2->Lftd->size;n--;){
            mapset[c2->Lftd->SSmC[n]].snd = c2->Lftd->gamma[n];
        }
    }
}

//----------------------------------------------------------------------------------

double
CoverCollection::cover_hasArc(const Cover * cover, int arc){
    if(!cover->hasArc(map[arc].fst, map[arc].snd)) return 0;
    if(cover->hasLftd){
        for(int n=cover->Lftd->size;n--;)
            if(cover->Lftd->SSmC[n]==arc) return cover->Lftd->gamma[n];
    }
    for(int n=cover->size;n--;)
        if(cover->C[n]==arc) return 1.0;
    
    return 0;
}


//====================================================================================
//====================================================================================
// Class Cover methods
//====================================================================================
//====================================================================================


int
Cover::at(int pos) const{
    if(pos >= get_total_sz()) return -1;
    
    if(pos >= size){
        if(hasLftd){
            return Lftd->SSmC[pos-size];
        }else return -1;
    }else return C[pos];
}

//----------------------------------------------------------------------------------

double
Cover::gamma_at(int pos) const{
    if(pos >= get_total_sz()) return -1;
    
    if(pos >= size){
        if(hasLftd){
            return Lftd->gamma[pos-size];
        }else return -1;
    }else return 1.0;
}

//----------------------------------------------------------------------------------


void
Cover::get_total_sz_rhs(int & sz, double &rhs) const{
    sz = size;
    rhs = rhs- rhs_dimsh;
    if(hasLftd){
        sz += Lftd->size;
        rhs += Lftd->ttgamma1;
    }
}

//----------------------------------------------------------------------------------

double
Cover::get_total_rhs() const{
    double rhs_ = rhs - rhs_dimsh;
    if(hasLftd){
        rhs_ += Lftd->ttgamma1;
    }
    
    return rhs_;
}

//----------------------------------------------------------------------------------


int
Cover::get_total_sz()const{
    int sz = size;
    if(hasLftd){
        sz += Lftd->size;
    }
    return sz;
}

//----------------------------------------------------------------------------------

void Cover::addArc(int iset, int arc){
    setBit(id_seq,iset,arc);
}

//----------------------------------------------------------------------------------


void Cover::removeArc(int iset, int arc){
    clearBit(id_seq,iset,arc);
}

//----------------------------------------------------------------------------------

bool Cover::hasArc(int iset, int arc) const{
    if(testBit(id_seq,iset,arc)) return true;
    else return false;
}

//----------------------------------------------------------------------------------

void Cover::addLftd(CoverL * l){
    if(l==0) return;
    
    if(l->size==0){
        delete l; Lftd=0; l=0;
        return;
    }
    
    Lftd=l;
    Lftd->owner = this;
    hasLftd=true;
}

//----------------------------------------------------------------------------------

void Cover::print(){
    double rhs_ = rhs;
    std::cout<<"C: ";
    for(int i=size;i--;)
        std::cout<<"("<<C[i]<<", mtl: 1) ";
    if(hasLftd){
        rhs_ += Lftd->ttgamma1;
        for(int i=Lftd->size;i--;)
            std::cout<<"("<<Lftd->SSmC[i]<<", mtl:"<<Lftd->gamma[i]<<") ";
    }
    std::cout<<" rhs: "<<rhs_<<std::endl;
}

//====================================================================================
//====================================================================================
// Class CoverVI methods
//====================================================================================
//====================================================================================


void
CoverL::addvar(int arc, double gam, bool lift_down){
    SSmC[size]=arc; gamma[size]=gam; ++size;
    if(lift_down)ttgamma1+=gam;
}
