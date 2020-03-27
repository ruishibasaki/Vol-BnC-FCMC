//
//  localcutcollection.cpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#include "localcutcollection.hpp"


void
LocalCutCollection::initialize(int M){
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

    discarted=0;
    end = begin =0;
}

//----------------------------------------------------------------------------------

LocalCut *
LocalCutCollection::createNewLocalCut(const std::vector<int>& c, int id_vi, int serial_num_, int sense_, double rhs_){
    int arc;
    int csize =(int) c.size();
    LocalCut * newC = new LocalCut(sizeOfIdSeq, csize,  id_vi, serial_num_, sense_);
    
    for(int n=csize;n--;){
        arc = c[n];
        newC->addArc(map[arc].fst, map[arc].snd);
        newC->vars[n]=arc;
        //std::cout<<"c: "<<arc<<std::endl;
    }
    newC->rhs = rhs_;
    return newC;
}

//----------------------------------------------------------------------------------

LocalCutCollection::~LocalCutCollection(){
    for(int i=0;i<sizeOfCollection;++i){
        LocalCut* next = begin->next;
        delete begin;
        begin = next;
    }
    if(map) delete [] map;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  VI methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int
LocalCutCollection::collected(LocalCut * tryC){
    int ret=0;
    int vi1sz=0;
    int vi2sz=0;
    bool equal;
    LocalCut* C = begin;
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
            
            if(equal){ return 2; }
        }
        C = C->next;
    }
    return 0;
}


//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  insert/del methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int 
LocalCutCollection::addLocalCut(LocalCut * tryC){
    int arc, ret;
    ret = collected(tryC);
    if(ret==2){++discarted; delete tryC;return 0;}
    else if(ret==0){
        //insert_end(tryC);
        std::cout<<"local add: "; tryC->print();
        //return 1;
    } 
    return 0;
}

//----------------------------------------------------------------------------------

void 
LocalCutCollection::insert_end(LocalCut * tryC){
	++sizeOfCollection;
	if(empty){
		end = begin = tryC;
		empty=false;
	}else{
		tryC->prev = end;
		end->next= tryC;
		end = tryC;
	}
}

//----------------------------------------------------------------------------------

void 
LocalCutCollection::insert_front(LocalCut * tryC){
	++sizeOfCollection;
	if(empty){
		end = begin = tryC;
		empty=false;
	}else{
		tryC->next = begin;
		begin->prev= tryC;
		begin = tryC;
	}
}

//----------------------------------------------------------------------------------

void
LocalCutCollection::replace(LocalCut * out, LocalCut * in){
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

//----------------------------------------------------------------------------------

void 
LocalCutCollection::pop_back_nodel(){
	if(sizeOfCollection==1){
		begin = 0;
		empty = true;
	}
	end = end->prev;
	if(end) end->next =0;
	--sizeOfCollection;

}

//----------------------------------------------------------------------------------

LocalCut * 
LocalCutCollection::remove_nodel(LocalCut * trgt){
	if(trgt == end){
		pop_back_nodel();
		return 0;
	}
	 
	--sizeOfCollection;
	if(trgt == begin){
		begin = trgt->next;
		if(begin) begin->prev =0;
		return 0;
	}
	LocalCut * ret  = trgt->next;
	if(trgt->next)trgt->next->prev = trgt->prev;
	if(trgt->prev)trgt->prev->next = trgt->next;
	trgt->next =0;
	trgt->prev =0;
	return ret;
}

//----------------------------------------------------------------------------------


int
LocalCutCollection::removeLocalCut(int lim, int * actvS, int & actvSSz, double * pstarv, double * dstaru, double * dualu){
    int idx_end, idx_out, id_end;
    LocalCut * ret=0;
    
    //std::cout<<"remove check "<<std::endl;
    //for(int h=2000;h<actvSSz;h++)std::cout<<"d: "<<h<<" pv: "<<pstarv[actvS[h]]<<" d*: "<<dstaru[actvS[h]]<<" d: "<<dualu[actvS[h]]<<std::endl;
    for(LocalCut *vi = begin; vi != 0;){
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
//------------------------------------------------------------------------------------
//   modifying methods
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

LocalCut*
LocalCutCollection::move_to_end(LocalCut * trgt){
	LocalCut * ret;
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

LocalCut*
LocalCutCollection::swap_to_end(LocalCut * trgt){
    if(trgt == end){return 0; }
    LocalCut * aux; LocalCut * ret;
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
LocalCutCollection::swap_toend_destruct(LocalCut * out, LocalCut *& ret, bool destruct){
    LocalCut* new_end=0;
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

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

LocalCut*
LocalCutCollection::operator[](int n){
    if(n>=sizeOfCollection){std::cout<<"PROBLEM::LocalCutCollection::operator[]   n>=sizeOfCollection "<<std::endl; return 0;}
    LocalCut* C = begin;
    for(int i=1;i<=n;++i)
        C = C->next;
    return C;
}

//----------------------------------------------------------------------------------

void
LocalCutCollection::advance(LocalCut*& C, int n){
    if(n>sizeOfCollection) return;
    C = begin;
    for(int i=1;i<=n;++i)
        C = C->next;
}

//----------------------------------------------------------------------------------

void
LocalCutCollection::print(){
    LocalCut* c = begin;
    for(int i=0;i<sizeOfCollection;++i){
        std::cout<<"LocalCut "<<i<<" "<<std::endl;
        for(int id=0;id<sizeOfIdSeq;++id){
            std::cout<<c->id_seq[id]<<" ";
        }
        std::cout<<std::endl;
        for(int id=0;id<c->size;++id){
            std::cout<<c->vars[id]<<" ";
        }
        std::cout<<std::endl;
        c = c->next;
    }
}

//----------------------------------------------------------------------------------

double
LocalCutCollection::LocalCut_hasArc(const LocalCut * LocalCut, int arc){
    if(!LocalCut->hasArc(map[arc].fst, map[arc].snd)) return 0;
    return 1.0;
}

//----------------------------------------------------------------------------------

void 
LocalCutCollection::map_collection(std::map<int, int>& mapd){
	int sz = sizeOfCollection;
	LocalCut *vi = begin;
	for(;sz--;){
		mapd.insert(std::pair<int,int>(vi->serial_nmbr,vi->id_vi));
		vi = vi->next;
	}
}

//====================================================================================
//====================================================================================
// Class LocalCut methods
//====================================================================================
//====================================================================================

bool 
LocalCut::check_updt_Viol(const double *y) {
	double sum=0;
    double comp= get_rhs();
    rhs_dimsh=0;
    for(int a=0;a<size;++a){
        sum+= y[vars[a]];
        if(sum>=comp){return false;}
    }
    rhs_dimsh = sum;
    return true;
}

//----------------------------------------------------------------------------------


double 
LocalCut::viol(const double *y)const {
	double sum=0;
    double comp= get_rhs();
     for(int a=0;a<size;++a){
        sum+= y[vars[a]];
        if(sum>=comp){return 0;}
    }
    return comp-sum;
}

//----------------------------------------------------------------------------------

int
LocalCut::at(int pos) const{
    if(pos >= get_total_sz()) return -1;
     return vars[pos];
}

//----------------------------------------------------------------------------------

double
LocalCut::gamma_at(int pos) const{
    if(pos >= size) return -1;
    return 1.0;
}

//----------------------------------------------------------------------------------


void
LocalCut::get_total_sz_rhs(int & sz, double &rhs) const{
    sz = size;
    rhs = rhs- rhs_dimsh;
}

//----------------------------------------------------------------------------------

double
LocalCut::get_total_rhs() const{ return  rhs - rhs_dimsh;}

//----------------------------------------------------------------------------------

double
LocalCut::get_rhs() const{  return rhs; }

//----------------------------------------------------------------------------------

int
LocalCut::get_total_sz()const{    return size;}

//----------------------------------------------------------------------------------

void LocalCut::addArc(int iset, int arc){
    setBit(id_seq,iset,arc);
}

//----------------------------------------------------------------------------------


void LocalCut::removeArc(int iset, int arc){
    clearBit(id_seq,iset,arc);
}

//----------------------------------------------------------------------------------

bool LocalCut::hasArc(int iset, int arc) const{
    if(testBit(id_seq,iset,arc)) return true;
    else return false;
}

//----------------------------------------------------------------------------------

void LocalCut::print(){
    double rhs_ = rhs;
    std::cout<<"T"<<sense<<": ";
    for(int i=size;i--;)
        std::cout<<"("<<vars[i]<<", mtl: 1) ";
    
    std::cout<<" rhs: "<<rhs_<<std::endl;
}

