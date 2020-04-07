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

    discarted=0;
    end = begin =0;
}

//----------------------------------------------------------------------------------

Cover *
CoverCollection::createNewCover(const std::deque<Pair2>& c,  int id_vi, int serial_num_){
    int arc;
    int csize =(int) c.size();
    Cover * newC = new Cover(sizeOfIdSeq, csize,  id_vi, serial_num_);
    
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


CoverCollection::~CoverCollection(){
     
    if(map) delete [] map;
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
                    //replace(C, tryC);
                    C->prgbl=true;
                    //abort();
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
    if(sz==0) return 0;
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
        //c1->print();
        //c2->print();
        //std::cout<<"v2dominate"<<std::endl;
        if(v1dominate){ std::cout<<"!!!!! PROBLEM::compLiftPart domdom!!!!!"<<std::endl; abort();}
        return 2;
    }else if(v1dominate){
        //c1->print();
        //c2->print();
        //std::cout<<"v1dominate"<<std::endl;
        return 1;
    }
    return 0;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  insert/del methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

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
    if(ret==2){++discarted; delete tryC;return 0;}
    else if(ret==0){
    	//if(xystar){
    		//ret = check_maximal(tryC , xystar);
        	//if(ret) return 0;
    	//}
        insert_end(tryC);
        return 1;
    }else return 0;

    
}

//----------------------------------------------------------------------------------

void 
CoverCollection::insert_end(Cover * tryC){
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
CoverCollection::insert_front(Cover * tryC){
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

//----------------------------------------------------------------------------------

void 
CoverCollection::pop_back_nodel(){
	if(sizeOfCollection==1){
		begin = 0;
		empty = true;
	}
	end = end->prev;
	if(end) end->next =0;
	--sizeOfCollection;

}

//----------------------------------------------------------------------------------

Cover * 
CoverCollection::remove_nodel(Cover * trgt){
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
	Cover * ret  = trgt->next;
	if(trgt->next)trgt->next->prev = trgt->prev;
	if(trgt->prev)trgt->prev->next = trgt->next;
	trgt->next =0;
	trgt->prev =0;
	return ret;
}

//----------------------------------------------------------------------------------


int
CoverCollection::desactvCover(int lim, int * actvS, int & actvSSz, int& num_actv, double * pstarv, double * dstaru, double * dualu){
    int idx_end, idx_out;
    int curr_num = num_actv;
    Cover * ret=0;
    
    //std::cout<<"remove check "<<std::endl;
    //for(int h=2000;h<actvSSz;h++)std::cout<<"d: "<<h<<" pv: "<<pstarv[actvS[h]]<<" d*: "<<dstaru[actvS[h]]<<" d: "<<dualu[actvS[h]]<<std::endl;
   	Cover *vi = begin;
   	Cover * last_actv;
   	//std::cout<<"reposition_covers "<<cover_manager->covers.sizeOfCollection <<" "<< cover_manager->num_actv<<std::endl;
    for(int i=curr_num;i--;){
        idx_out = actvS[vi->id_vi];
        //std::cout<<"try vi: "<<vi->serial_nmbr<<std::endl;
        if((vi->n_nviol >= lim /*&& dualu[idx_out]<=0 */&& dstaru[idx_out]<=0)|| vi->prgbl){
            //std::cout<<"remove vi: "<<vi->id_vi<<" "<<" srial: "<<vi->serial_nmbr<<" id: "
            //<<idx_out<<" purgbl: "<<vi->prgbl<<" nviol: "<<vi->n_nviol<<" d* "<<dstaru[idx_out]<<std::endl;;
            
            last_actv = begin ;
			advance(last_actv, num_actv-1);
			//std::cout<<"last actv: "<<last_actv->id_vi<<" serial: "<<last_actv->serial_nmbr<<std::endl;
			idx_end = actvS[last_actv->id_vi];
            if(swap(vi,last_actv)){
            	pstarv[idx_out] = pstarv[idx_end];
				dstaru[idx_out] = dstaru[idx_end];
				dualu[idx_out] = dualu[idx_end];
				dstaru[idx_end] = 0;
				dualu[idx_end] = 0;
				pstarv[idx_end] = 0;
				actvS[last_actv->id_vi] = idx_out;
            }
            actvS[vi->id_vi] = -1;
            --actvSSz;
            --num_actv;
            //std::cout<<" size: "<<actvSSz<<std::endl;
            //std::cout<<" idx: "<<idx_out<<" new end: "<<end->id_vi<<" sz: "<<sizeOfCollection<<" idx_end: "<<idx_end<<" id_end: "<<id_end<<" now: "<<ret->id_vi<<std::endl;
            
            vi = last_actv;
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

bool
CoverCollection::swap(Cover * c1, Cover * c2){
	Cover* aux;
	if(c1==c2) return false;
	if(c1 == begin) begin =c2;
	if(c2 == end) end = c1;

	if(c1->prev){c1->prev->next = c2;	}
	if(c2->next){c2->next->prev = c1; 	}
	
	if(c1->next == c2){
		c1->next=c2->next;
		c2->next=c1;
		c2->prev=c1->prev;
		c1->prev=c2;
		return true;
	}
	
	if(c1->next){c1->next->prev = c2;}
	if(c2->prev){c2->prev->next = c1;}
	
	aux = c1->next;
	c1->next = c2->next;
	c2->next = aux;
	
	aux = c1->prev;
	c1->prev = c2->prev;
	c2->prev = aux;
 	return true;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

Cover*
CoverCollection::operator[](int n){
    if(n>=sizeOfCollection){std::cout<<"PROBLEM::CoverCollection::operator[]   n>=sizeOfCollection "<<std::endl; return 0;}
    Cover* C = begin;
    for(int i=1;i<=n;++i)
        C = C->next;
    return C;
}

//----------------------------------------------------------------------------------

void
CoverCollection::advance(Cover*& C, int n){
    if(n>sizeOfCollection) return;
    C = begin;
    for(int i=1;i<=n;++i)
        C = C->next;
}

//----------------------------------------------------------------------------------

void
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

//----------------------------------------------------------------------------------

void 
CoverCollection::map_collection(std::map<int, int>& mapd){
	int sz = sizeOfCollection;
	Cover *vi = begin;
	for(;sz--;){
		mapd.insert(std::pair<int,int>(vi->serial_nmbr,vi->id_vi));
		vi = vi->next;
	}
}

//====================================================================================
//====================================================================================
// Class Cover methods
//====================================================================================
//====================================================================================

bool 
Cover::check_updt_Viol(const double *y, bool & infeas){
	double sum=0;
    double rhs_= get_rhs();
    int sz = get_total_sz();
    int arc;
     
    rhs_dimsh=0;
    for(int a=0;a<sz;++a){
     	arc = at(a);
		if(y[arc] ==1){
        	rhs_dimsh+= gamma_at(a);
		}else if(y[arc]==-1){
			sum+= gamma_at(a);
 		}
    }
	rhs_ -= rhs_dimsh;
 	infeas=false;
	if(sum<rhs_){ infeas= true; return false;} //abort;
	if(rhs_ <= 0){return false;}
	 
    return true;
}

//----------------------------------------------------------------------------------

double 
Cover::viol(const double *y)const {
	double sum=0;
    double comp= get_rhs();
    int sz = get_total_sz();
    for(int a=0;a<sz;++a){
        sum+= gamma_at(a)*y[at(a)];
        if(sum>=comp){return 0;}
    }
    return comp-sum;
}

//----------------------------------------------------------------------------------

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

double
Cover::get_rhs() const{
    double rhs_ = rhs;
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
