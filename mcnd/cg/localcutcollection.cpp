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
LocalCutCollection::createNewLocalCut(int sz,  int* vars_, int id_vi, int serial_num_){
    int arc;
    LocalCut * newC = new LocalCut(sizeOfIdSeq, sz,  id_vi, serial_num_, 1, 0);
    newC->vars = vars_;  
    
	for(int n=sz;n--;){
		arc = vars_[n];
		newC->addArc(map[arc].fst, map[arc].snd);	 
	}
    newC->rhs = 1.0;
    vars_=0;
    return newC;
}

//----------------------------------------------------------------------------------

LocalCut *
LocalCutCollection::createNewLocalCut(const std::vector<int>& c, int id_vi, int serial_num_, int sense_, double rhs_){
    int arc;
    int csize =(int) c.size();
    LocalCut * newC = new LocalCut(sizeOfIdSeq, csize,  id_vi, serial_num_, sense_, 1);
    newC->vars = new int[csize];
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


LocalCut * 
LocalCutCollection::createNewLocalCut(int sz,  int* vars_, double * coef_, double rhs_, int id_vi, 
    								int serial_num_){
    int arc;
    LocalCut * newC = new LocalCut(sizeOfIdSeq, sz,  id_vi, serial_num_, 1, 2);
    newC->vars = vars_;  
    newC->coef  = coef_;
    newC->rhs = rhs_;
	for(int n=sz;n--;){
		if(coef_[n]==-1) continue;
		arc = vars_[n];
		newC->addArc(map[arc].fst, map[arc].snd);	 
	}
    vars_=0;
    coef_=0;
    return newC;						
    								
}
    								
//----------------------------------------------------------------------------------

LocalCutCollection::~LocalCutCollection(){
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
    //std::cout<<"coll: "<<sizeOfCollection<<std::endl;
    for(int i=0;i<sizeOfCollection;++i){
        vi1sz = tryC->size;
        vi2sz = C->size;            
        if( vi1sz == vi2sz && tryC->sense == C->sense && tryC->rhs == C->rhs){
            equal = true;
            for(int id=0;id<sizeOfIdSeq;++id){
                if(tryC->id_seq[id]!=C->id_seq[id]){
                    equal=false;
                    break;
                }
            }
            
            if(equal){  return 1; }
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
    if(ret==1){  delete tryC; return 0;}
    else if(ret==0){
        //std::cout<<"local add: "; tryC->print();
        insert_end(tryC);
        return 1;
    }
    delete tryC; 
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

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//   modifying methods
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int
LocalCutCollection::desactvLocalc(int lim, int * actvS, int & actvSSz, int& num_actv, double * pstarv, double * dstaru, double * dualu){
    int idx_end, idx_out;
    int curr_num = num_actv;
    LocalCut * ret=0;
    
    //std::cout<<"remove check "<<std::endl;
   	LocalCut *vi = begin;
   	LocalCut * last_actv;
    for(int i=curr_num;i--;){
        idx_out = actvS[vi->id_vi];
        //std::cout<<"try vi: "<<vi->serial_nmbr<<std::endl;
        if((vi->n_nviol >= lim /*&& dualu[idx_out]<=0 */&& dstaru[idx_out]<=0)){
            //std::cout<<"remove vi: "<<vi->id_vi<<" "<<" srial: "<<vi->serial_nmbr<<std::endl;
            //<<idx_out<<" nviol: "<<vi->n_nviol<<" d* "<<dstaru[idx_out]<<std::endl;;
            
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
            
            vi = last_actv;
        }else vi = vi->next;
    }
   
    return 0;
}

//----------------------------------------------------------------------------------

bool
LocalCutCollection::swap(LocalCut * c1, LocalCut * c2){
	LocalCut* aux;
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
LocalCutCollection::LocalCut_hasArc(const LocalCut * localCut, int arc){
    if(!localCut->hasArc(map[arc].fst, map[arc].snd)) return 0;
    if(localCut->coef){
		for(int n=localCut->size;n--;)
			if(localCut->vars[n]==arc) return localCut->coef[n];
		return 0;
	} return 1.0;
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
LocalCut::check_updt_Viol2(const double *y, bool & infeas) {
	int arc;
    int ntofx=0;
    double coef_;
	double viol=true;
	rhs_dimsh=0;

	for(int a=0;a<size;++a){
		arc = vars[a];
		coef_ = coef_at(a);
		if(y[arc]==1){
			rhs_dimsh+= coef[a];
			if(coef[a]==1){
				viol= false;
			}
		}else if(y[arc]==-1){
			++ntofx;
		}else if(coef[a]==-1){
			viol= false;
		}
	}
	if(!viol) return false;
	if(viol && ntofx==0){/*std::cout<<"GlobalCut::check_updt_Viol:: infeas"<<std::endl;*/ infeas= true; return false;} //abort;
	return true;

}

//----------------------------------------------------------------------------------

bool 
LocalCut::check_updt_Viol(const double *y, bool & infeas) {
	if(type==2)check_updt_Viol2(y,infeas);
	int arc;
	int sumzro=0;
	double sum=0; 
	double rhs_ = rhs;
    rhs_dimsh=0;
    for(int a=0;a<size;++a){
    	arc = vars[a];
    	if( y[arc]==1){
        	rhs_dimsh+= 1.0;
		}else if( y[arc]==-1){
			sum+= 1.0;
		}else sumzro+=1;
    }
	rhs_ -= rhs_dimsh;
	infeas=false;
	if(sense==1){
		if(sum<rhs_){ infeas=true; return false;}  
		if(rhs_ <= 0){ return false;}
    }else{
 		if(rhs_ < 0){infeas=true; return false;}
 		if(sumzro>=(size-1)){return false;}
    }
    return true;
}


//----------------------------------------------------------------------------------

double
LocalCut::coef_at(int pos) const{
    if(pos >= size) return -1;
    if(coef) return coef[pos];
    return 1.0;
}


//----------------------------------------------------------------------------------

double
LocalCut::get_total_rhs() const{ return  rhs - rhs_dimsh;}

//----------------------------------------------------------------------------------

double
LocalCut::get_rhs() const{  return rhs; }

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
    std::cout<<"T"<<type<<" sens: "<<sense<<": ";
    for(int i=size;i--;)
        std::cout<<"("<<vars[i]<<", mtl: 1) ";
    
    std::cout<<" rhs: "<<rhs_<<std::endl;
}

