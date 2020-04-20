//
//  globalcutcollection.cpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#include "globalcutcollection.hpp"


void
GlobalCutCollection::initialize(int M){
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

     end = begin =0;
}

//----------------------------------------------------------------------------------

GlobalCut *
GlobalCutCollection::createNewGlobalCut(int sz,  int* vars_,  double * coef_, int id_vi, int serial_num_, int sense_, double rhs_, int type_){
    int arc;
     GlobalCut * newC = new GlobalCut(sizeOfIdSeq, sz,  id_vi, serial_num_, sense_, type_);
    newC->vars = vars_; 
    newC->coef = coef_;
    newC->rhs = rhs_;
    
    if(type_==1){
    	for(int n=sz;n--;){
			arc = vars_[n];
			if(coef_[n]==-1)newC->addArc(map[arc].fst, map[arc].snd);
 			 
		}
    }else{
		for(int n=sz;n--;){
			arc = vars_[n];
			newC->addArc(map[arc].fst, map[arc].snd);
			 
		}
    }
    vars_=0;
    coef_=0;
    return newC;
}

//----------------------------------------------------------------------------------

GlobalCutCollection::~GlobalCutCollection(){
     for(int i=track.size();i--;){
     	delete track[i];
     }
     track.clear();
    if(map) delete [] map;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  VI methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int
GlobalCutCollection::collected(GlobalCut * tryC){
    int ret=0;
    int vi1sz=0;
    int vi2sz=0;
    bool equal;
    const GlobalCut* C;
    //std::cout<<"coll: "<<sizeOfCollection<<std::endl;
    for(int i=track.size();i--;){
    	C = track[i];
        vi1sz = tryC->get_total_sz();
        vi2sz = C->get_total_sz();            
        if( vi1sz == vi2sz && tryC->sense == C->sense && tryC->type == C->type){
            equal = true;
            for(int id=0;id<sizeOfIdSeq;++id){
                if(tryC->id_seq[id]!=C->id_seq[id]){
                    equal=false;
                    break;
                }
            }
            
            if(equal){  return 2; }
        }
    }
    return 0;
}


//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  insert/del methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int 
GlobalCutCollection::addGlobalCut(GlobalCut * tryC){
    int arc, ret;
            

    ret = collected(tryC);
    if(ret==2){ delete tryC; return 0;}
    else if(ret==0){
        //std::cout<<"local add: "; tryC->print();
        track.push_back(tryC);
        return 1;
    } 
    delete tryC; 
    return 0;
}

//----------------------------------------------------------------------------------

void 
GlobalCutCollection::insert_end(GlobalCut * tryC){
	++sizeOfCollection;
	//std::cout<<"ok"<<std::endl;
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
GlobalCutCollection::insert_front(GlobalCut * tryC){
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
GlobalCutCollection::replace(GlobalCut * out, GlobalCut * in){
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
GlobalCutCollection::pop_back_nodel(){
	if(sizeOfCollection==1){
		begin = 0;
		empty = true;
	}
	end = end->prev;
	if(end) end->next =0;
	--sizeOfCollection;

}

//----------------------------------------------------------------------------------

GlobalCut * 
GlobalCutCollection::remove_nodel(GlobalCut * trgt){
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
	GlobalCut * ret  = trgt->next;
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
GlobalCutCollection::desactvGlobalc(int lim, int * actvS, int & actvSSz, int& num_actv, double * pstarv, double * dstaru, double * dualu){
    int idx_end, idx_out;
    int curr_num = num_actv;
    GlobalCut * ret=0;
    
    //std::cout<<"remove check "<<std::endl;
   	GlobalCut *vi = begin;
   	GlobalCut * last_actv;
    for(int i=curr_num;i--;){
        idx_out = actvS[vi->id_vi];
        //std::cout<<"try vi: "<<vi->serial_nmbr<<std::endl;
        if((vi->n_nviol >= lim /*&& dualu[idx_out]<=0 */&& dstaru[idx_out]<=0)){
            //std::cout<<"remove vi: "<<vi->id_vi<<" "<<" srial: "<<vi->serial_nmbr<<std::endl;
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
            
            vi = last_actv;
        }else vi = vi->next;
    }
   
    return 0;
}

//----------------------------------------------------------------------------------

bool
GlobalCutCollection::swap(GlobalCut * c1, GlobalCut * c2){
	GlobalCut* aux;
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

GlobalCut*
GlobalCutCollection::move_to_end(GlobalCut * trgt){
	GlobalCut * ret;
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

GlobalCut*
GlobalCutCollection::swap_to_end(GlobalCut * trgt){
    if(trgt == end){return 0; }
    GlobalCut * aux; GlobalCut * ret;
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
GlobalCutCollection::swap_toend_destruct(GlobalCut * out, GlobalCut *& ret, bool destruct){
    GlobalCut* new_end=0;
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
     
    if(destruct) delete end;
    end = new_end;
    return id;
}

//------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//  auxiliary methods
//----------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

GlobalCut*
GlobalCutCollection::operator[](int n){
    if(n>=sizeOfCollection){std::cout<<"PROBLEM::GlobalCutCollection::operator[]   n>=sizeOfCollection "<<std::endl; return 0;}
    GlobalCut* C = begin;
    for(int i=1;i<=n;++i)
        C = C->next;
    return C;
}

//----------------------------------------------------------------------------------

void
GlobalCutCollection::advance(GlobalCut*& C, int n){
    if(n>sizeOfCollection) return;
    C = begin;
    for(int i=1;i<=n;++i)
        C = C->next;
}

//----------------------------------------------------------------------------------

void
GlobalCutCollection::print(){
    GlobalCut* c = begin;
    for(int i=0;i<sizeOfCollection;++i){
        std::cout<<"GlobalCut "<<i<<" "<<std::endl;
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
GlobalCutCollection::GlobalCut_hasArc(const GlobalCut * globalCut, int arc){
    if(!globalCut->hasArc(map[arc].fst, map[arc].snd)) return 0;
    
    if(globalCut->type==1){
    	return globalCut->coef[arc];
    }else{
    	for(int n=globalCut->size;n--;)
    		if(globalCut->vars[n]==arc) return globalCut->coef[n];
    }
    
    
    return 0.0;
}

//----------------------------------------------------------------------------------

void 
GlobalCutCollection::map_collection(std::map<int, int>& mapd){
	int sz = sizeOfCollection;
	GlobalCut *vi = begin;
	for(;sz--;){
		mapd.insert(std::pair<int,int>(vi->serial_nmbr,vi->id_vi));
		vi = vi->next;
	}
}

//====================================================================================
//====================================================================================
// Class GlobalCut methods
//====================================================================================
//====================================================================================

bool 
GlobalCut::check_viol_updt_fix(const BCP_vec<BCP_var*>& vbd, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd){
 	double dimsh=0;
    int sz = size;
    int arc;
    int ntofx=0;
	int tofix;
    
    viol=true;
    rhs_dimsh=0;
    
    for(int a=0;a<sz;++a){
     	arc = vars[a];
		if(vbd[arc]->lb() > 0.5 || fixd[arc]==1){
        	dimsh+= coef[a];
        	if(coef[a]==1){
        		viol= false;
        	}
		}else if(vbd[arc]->ub() > 0.5 && fixd[arc]==-1){
 			tofix=a;
			++ntofx;
		}else if(coef[a]==-1){
        	viol= false;
        }
    }
 	rhs_dimsh = dimsh;
  	if(viol && ntofx==0){ return false;} //abort;
	if(viol && ntofx==1){
		viol= false;
		if(coef[ntofx]==1){
 			arc =vars[tofix]; 
			fixd[arc]=1;
			std::cout<<"fix "<<arc<<" to 1 " <<std::endl;
 			var_changed_pos.push_back(arc);
			var_new_bd.push_back(1.0);
			var_new_bd.push_back(1.0);
		}else if(coef[ntofx]==-1){
			arc =vars[tofix]; 
			fixd[arc]=0;
			std::cout<<"fix "<<arc<<" to 0 "<<std::endl;
			zrofx=true;
 			var_changed_pos.push_back(arc);
			var_new_bd.push_back(0.0);
			var_new_bd.push_back(0.0);
		}
	}
	return true;
}

//----------------------------------------------------------------------------------

bool 
GlobalCut::check_updt_Viol(const double *y, bool & infeas)  {
	double dimsh=0;
    int sz = size;
    int arc;
    int ntofx=0;
 	
    double viol=true;
    rhs_dimsh=0;
    
    for(int a=0;a<sz;++a){
     	arc = vars[a];
     	//std::cout<<arc<<" "<<y[arc]<<std::endl;
		if(y[arc]==1){
        	dimsh+= coef[a];
        	if(coef[a]==1){
        		viol= false;
        	}
		}else if(y[arc]==-1){
 			++ntofx;
		}else if(coef[a]==-1){
        	viol= false;
        }
    }
 	rhs_dimsh = dimsh;
	if(!viol) return false;
 	if(viol && ntofx==0){/*std::cout<<"GlobalCut::check_updt_Viol:: infeas"<<std::endl;*/ infeas= true; return false;} //abort;
	return true;
}

//----------------------------------------------------------------------------------

bool 
GlobalCut::check_viol(const BCP_vec<BCP_var*>& vbd) {
 	int arc=0;
    for(int a=0;a<size;++a){
      	arc = vars[a];
      	if(vbd[arc]->lb()>0.5 && coef[a]==1)
      		 return false;
        if(vbd[arc]->ub()<0.5 && coef[a]==-1)
        	return false;
    }
    return true;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

int
GlobalCut::at(int pos) const{
    if(pos >= get_total_sz()) return -1;
     return vars[pos];
}

//----------------------------------------------------------------------------------

double
GlobalCut::gamma_at(int pos) const{
    if(pos >= size) return -1;
    return coef[pos];
}

//----------------------------------------------------------------------------------


void
GlobalCut::get_total_sz_rhs(int & sz, double &rhs) const{
    sz = size;
    rhs = rhs- rhs_dimsh;
}

//----------------------------------------------------------------------------------

double
GlobalCut::get_total_rhs() const{ return  rhs - rhs_dimsh;}

//----------------------------------------------------------------------------------

double
GlobalCut::get_rhs() const{  return rhs; }

//----------------------------------------------------------------------------------

int
GlobalCut::get_total_sz()const{    return size;}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

void 
GlobalCut::addArc(int iset, int arc){
    setBit(id_seq,iset,arc);
}

//----------------------------------------------------------------------------------


void 
GlobalCut::removeArc(int iset, int arc){
    clearBit(id_seq,iset,arc);
}

//----------------------------------------------------------------------------------

bool 
GlobalCut::hasArc(int iset, int arc) const{
    if(testBit(id_seq,iset,arc)) return true;
    else return false;
}

//----------------------------------------------------------------------------------

void 
GlobalCut::print(){
    double rhs_ = rhs;
    std::cout<<"G"<<type<<": ";
    for(int i=size;i--;)
        std::cout<<"("<<vars[i]<<", mtl: "<<coef[i]<<") ";
    
    std::cout<<" rhs: "<<rhs_<<" nmbr: "<<serial_nmbr<<std::endl;
}

