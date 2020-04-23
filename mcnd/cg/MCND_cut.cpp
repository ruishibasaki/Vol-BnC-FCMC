//
//  MCND_cut.hpp
//  
//
//  Created by Rui Shibasaki on 08/01/2020.
//


#include "MCND_cut.hpp"

#include "covercollection.hpp"
#include "localcutcollection.hpp"

int CutManager::ttgend = 0;

bool 
CoverCut::check_viol_updt_fix(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol,bool & zrofx, int* fixd){
	double sum=0;
	double dimsh=0;
    double rhs= cover->get_rhs();
    int sz = cover->get_total_sz();
    int arc;
    int ntofx=0;
    std::vector<int> tofix(sz);
    
    viol=true;
    cover->rhs_dimsh=0;
    
    for(int a=0;a<sz;++a){
     	arc = cover->at(a);
		if(vars[arc]->lb() > 0.5 || fixd[arc]==1){
        	dimsh+= cover->gamma_at(a);
		}else if(vars[arc]->ub() > 0.5 && fixd[arc]==-1){
			sum+= cover->gamma_at(a);
			tofix[ntofx++]=arc;
		}
    }
	rhs -= dimsh;
	cover->rhs_dimsh = dimsh;
	if(sum<rhs){ tofix.clear();return false;} //abort;
	if(rhs <= 0){viol= false;}
	else if(sum==rhs){
		viol= false;
		for(;ntofx--;){
			arc = tofix[ntofx];
			fixd[arc]=1;
			std::cout<<"fix "<<arc<<" to 1 "<<std::endl;
 			var_changed_pos.push_back(arc);
			var_new_bd.push_back(1.0);
			var_new_bd.push_back(1.0);
		}
	}
	return true;
}

//-------------------------------------------------------------------------------------------

bool
CoverCut::check_viol(const BCP_vec<BCP_var*>& vars){
	double sum=0;
    double comp= cover->get_rhs();
    int sz = cover->get_total_sz();
    cover->rhs_dimsh=0;
    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	if(vars[cover->at(a)]->lb() > 0.5)
        	sum+= cover->gamma_at(a);
        if(sum>=comp){return false;}
    }
    cover->rhs_dimsh = sum;
    return true;
}

//-------------------------------------------------------------------------------------------

double 
CoverCut::check_viol(const double* vars){
	return cover->viol(vars);
}

//-------------------------------------------------------------------------------------------
// 2 fix open
// 3 fix zero
bool 
CoverCut::check_logical_fix(const BCP_vec<BCP_var*>& vars, int* yarcs){
	double sum=0;
    double rhs= cover->get_rhs();
    int sz = cover->get_total_sz();
    int dimsh=0;
    int arc;
    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	arc = cover->at(a);
		if(/*vars[arc]->lb() > 0.5 ||*/ yarcs[arc]==1 || yarcs[arc]==2){
        	dimsh+= cover->gamma_at(a);
		}else if(/*vars[arc]->ub() > 0.5 &&*/ yarcs[arc]==-1){
			sum+= cover->gamma_at(a);
		}
    }
	rhs -= dimsh;
	if(sum == rhs && rhs>0){
		for(int a=0;a<sz;++a){
			arc = cover->at(a);
			if(vars[arc]->lb() < 0.5 && vars[arc]->ub() > 0.5){
				 yarcs[arc] = 2;
			}
		}
		std::cout<<"serial: "<<cover->serial_nmbr<<" "; cover->print();
		return true;
	}else if(sum < rhs){std::cout<<"check_logical_fix::PRUNE: "<<sum<<" < "<<rhs<<" "; cover->print(); return false;}
	return true;
}

//======================================================================
//======================================================================
// LocalCCut methods
//======================================================================
//======================================================================

bool
LocalCCut::check_viol_updt_fix2(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd){
	double dimsh=0;
	int sz = localc->size;
	int arc;
	int ntofx=0;
	int tofix;
	double coef_;
	viol=true;
	localc->rhs_dimsh=0;

	for(int a=0;a<sz;++a){
		arc = localc->vars[a];
		coef_ = localc->coef[a];
		if(vars[arc]->lb() > 0.5 || fixd[arc]==1){
			dimsh+= coef_;
			if(coef_==1){
				viol= false;
			}
 		}else if(vars[arc]->ub() > 0.5 && fixd[arc]==-1){
 			tofix= a;
			++ntofx;
		}else if(coef_==-1){
			viol= false;
		}
	} 
	localc->rhs_dimsh = dimsh;
	if(viol && ntofx==0){ return false;} //abort;
	if(viol && ntofx==1){
		viol= false;
		coef_ = localc->coef[tofix];
		if(coef_==1){
			arc = localc->vars[tofix]; 
			fixd[arc]=1;
			std::cout<<"localfix "<<arc<<" to 1 " <<std::endl;
			var_changed_pos.push_back(arc);
			var_new_bd.push_back(1.0);
			var_new_bd.push_back(1.0);
		}else if(coef_==-1){
			arc = localc->vars[tofix]; 
			fixd[arc]=0;
			std::cout<<"localfix "<<arc<<" to 0 "<<std::endl;
			zrofx=true;
			var_changed_pos.push_back(arc);
			var_new_bd.push_back(0.0);
			var_new_bd.push_back(0.0);
		}
	}
	return true;
}


//-------------------------------------------------------------------------------------------

bool
LocalCCut::check_viol_updt_fix(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd){
    
    if(localc->type==2)
    	return check_viol_updt_fix2( vars, var_changed_pos, var_new_bd,  viol,  zrofx, fixd);
    	
	double sum=0;
	double dimsh=0;
    double rhs= localc->rhs;
    int sz = localc->size;
	int sumzro =0;
    int ntofx=0;
    int arc;
    std::vector<int> tofix(sz);

	viol=true;
    localc->rhs_dimsh=0;

    for(int a=0;a<sz;++a){
    	//std::cout<<"c: "<<cover->at(a)<<" lb: "<<vars[cover->at(a)]->lb()<<std::endl;
    	arc = localc->vars[a];
    	if(vars[arc]->lb() > 0.5 || fixd[arc]==1){
        	dimsh+= 1.0;
		}else if(vars[arc]->ub() > 0.5 && fixd[arc]==-1){
			sum+= 1.0;
			tofix[ntofx++] = arc;
		}else{ sumzro+=1; }
    }
	rhs -= dimsh;
	localc->rhs_dimsh = dimsh;
    if(localc->sense==1){
		if(sum<rhs){ tofix.clear();return false;} //abort;
		if(rhs <= 0){viol= false;}
		else if(sum==rhs){
			viol= false;
			for(;ntofx--;){
				arc = tofix[ntofx];
				fixd[arc] = 1;
				std::cout<<"localfix "<<arc<<" to 1 "<<std::endl;
				var_changed_pos.push_back(arc);
				var_new_bd.push_back(1.0);
				var_new_bd.push_back(1.0);
			}
		}
    }else{
 		if(rhs < 0){tofix.clear();return false;}
 		if(rhs==0){
 			viol= false;
 			for(;ntofx--;){
				arc = tofix[ntofx];
				fixd[arc] = 0;
				zrofx=true;
				std::cout<<"localfix "<<arc<<" to 0"<<std::endl;
				var_changed_pos.push_back(arc);
				var_new_bd.push_back(0.0);
				var_new_bd.push_back(0.0);
			}
 		}else if(sumzro>=(sz-1)){viol= false;}
    }
    tofix.clear();
    return true;
}




