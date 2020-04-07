//
//  globalcutcollection.hpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#ifndef globalcutcollection_hpp
#define globalcutcollection_hpp

#include <stdio.h>
#include <algorithm>
#include <map>
#include "MCND_cut.hpp"
#include "BCP_buffer.hpp"
#include "BCP_var.hpp"
#include "BCP_vector.hpp"
#include "UtilsMethods.hpp"

class GlobalCut;

class GlobalCutCollection{
public:
    int sizeOfIdSeq;
    int sizeOfMap;
    int sizeOfCollection;
    bool empty;
    Pair2 *map;
    GlobalCut *begin;
    GlobalCut *end;
    
    std::deque<GlobalCut *> track;
    //--------------------------------------
    //  Construction/Destruction methods
    //---------------------------------------
    
    void initialize(int M);
    inline GlobalCutCollection():map(0),begin(0),end(0){
         sizeOfIdSeq = sizeOfMap = sizeOfCollection = 0; empty=true;}
    
    ~GlobalCutCollection();
    
    GlobalCut * createNewGlobalCut(int sz, int* vars_,  double * coef_, int id_vi, 
    								int serial_num_, int sense_, double rhs_, int type_);

    //--------------------------------------
    //  insert/del methods
    //---------------------------------------
    
    int addGlobalCut(GlobalCut * tryC);
    void insert_front(GlobalCut * tryC);
    void insert_end(GlobalCut * tryC);
    void replace(GlobalCut * out, GlobalCut * in);
    void pop_back_nodel();
    GlobalCut * remove_nodel(GlobalCut * trgt);

    //--------------------------------------
    //  VI methods
    //---------------------------------------
    int collected(GlobalCut * tryC);
    
    //--------------------------------------
    //  auxiliary methods
    //---------------------------------------
    void print();
    void advance(GlobalCut*& C, int n); 
    GlobalCut* operator[](int n);
    void mapGlobalCut(GlobalCut * c1, GlobalCut * c2, std::vector<PairF> & mapset);
    double GlobalCut_hasArc(const GlobalCut * localc, int arc);
	void map_collection(std::map<int, int>& mapd);
    //--------------------------------------
    // modifying methods
    //--------------------------------------
    
    int desactvGlobalc(int lim, int * actvS, int & actvSSz, int& num_actv, double * pstarv, double * dstaru, double * dualu);
    bool swap(GlobalCut * c1, GlobalCut * c2);
    int swap_toend_destruct(GlobalCut * out, GlobalCut *& ret, bool destruct);
    GlobalCut* swap_to_end(GlobalCut * trgt);
    GlobalCut* move_to_end(GlobalCut * trgt);

};

//=================================================================================================
//=================================================================================================
//=================================================================================================


class GlobalCut: public MCND_CutUnit{
public:
    GlobalCut* next;
    GlobalCut* prev;

    int id_vi;
    int serial_nmbr;
    int n_zerom;
    int n_nviol;
    int type;
    
    int size;
    int *vars;
    double * coef;
    int sense;
    bool purgbl;
    
    unsigned int *id_seq;
    double hs;
    double rhs;
    double rhs_dimsh;
	
	//--------------------------------------

    void addArc(int iset, int arc);
    void removeArc(int iset, int arc);
    bool hasArc(int iset, int arc) const ;
    int at(int pos)const ;
    double gamma_at(int pos)const;

    //--------------------------------------
 
    void print();
    void get_total_sz_rhs(int & sz, double &rhs)const;
    int get_total_sz()const;
    double get_total_rhs() const;
    double get_rhs() const;
  
    //--------------------------------------

    bool check_viol(const BCP_vec<BCP_var*>& vbd);
    bool check_updt_Viol(const double *y, bool & infeas);
	bool check_viol_updt_fix(const BCP_vec<BCP_var*>& vbd, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd);
    //implemetn map
    inline GlobalCut(int M, int sz, int id_vi_, int serial_nmbr_, int sense_, int type_):size(sz), next(0), prev(0){
        id_seq = new unsigned int[M];
        std::fill(id_seq, id_seq+M, 0);
        rhs=0.0;
        n_nviol = 0;
        n_zerom = 0;
        id_vi = id_vi_;
        rhs_dimsh = 0.0;
        hs=0;
        purgbl=false;
        sense=sense_;
        type = type_;
		serial_nmbr = serial_nmbr_;
    }
    inline ~GlobalCut(){
        delete [] id_seq;
        if(size>0){ delete [] vars; delete []  coef;}
    }
};


#endif /* localcutcollection_hpp */
