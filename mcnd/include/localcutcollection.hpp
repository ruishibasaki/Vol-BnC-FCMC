//
//  localcutcollection.hpp
//  
//
//  Created by Rui Shibasaki on 26/03/2020.
//

#ifndef localcutcollection_hpp
#define localcutcollection_hpp

#include <stdio.h>
#include <algorithm>
#include <map>
#include "MCND_cut.hpp"
#include "BCP_buffer.hpp"
#include "UtilsMethods.hpp"

class LocalCut;

class LocalCutCollection{
public:
    int sizeOfIdSeq;
    int sizeOfMap;
    int sizeOfCollection;
    int discarted;
    bool empty;
    Pair2 *map;
    LocalCut *begin;
    LocalCut *end;
    
    //--------------------------------------
    //  Construction/Destruction methods
    //---------------------------------------
    
    void initialize(int M);
    inline LocalCutCollection():map(0),begin(0),end(0){
        discarted = sizeOfIdSeq = sizeOfMap = sizeOfCollection = 0; empty=true;}
    
    ~LocalCutCollection();
    
    LocalCut * createNewLocalCut(const std::vector<int>& c, int id_vi, int serial_num_, int sense_, double rhs_);
    LocalCut * createNewLocalCut(int sz,  int* vars_, int id_vi, 
    								int serial_num_);
    LocalCut * createNewLocalCut(int sz,  int* vars_, double * coef_, double rhs_, int id_vi, 
    								int serial_num_);

    //--------------------------------------
    //  insert/del methods
    //---------------------------------------
    
    int addLocalCut(LocalCut * tryC);
    void insert_front(LocalCut * tryC);
    void insert_end(LocalCut * tryC);
    void replace(LocalCut * out, LocalCut * in);
    void pop_back_nodel();
    LocalCut * remove_nodel(LocalCut * trgt);

    //--------------------------------------
    //  VI methods
    //---------------------------------------
    int collected(LocalCut * tryC);
    
    //--------------------------------------
    //  auxiliary methods
    //---------------------------------------
    void print();
    void advance(LocalCut*& C, int n); 
    LocalCut* operator[](int n);
    void mapLocalCut(LocalCut * c1, LocalCut * c2, std::vector<PairF> & mapset);
    double LocalCut_hasArc(const LocalCut * localc, int arc);
	void map_collection(std::map<int, int>& mapd);
    //--------------------------------------
    // modifying methods
    //--------------------------------------
    
    int desactvLocalc(int lim, int * actvS, int & actvSSz, int& num_actv, double * pstarv, double * dstaru, double * dualu);
    bool swap(LocalCut * c1, LocalCut * c2);
    int swap_toend_destruct(LocalCut * out, LocalCut *& ret, bool destruct);
    LocalCut* swap_to_end(LocalCut * trgt);
    LocalCut* move_to_end(LocalCut * trgt);

};

//=================================================================================================
//=================================================================================================
//=================================================================================================


class LocalCut: public MCND_CutUnit{
public:
    LocalCut* next;
    LocalCut* prev;

    int id_vi;
    int serial_nmbr;
    int n_zerom;
    int n_nviol;
    bool toadd;
    bool prgbl;
    
    int size;
    int *vars;
    double *coef;
    
    int sense;
    int type;
    
    unsigned int *id_seq;
    double hs;
    double rhs;
    double rhs_dimsh;
	
	//--------------------------------------

    void addArc(int iset, int arc);
    void removeArc(int iset, int arc);
    bool hasArc(int iset, int arc) const ;
    double coef_at(int pos)const;

    //--------------------------------------
 
    void print();
    double get_total_rhs() const;
    double get_rhs() const;
  
    //--------------------------------------

    bool check_updt_Viol(const double *y, bool & infeas);
	bool check_updt_Viol2(const double *y, bool & infeas);
    //implemetn map
    inline LocalCut(int M, int sz, int id_vi_, int serial_nmbr_, int sense_, int type_):
    				size(sz), next(0), prev(0), coef(0), vars(0){
        id_seq = new unsigned int[M];
        std::fill(id_seq, id_seq+M, 0);
        rhs=0.0;
        n_nviol = 0;
        n_zerom = 0;
        id_vi = id_vi_;
        rhs_dimsh = 0.0;
        hs=0;
        prgbl=false;
        toadd = true;
        sense=sense_;
		serial_nmbr = serial_nmbr_;
		type = type_;
    }
    inline ~LocalCut(){
        delete [] id_seq;
        if(size>0){ delete [] vars; }
        if(coef)delete [] coef;
    }
};


//=================================================================================================
//=================================================================================================
//=================================================================================================

//MCND_cut.cpp
class LocalCCut : public MCND_Cut {

	LocalCut* localc;
	
public:
	
	inline LocalCCut(LocalCut* c): localc(c){ type = 2;}
	inline LocalCCut(): localc(0){ type =2;}
	inline ~LocalCCut(){}

	
	inline void pack(BCP_buffer& buf) const{ buf.pack(localc); }
	inline void unpack(BCP_buffer& buf){ buf.unpack(localc);}
	
	inline LocalCut* get_localc(){ return localc;}
	inline void set_localc(LocalCut* c){ localc = c;}
	
	bool check_viol_updt_fix2(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd); 
 	bool check_viol_updt_fix(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd);
	bool purgbl(){return localc->prgbl;}
	void mark_unpurgbl(){localc->prgbl=false; }
	void mark_purgbl(){localc->prgbl=true; }

	int id_vi(){return localc->id_vi;}
	int serial_nmbr(){return localc->serial_nmbr;}
 	void print(){localc->print();}


};

#endif /* localcutcollection_hpp */
