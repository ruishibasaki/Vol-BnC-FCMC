//
//  covercollection_hpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#ifndef covercollection_hpp
#define covercollection_hpp

#include <stdio.h>
#include "Structures.hpp"
#include "UtilsMethods.hpp"
#include "MCND_cut.hpp"
#include "BCP_buffer.hpp"
#include <map>

class CoverVI;
class Cover;

class CoverCollection{
public:
    int sizeOfIdSeq;
    int sizeOfMap;
    int sizeOfCollection;
    int discarted;
    bool empty;
    Pair2 *map;
    Cover *begin;
    Cover *end;
    
    //--------------------------------------
    //  Construction/Destruction methods
    //---------------------------------------
    
    void initialize(int M);
    inline CoverCollection():map(0),begin(0),end(0){
        discarted = sizeOfIdSeq = sizeOfMap = sizeOfCollection = 0; empty=true;}
    
    ~CoverCollection();
    
    Cover * createNewCover(const std::deque<Pair2>& c, int id_vi, int serial_num_);

    //--------------------------------------
    //  insert/del methods
    //---------------------------------------
    
    int addCover(Cover * tryC, const double * xystar);
    void insert_front(Cover * tryC);
    void insert_end(Cover * tryC);
    void replace(Cover * out, Cover * in);
    void pop_back_nodel();
    Cover * remove_nodel(Cover * trgt);
    int desactvCover(int lim, int * actvS, int & actvSSz, int& num_actv, double * pstarv, double * dstaru, double * dualu);

    //--------------------------------------
    //  VI methods
    //---------------------------------------
    int collected(Cover * tryC);
    int compScalar(Cover* c1, Cover* c2);
    
    //--------------------------------------
    //  auxiliary methods
    //---------------------------------------
    void print();
    void advance(Cover*& C, int n); 
    Cover* operator[](int n);
    void mapCover(Cover * c1, Cover * c2, std::vector<PairF> & mapset);
    double cover_hasArc(const Cover * cover, int arc);
	void map_collection(std::map<int, int>& mapd);
    //--------------------------------------
    // modifying methods
    //--------------------------------------
    
    int swap_toend_destruct(Cover * out, Cover *& ret, bool destruct);
    Cover* swap_to_end(Cover * trgt);
    Cover* move_to_end(Cover * trgt);
    bool swap(Cover * c1, Cover * c2);

};

//=================================================================================================
//=================================================================================================
//=================================================================================================

class CoverL{
public:
    const Cover* owner;
    int size;
    int maxsize;
    int* SSmC;
    double* gamma;
    double ttgamma1;
    
    void addvar(int arc, double gam, bool lift_down);
    inline CoverL(int sz, const Cover* _owner):size(0),maxsize(sz), owner(_owner){
        SSmC = new int[sz];
        gamma = new double[sz];
        ttgamma1=0;
    }
    inline CoverL(int sz):size(0),maxsize(sz), owner(0){
        SSmC = new int[sz];
        gamma = new double[sz];
        ttgamma1=0;
    }
    inline ~CoverL(){if(maxsize>0){ delete [] SSmC; delete [] gamma;}}
};

//=================================================================================================
//=================================================================================================
//=================================================================================================


class Cover: public MCND_CutUnit{
public:
    Cover* next;
    Cover* prev;
    CoverL * Lftd;
    int id_vi;
    int serial_nmbr;
    int n_zerom;
    int n_nviol;
    bool prgbl;
    bool toadd;
    
    int size;
    int maxsize;
    int *C;
    const int * owner;
    
    unsigned int *id_seq;
    bool hasLftd;
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
    void addLftd(CoverL * l);
    void get_total_sz_rhs(int & sz, double &rhs)const;
    int get_total_sz()const;
    double get_total_rhs() const;
    double get_rhs() const;
  
    //--------------------------------------

    double viol(const double *y)const;
    bool check_updt_Viol(const double *y, bool & infeas);

    //implemetn map
    inline Cover(int M, int sz, int id_vi_, int serial_nmbr_):size(sz), Lftd(0), next(0), prev(0){
        C = new int[sz];
        id_seq = new unsigned int[M];
        for(;M--;)id_seq[M]=0;
        hasLftd=false;
        rhs=1.0;
        n_nviol = 0;
        n_zerom = 0;
        id_vi = id_vi_;
        rhs_dimsh = 0.0;
        hs=0;
        prgbl=false;
        toadd=true;
		serial_nmbr = serial_nmbr_;
		maxsize=0;
    }
    inline ~Cover(){
        delete [] id_seq;
        if(size>0){ delete [] C; }
        if(hasLftd) delete Lftd;
    }
};


//=================================================================================================
//=================================================================================================
//=================================================================================================

//MCND_cut.cpp
class CoverCut : public MCND_Cut{

	Cover* cover;
	
public:
	
	inline CoverCut(Cover* c): cover(c){type =1;}
	inline CoverCut():cover(0){type =1;}
	inline ~CoverCut(){}

	
	inline void pack(BCP_buffer& buf) const{ buf.pack(cover); }
	inline void unpack(BCP_buffer& buf){ buf.unpack(cover);}
	
	inline Cover* get_cover(){ return cover;}
	inline void set_cover(Cover* c){ cover = c;}
	
	bool check_viol(const BCP_vec<BCP_var*>& vars);
	double check_viol(const double* vars);
	bool check_logical_fix(const BCP_vec<BCP_var*>& vars, int* yarcs);
	bool check_viol_updt_fix(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd);
	bool purgbl(){return cover->prgbl;}
	void mark_unpurgbl(){cover->prgbl=false;}
	void mark_purgbl(){cover->prgbl=true; }

	int id_vi(){return cover->id_vi;}
	int serial_nmbr(){return cover->serial_nmbr;}
 	void print(){cover->print();}


};



#endif /* covervi_hpp */
