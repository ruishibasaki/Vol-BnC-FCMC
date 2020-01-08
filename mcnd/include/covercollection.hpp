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
#include "BCP_cut.hpp"


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
    Cover *begin_actv;
    Cover *end_actv;
    int newly_added;
    //--------------------------------------
    //  Construction/Destruction methods
    //---------------------------------------
    
    void initialize(int M);
    inline CoverCollection():map(0),begin(0),end(0), begin_actv(0), end_actv(0){
        discarted = sizeOfIdSeq = sizeOfMap = sizeOfCollection = newly_added=0; empty=true;}
    
    ~CoverCollection();
    
    //--------------------------------------
    //  setter/getter methods
    //---------------------------------------
    Cover * createNewCover(const std::deque<Pair2>& c, double mu, int id_vi, int id_owner_, int serial_num_);
    int addCover(Cover * tryC, const double * xystar);
    void replace(Cover * out, Cover * in);
    //--------------------------------------
    //  VI methods
    //---------------------------------------
    int collected(Cover * tryC);
    int compScalar(Cover* c1, Cover* c2);
    //--------------------------------------
    //  auxiliary methods
    //---------------------------------------
    const void print();
    const void advance(Cover*& C, int n);
    Cover* operator[](int n);
    void mapCover(Cover * c1, Cover * c2, std::vector<PairF> & mapset);
    double cover_hasArc(const Cover * cover, int arc);

    //--------------------------------------
    // modifying methods
    //--------------------------------------
    
    int swap_toend_destruct(Cover * out, Cover *& ret, bool destruct);
    int removeCover(int lim, int * actvS, int & actvSSz, double * pstarv, double * dstaru, double * dualu);
    Cover* swap_to_end(Cover * trgt);
    Cover* move_to_end(Cover * trgt);

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


class Cover : public BCP_cut_algo{
public:
    Cover* next;
    Cover* prev;
    CoverL * Lftd;
    int id_vi;
    int serial_nmbr;
    int n_zerom;
    int n_nviol;
    
    int size;
    int maxsize;
    int *C;
    const int * owner;
    int id_owner;
    
    unsigned int *id_seq;
    bool hasLftd;
    double hs;
    double rhs;
    double rhs_dimsh;

    void addArc(int iset, int arc);
    void removeArc(int iset, int arc);
    bool hasArc(int iset, int arc) const ;
    int at(int pos)const ;
    double gamma_at(int pos)const;
    
    void print();
    void addLftd(CoverL * l);
    void get_total_sz_rhs(int & sz, double &rhs)const;
    int get_total_sz()const;
    double get_total_rhs() const;
    //implemetn map
    inline Cover(int M, int sz, int id_vi_, int id_owner_, int serial_nmbr_):size(sz), Lftd(0), next(0), prev(0), BCP_cut_algo(0, 1e40){
        C = new int[sz];
        id_seq = new unsigned int[M];
        for(;M--;)id_seq[M]=0;
        hasLftd=false;
        rhs=1.0;
        n_nviol = 0;
        n_zerom = 0;
        id_vi = id_vi_;
        id_owner = id_owner_;
        rhs_dimsh = 0.0;
        hs=0;
		serial_nmbr = serial_nmbr_;
    }
    inline ~Cover(){
        delete [] id_seq;
        if(size>0){ delete [] C; }
        if(hasLftd) delete Lftd;
    }
};



#endif /* covervi_hpp */
