//
//  cutsetcollection.hpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#ifndef cutsetcollection_hpp
#define cutsetcollection_hpp

#include <stdio.h>
#include "Structures.hpp"
#include "UtilsMethods.hpp"

class CutSet{
public:
    CutSet* next;
    int id;
    int s_size;
    int ss_size;
    int s_ssize;
    
    int *SS_arcs;
    int *S_Sarcs;

    unsigned int *SS_nodes;
    
    double uss, us_s;
    double dss, ds_s;
    
    
    void addNode(int iset, int node);
    void removeNode(int iset, int node);
    bool nodeInS(int iset, int node) const;
    void copyref(bool stosb, int& ssb_sz, int *& ssb_arcs, double &dss_ ) const;
    void copyref(bool stosb, int& ssb_sz, int *& ssb_arcs, double &dss_, double &uss_ ) const;
    void copyref(bool stosb, int& sb_ssz, int& ssb_sz, int *& sb_sarcs, int *& ssb_arcs, double &dss) const;

    
    inline CutSet(int M): next(0), ss_size(0), s_size(0), SS_arcs(0), S_Sarcs(0), s_ssize(0){
        SS_nodes = new unsigned int[M];
        id = -1;
        for(;M--;)SS_nodes[M]=0;
    }
    inline ~CutSet(){ delete [] SS_nodes;
        if(ss_size>0) delete [] SS_arcs;
        if(s_ssize>0) delete [] S_Sarcs;
    }
};

//=================================================================================================
//=================================================================================================
//=================================================================================================

class CutSetCollection{
public:
    unsigned int last_mirror;
    int sizeOfContainer;
    int sizeOfCollection;
    int sizeOfMap;
    int discarted;
    bool empty;
    Pair2 *map;
    CutSet *begin;
    CutSet *end;
    
    const bool collected(const CutSet * trySS_); 
    CutSet * collect(const std::vector<int>& s_nodes);

    
    int addCutSet(CutSet * newSS_, const std::deque<int>& ss_arcs, double uss, double dss, const std::deque<int>& s_sarcs, double us_s, double ds_s);
    void advance(CutSet*& SS_, int n);
    void initialize(int M);
    
    inline bool inS(CutSet* SS_, int node) const{
        return SS_->nodeInS(map[node].fst,map[node].snd);
    }
    
    void print() const;
    inline CutSetCollection(): map(0), begin(0), end(0){
        sizeOfContainer = sizeOfCollection = sizeOfMap = discarted = 0; empty=true;}
    CutSetCollection(int M){initialize(M);}
    ~CutSetCollection();
    
};

//=================================================================================================
//=================================================================================================
//=================================================================================================

class MinCardCS{
public:
    MinCardCS * next;
    MinCardCS * prev;
    
    int id, csid;
    int size;
    int card;
    int actv;

    double uss, dss;
    const int * arcs;
    
    unsigned int * bin;
    
    //--------------------------------------
    //  card methods
    //---------------------------------------
    void mark(int binpos, int arcpos);

    MinCardCS(double uss_, double dss_, int id_, int rhs, int ss_size, const int * SS_arcs, int id_owner);
    ~MinCardCS(){if(bin!=0){ delete [] bin; bin =0;}}
};

//=================================================================================================
//=================================================================================================
//=================================================================================================

class MinCardCSCollection{
public:
    MinCardCS *begin;
    MinCardCS *end;
    
    int sizeOfContainer;
    int sizeOfCollection;
    int sizeOfMap;
    Pair2 *map;
    bool empty;
    
    void initialize(int M);
    inline MinCardCSCollection(){}
    inline MinCardCSCollection(int M){initialize(M);}
    ~MinCardCSCollection();
    
    
    //--------------------------------------
    //  aux methods
    //---------------------------------------
    int equalbin(const MinCardCS * card1, const MinCardCS * card2);
    void clearbin(MinCardCS * card);
    void mark(int arc, MinCardCS * card);

    //--------------------------------------
    //  VI methods
    //---------------------------------------
    int addMinCardCS(MinCardCS * newvi);
    int collected(MinCardCS * tryvi);

};


#endif /* cutsetcollection_hpp */
