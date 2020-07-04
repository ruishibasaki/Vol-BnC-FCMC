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
    int ss_ksize;
    int s_sksize;
    
    int *SS_arcs;
    int *S_Sarcs;
    int *SS_comm;
    int *S_Scomm;

    unsigned int *SS_nodes;
    
    double uss, us_s;
    double dss, ds_s;
    
    
    void addNode(int iset, int node);
    void removeNode(int iset, int node);
    bool nodeInS(int iset, int node) const;
    
    void print(bool reverse ) const;
    
    inline CutSet(int M): next(0), ss_size(0), s_size(0), SS_arcs(0), S_Sarcs(0), s_ssize(0), SS_comm(0), S_Scomm(0){
        SS_nodes = new unsigned int[M];
        id = -1;
        ss_ksize = s_sksize=0;
        for(;M--;)SS_nodes[M]=0;
    }
    inline ~CutSet(){ delete [] SS_nodes;
        if(ss_size>0) delete [] SS_arcs;
        if(s_ssize>0) delete [] S_Sarcs;
        if(ss_ksize>0) delete [] SS_comm;
        if(s_sksize>0) delete [] S_Scomm;
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

    
    int addCutSet(CutSet * newSS_, const std::deque<int>& ss_arcs, const std::deque<int>& ss_k, double uss, double dss, const std::deque<int>& s_sarcs, const std::deque<int>& s_sk, double us_s, double ds_s);
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



#endif /* cutsetcollection_hpp */
