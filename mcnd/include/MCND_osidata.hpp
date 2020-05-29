// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_OSIDATA_H
#define _MCND_OSIDATA_H

#include <iostream>
#include <vector>
#include <deque>
#include "covermanager.hpp"
#include "cutsetmanager.hpp"
#include "localcutmanager.hpp"
#include "globalcutmanager.hpp"
#include "VolVolume.hpp"

#include "Structures.hpp"
#include "OsiAuxInfo.hpp"

#include "MCND_checklp.hpp"


//-----------------------------------------------------
// OsiAuxInfo
//-----------------------------------------------------

class OsiVolAuxInfo: public OsiAuxInfo{
public:
    const Data* data;
    CoverManager * cover_manager;
    CutSetManager * ss_manager;
    LocalCutManager* localc_manager;
    GlobalCutManager* globalc_manager;
    LPChecker* lpchecker;
    VOL_problem volprob;
    
    //OsiSolver attributes record;
    double  *colub;
    double  *collb;
    double  *rowub;
    double  *rowlb;
    
    double  *solution;
    double  *dual;
    double  *rc_;
    double  *lhs_;
    
    double * VItopo;
    double * yhit;
    int * arc_map;
    int * actv; 
 
     //Cut generation attributes
    int maxNumVI, intvlVI;
    int maxPos;
    
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
    
	void init_solver(int corevi){
		maxPos+= corevi;
		volprob.ext_initializer(data->narcs+data->narcs*data->ndemands, maxPos);
	}

//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------

    OsiVolAuxInfo():volprob("volmcnd.par") {
        maxNumVI =1000;
		maxPos = 10000;
        intvlVI =50;
        appData_ = this;
        colub= 0;
    	collb=0;
    	rowub=0;
    	rowlb=0;
    
    	solution=0;
    	dual=0;
    	rc_=0;
    	lhs_=0;
    
    	VItopo=0;
    	yhit=0;
    	arc_map=0;
    	actv=0; 
    	data =(0); cover_manager =(0); ss_manager =(0); localc_manager =(0); globalc_manager =(0);
    }  

//-----------------------------------------------------

    ~OsiVolAuxInfo(){
		if(yhit){ delete [] yhit;  yhit=0;}
		if(arc_map){ delete [] arc_map;  arc_map=0;}
		if(VItopo){ delete [] VItopo;  VItopo=0;}
		if(rowub){delete[] rowub;    rowub = 0;}
		if(rowlb){delete[] rowlb;    rowlb = 0;}
		if(colub){delete[] colub; colub = 0;}
		if(collb){delete[] collb;    collb = 0;}
		if(solution){delete[] solution;	        solution = 0;}
		if(dual){delete[] dual;	        dual = 0;}
		if(rc_){delete[] rc_;     rc_ = 0;}
		if(lhs_){delete[] lhs_;    lhs_ = 0;}
		if(actv){delete[] actv; actv=0;}
    } 
    
    
   
 
};




#endif
