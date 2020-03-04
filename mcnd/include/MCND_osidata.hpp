// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_OSIDATA_H
#define _MCND_OSIDATA_H

#include <iostream>
#include <vector>
#include <deque>
#include "covermanager.hpp"
#include "cutsetmanager.hpp"

#include "Structures.hpp"
#include "OsiAuxInfo.hpp"


//-----------------------------------------------------
// OsiAuxInfo
//-----------------------------------------------------

class OsiVolAuxInfo: public OsiAuxInfo{
public:
    const Data* data;
    CoverManager * cover_manager;
    CutSetManager * ss_manager;
    
    int maxNumVI, intvlVI;
    
    inline OsiVolAuxInfo():data(0), cover_manager(0), ss_manager(0){
        maxNumVI =1000;
        intvlVI =50;
        appData_ = this;
    }  
    
};




#endif
