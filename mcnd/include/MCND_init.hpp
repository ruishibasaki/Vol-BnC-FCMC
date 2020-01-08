// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_INIT_H
#define _MCND_INIT_H

#include "BCP_USER.hpp"
#include <BCP_warmstart_dual.hpp> 
#include <CoinWarmStartDual.hpp>
#include "UtilsMethods.hpp"

class MCND_packer : public BCP_user_pack {

public:

    virtual void 
    pack_user_data(const BCP_user_data* ud, BCP_buffer& buf);
    virtual BCP_user_data* 
    unpack_user_data(BCP_buffer& buf);
   
};


class MCND_initialize : public USER_initialize {
  // Declare this function if not the default single process communication is
  // wanted
  //   BCP_message_environment * msgenv_init();
  BCP_user_pack* packer_init(BCP_user_class* p);
  BCP_tm_user* tm_init(BCP_tm_prob& p,
		       const int argnum, const char * const * arglist);
  BCP_lp_user* lp_init(BCP_lp_prob& p);
};

#endif
