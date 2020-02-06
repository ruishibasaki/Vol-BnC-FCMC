// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#include <fstream>
#include <algorithm>

#include "BCP_lp_param.hpp"
#include "BCP_parameters.hpp"

#include "BCP_tm.hpp"
#include "BCP_lp.hpp"

#include "BCP_USER.hpp"

#include "MCND_branch_score.hpp"
#include "MCND_init.hpp"
#include "MCND_tm.hpp"
#include "MCND_lp.hpp"
#include "MCND_data.hpp"
//#############################################################################


void 
MCND_packer::pack_user_data(const BCP_user_data* ud, BCP_buffer& buf){
	//std::cout<<"pack_user_data"<<std::endl;
	const MCND_node_branch_data * user = dynamic_cast<const MCND_node_branch_data *>(ud);
	if(user) user->pack(buf);

}
    
//-----------------------------------------------------------------------------

BCP_user_data* 
MCND_packer::unpack_user_data(BCP_buffer& buf){
	//std::cout<<"unpack_user_data"<<std::endl;

	 MCND_node_branch_data*  user = new MCND_node_branch_data;
	 user->unpack(buf);
	 return user;
}

//-----------------------------------------------------------------------------

void 
MCND_packer::pack_cut_algo(const BCP_cut_algo* cut, BCP_buffer& buf){
	//std::cout<<"pack_cut_algo:"<<std::endl;
	const CoverCut * c = dynamic_cast<const CoverCut *>(cut);
	c->pack(buf);
}

//-----------------------------------------------------------------------------

BCP_cut_algo* 
MCND_packer::unpack_cut_algo(BCP_buffer& buf){
	//std::cout<<"unpack_cut_algo:"<<std::endl;
	CoverCut* cut= new CoverCut();
	cut->unpack(buf);
	return cut;
}


//#############################################################################


BCP_user_pack*
MCND_initialize::packer_init(BCP_user_class* p){
  return new MCND_packer;
}

//-----------------------------------------------------------------------------


BCP_tm_user *
MCND_initialize::tm_init(BCP_tm_prob& p,
		       const int argnum, const char * const * arglist){
   MCND_tm* tm = new MCND_tm;

  // MCND_read_parameters(*tm, arglist[1]);
   MCND_read_data(arglist[2], tm->data);
   p.par.set_entry(BCP_tm_par::Granularity, 1e-2);
   p.par.set_entry(BCP_tm_par::ReportWhenDefaultIsExecuted, false);
   return tm;
}

//-----------------------------------------------------------------------------

BCP_lp_user *
MCND_initialize::lp_init(BCP_lp_prob& p){
	
   MCND_lp* lp = new MCND_lp;
   p.par.set_entry(BCP_lp_par::LpVerb_GeneratedCutCount, true);
   p.par.set_entry(BCP_lp_par::LpVerb_CutsToCutPoolCount, true);
   p.par.set_entry(BCP_lp_par::LpVerb_PresolveResult, false);
   p.par.set_entry(BCP_lp_par::LpVerb_StrongBranchResult, false);
   p.par.set_entry(BCP_lp_par::LpVerb_ChildrenInfo, false);
   p.par.set_entry(BCP_lp_par::DoReducedCostFixingAtZero ,false);
   p.par.set_entry(BCP_lp_par::DoReducedCostFixingAtAnything, false);
   p.par.set_entry(BCP_lp_par::LpVerb_NodeTime, false);
   p.par.set_entry(BCP_lp_par::ReportWhenDefaultIsExecuted , false);
   p.par.set_entry(BCP_lp_par::WarmstartInfo, BCP_WarmstartParent);
   p.par.set_entry(BCP_lp_par::MessagePassingIsSerial, true);
   p.par.set_entry(BCP_lp_par::MaxCutsAddedPerIteration, 1000);
   p.par.set_entry(BCP_lp_par::MaxResolveIter, 1 );
   p.par.set_entry(BCP_lp_par::IneffectiveConstraints, BCP_IneffConstr_None);
   return lp;
}

//#############################################################################
/*
void MCND_read_parameters(MCND_tm& tm, const char * paramfile)
{
  // tm.tm_par.read_from_file(paramfile);
}*/

//#############################################################################
