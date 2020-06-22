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
	 MCND_node_branch_data*  user = new MCND_node_branch_data;
	 user->unpack(buf);
	 return user;
}

//-----------------------------------------------------------------------------

void 
MCND_packer::pack_cut_algo(const BCP_cut_algo* cut, BCP_buffer& buf){
	//std::cout<<"pack_cut_algo:"<<std::endl;
	const MCND_Cut * c = dynamic_cast<const MCND_Cut *>(cut);
	buf.pack(c->cut_type);
	c->pack(buf);
}

//-----------------------------------------------------------------------------

BCP_cut_algo* 
MCND_packer::unpack_cut_algo(BCP_buffer& buf){
	//std::cout<<"unpack_cut_algo:"<<std::endl;
	int type;
	buf.unpack(type);
	switch(type){
		case 1:{
			CoverCut* cut= new CoverCut();
			cut->unpack(buf);
			return cut;
		}
		case 2:{
			LocalCCut* cut= new LocalCCut();
			cut->unpack(buf);
			return cut;
		}
	}
	
	return 0;
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

   MCND_read_data(arglist[1], tm->data);
   std::string inst(arglist[1]);
   inst = inst.substr(inst.find("instances")+10); 
   tm->instance = inst;
   		 
  
   tm->init_sol.size = tm->data.narcs;
   tm->init_sol.xy = new double [tm->data.narcs];
   std::fill(tm->init_sol.xy, tm->init_sol.xy+tm->data.narcs,0.0);
   int ret = read_init_sol("optimal.txt", inst, tm->init_sol.xy, tm->init_sol.cost);
   if(ret>=0){
   		//tm->init_sol.cost = 1e30;
    	if(ret==0) tm->init_sol.onlyvalue=true;
		p.upper_bound = tm->init_sol.cost;
		MCND_solution * sol = new MCND_solution();
		sol->copy(tm->init_sol, tm->data.narcs);
   		p.feas_sol = sol;
   		
   }else{tm->init_sol=0;}

   p.par.set_entry(BCP_tm_par::VerbosityShutUp, false);
   p.par.set_entry(BCP_tm_par::TmVerb_AllFeasibleSolutionValue, false);
    p.par.set_entry(BCP_tm_par::TmVerb_AllFeasibleSolution, false);
    p.par.set_entry(BCP_tm_par::TmVerb_BetterFeasibleSolutionValue, false);
    p.par.set_entry(BCP_tm_par::TmVerb_BetterFeasibleSolution, false);
    p.par.set_entry(BCP_tm_par::TmVerb_BestFeasibleSolution, false);
    p.par.set_entry(BCP_tm_par::TmVerb_NewPhaseStart, false);
    p.par.set_entry(BCP_tm_par::TmVerb_TrimmedNum, false);
    p.par.set_entry(BCP_tm_par::TmVerb_TimeOfImprovingSolution, false);
    p.par.set_entry(BCP_tm_par::TmVerb_PrunedNodeInfo, false);
    p.par.set_entry(BCP_tm_par::TmVerb_FinalStatistics, false);
    p.par.set_entry(BCP_tm_par::TmVerb_ReportDefault, false);
    
   p.par.set_entry(BCP_tm_par::Granularity, 1e-2);
   p.par.set_entry(BCP_tm_par::ReportWhenDefaultIsExecuted, false);
   p.par.set_entry(BCP_tm_par::WarmstartInfo, BCP_WarmstartNone); //BCP_WarmstartParent;
   p.par.set_entry(BCP_tm_par::TreeSearchStrategy, BCP_BestFirstSearch); 
   p.par.set_entry(BCP_tm_par::UnconditionalDiveProbability, 1.0);
   p.par.set_entry(BCP_tm_par::QualityRatioToAllowDiving_HasUB, -1.0);
   p.par.set_entry(BCP_tm_par::QualityRatioToAllowDiving_NoUB, -1.0);
   p.par.set_entry(BCP_tm_par::RemoveExploredBranches, true);
   p.par.set_entry(BCP_tm_par::MessagePassingIsSerial, true);
	tm->time_lim = 10*3600.0;
   p.par.set_entry(BCP_tm_par::MaxRunTime, tm->time_lim);

   
   return tm; 
}

//-----------------------------------------------------------------------------

BCP_lp_user *
MCND_initialize::lp_init(BCP_lp_prob& p){
	
   MCND_lp* lp = new MCND_lp;
    
    p.par.set_entry(BCP_lp_par::LpVerb_AddedCutCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_AddedVarCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ChildrenInfo, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ColumnGenerationInfo, false);
	p.par.set_entry(BCP_lp_par::LpVerb_CutsToCutPoolCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_VarsToVarPoolCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_FathomInfo, false);
	p.par.set_entry(BCP_lp_par::LpVerb_IterationCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_RelaxedSolution, false);
	p.par.set_entry(BCP_lp_par::LpVerb_FinalRelaxedSolution, false);
	p.par.set_entry(BCP_lp_par::LpVerb_LpMatrixSize, false);
	p.par.set_entry(BCP_lp_par::LpVerb_LpSolutionValue, false);
	p.par.set_entry(BCP_lp_par::LpVerb_MatrixCompression, false);
	p.par.set_entry(BCP_lp_par::LpVerb_NodeTime, false);
	p.par.set_entry(BCP_lp_par::LpVerb_PresolvePositions, false);
	p.par.set_entry(BCP_lp_par::LpVerb_PresolveResult, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ProcessedNodeIndex, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ReportCutGenTimeout, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ReportVarGenTimeout, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ReportLocalCutPoolSize, false);
	p.par.set_entry(BCP_lp_par::LpVerb_ReportLocalVarPoolSize, false);
	p.par.set_entry(BCP_lp_par::LpVerb_RepricingResult, false);
	p.par.set_entry(BCP_lp_par::LpVerb_RowEffectivenessCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_VarTightening, false);
	p.par.set_entry(BCP_lp_par::LpVerb_StrongBranchPositions, false);
	p.par.set_entry(BCP_lp_par::LpVerb_StrongBranchResult, false);
	p.par.set_entry(BCP_lp_par::LpVerb_GeneratedCutCount, false);
	p.par.set_entry(BCP_lp_par::LpVerb_GeneratedVarCount, false);
    
   p.par.set_entry(BCP_lp_par::DoReducedCostFixingAtZero ,false);
   p.par.set_entry(BCP_lp_par::DoReducedCostFixingAtAnything, false);
   p.par.set_entry(BCP_lp_par::ReportWhenDefaultIsExecuted, false);
   p.par.set_entry(BCP_lp_par::MessagePassingIsSerial, true);
   p.par.set_entry(BCP_lp_par::MaxCutsAddedPerIteration, 10000);
   p.par.set_entry(BCP_lp_par::MaxPresolveIter, 100 );
   p.par.set_entry(BCP_lp_par::StrongBranchNum, 100 );
   p.par.set_entry(BCP_lp_par::IneffectiveConstraints, BCP_IneffConstr_None);
   p.par.set_entry(BCP_lp_par::WarmstartInfo, BCP_WarmstartNone); //BCP_WarmstartParent;
   p.par.set_entry(BCP_lp_par::DeletedRowToCompress_Frac, 0.0); //BCP_WarmstartParent;
   p.par.set_entry(BCP_lp_par::DeletedRowToCompress_Min, 0); //BCP_WarmstartParent;
   p.par.set_entry(BCP_lp_par::SendFathomedNodeDesc, true);



   return lp;
}

//#############################################################################
/*
void MCND_read_parameters(MCND_tm& tm, const char * paramfile)
{
  // tm.tm_par.read_from_file(paramfile);
}*/

//#############################################################################
