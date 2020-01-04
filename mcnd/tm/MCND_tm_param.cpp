// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#include "MCND_tm_param.hpp"
#include "BCP_parameters.hpp"

using std::make_pair;

template <>
void BCP_parameter_set<MCND_tm_par>::create_keyword_list() {
  
   //--------------------------------------------------------------------------
   // StringPar
   keys.push_back(make_pair(BCP_string("MCND_FeasSolFile"),
			    BCP_parameter(BCP_StringPar,
					  FeasSolFile)));
   keys.push_back(make_pair(BCP_string("MCND_InputFile"),
			    BCP_parameter(BCP_StringPar,
					  InputFile)));
   keys.push_back(make_pair(BCP_string("MCND_SolutionFile"),
			    BCP_parameter(BCP_StringPar,
					  SolutionFile)));

}

//#############################################################################

template <>
void BCP_parameter_set<MCND_tm_par>::set_default_entries(){
   
   // StringPar
   set_entry(FeasSolFile, "");
   set_entry(InputFile, "c33.dow");
   set_entry(SolutionFile, "");
   
}
