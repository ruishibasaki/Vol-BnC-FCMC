// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

// #include <cfloat>
#include <cstdio>
#include <fstream>

#include "BCP_math.hpp"
#include "BCP_string.hpp"
#include "BCP_buffer.hpp"
#include "BCP_vector.hpp"

#include "MCND_data.hpp"
#include "MCND_solution.hpp"



//#############################################################################

BCP_buffer&
MCND_solution::pack(BCP_buffer& buf) const {
		//std::cout<<"packing solution: cost "<<cost<<" size: "<<size<<std::endl;

  buf.pack(cost).pack(size).pack(onlyvalue);
    for(int i=size;i--;){
		buf.pack(xy[i]);
	}
	
	return buf;
}

//#############################################################################

BCP_buffer&
MCND_solution::unpack(BCP_buffer& buf) {
	//std::cout<<"unpack sol"<<std::endl;
  bool alloc=false;
  if(size==0)alloc=true;
  
  buf.unpack(cost).unpack(size).unpack(onlyvalue);
  
  if(alloc)xy=new double[size];
  for(int i=size;i--;)
	buf.unpack(xy[i]);
	
	return buf;
}

//#############################################################################

void 
MCND_solution::copy(const MCND_solution& sol, int sz){
  bool alloc=false;
  if(size==0)alloc=true;
  
  sz = sol.size < sz? sol.size: sz;
  cost = sol.cost;
  size = sz;
  onlyvalue = sol.onlyvalue;
   
  if(alloc)xy=new double[sz];
    for(int i=sz;i--;)
		xy[i] = sol.xy[i];
		
}


//#############################################################################

MCND_solution&
MCND_solution::operator=(const MCND_solution& sol) {
	//std::cout<<"copy solution"<<std::endl;
  bool alloc=false;
  if(size==0)alloc=true;
  
  cost = sol.cost;
  size = sol.size;  
  onlyvalue = sol.onlyvalue;
  
  if(alloc)xy=new double[size];
    for(int i=size;i--;)
		xy[i] = sol.xy[i];
		
  return *this;
}

