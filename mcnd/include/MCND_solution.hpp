// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef _MCND_SOLUTION_H
#define _MCND_SOLUTION_H
#include <iostream>

#include "BCP_vector.hpp"
#include "BCP_solution.hpp"

class BCP_buffer;
class BCP_string;
class MCND_problem;

class MCND_solution : public BCP_solution {
public:
  double cost;
  int size;
  double* xy;
  
public:
   MCND_solution():cost(1e30), size(0), xy(0){}
		
   MCND_solution(const int sz) :
		size(sz),cost(1e30){
		xy =new double[sz];
		for(int i= sz;i--;)xy[i]=0;}
	
	
    MCND_solution(const int sz, const double c, double* sol) :
		size(sz),
		cost(c), 
		xy(sol) {}
		
  ~MCND_solution() {if(size) delete[] xy;}

  virtual double objective_value() const { return cost; }


  MCND_solution& operator=(const MCND_solution& sol);


  void resize(int sz){
	 if(size == 0){
		 size = sz;
		 xy = new double[sz];
	 }else if(size > 0 && (size != sz)){
		 delete [] xy;
		 size = sz;
		 xy = new double[sz];
	 }
  }
  
  BCP_buffer& pack(BCP_buffer& buf) const;
  BCP_buffer& unpack(BCP_buffer& buf);

  
};

#endif
