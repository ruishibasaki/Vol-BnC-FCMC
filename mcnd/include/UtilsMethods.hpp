//
//  UtilsMethods.hpp
//  
//
//  Created by Rui Shibasaki on 30/11/2018.
//
//

#ifndef UtilsMethods_hpp
#define UtilsMethods_hpp

#include <stdio.h>
#include "MCND_data.hpp"
#include "Structures.hpp"


//----------------------------------------------------------------------------------
// Read/Translate methods
//----------------------------------------------------------------------------------
void MCND_read_data(std::string fname, Data & data);
bool getScalar(int id,const std::vector<PairF> & mapset, double & gam1, double & gam2);
//----------------------------------------------------------------------------------
// Flow methods
//----------------------------------------------------------------------------------
double dijkstra(int source, int sink, const std::vector<int> &grid, std::vector<Pair1> & preced, const std::vector<double>& cij);

//Auxiliary
void remove_i(int p, std::list<int> &N);
Pair1 arg_min(const std::list<int> & nodes, const std::vector<double> & costs);


//----------------------------------------------------------------------------------
// Knapsack dp
//----------------------------------------------------------------------------------

int MaxKnapsack(int n, int W, const int * wt, const  int *val, int * kws);

//----------------------------------------------------------------------------------
// Math methods
//----------------------------------------------------------------------------------

unsigned int min(unsigned int a, unsigned int b);
double fmin(double a, double b);
double max(double a, double b);
int max(int a, int b);

int factorial(unsigned int n, unsigned int stop);
int combination(unsigned int n, unsigned int p);
void  setBit(unsigned int *A,  int iA , int bit);
void  clearBit(unsigned int *A,   int iA , int bit);
int testBit(unsigned int *A,   int iA , int bit);
int countSetBits(unsigned int n);

#endif /* UtilsMethods_hpp */
