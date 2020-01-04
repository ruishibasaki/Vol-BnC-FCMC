#ifndef _KSOLVER_H
#define _KSOLVER_H

#include <vector>
#include <deque>
#include "Structures.hpp"
#include "UtilsMethods.hpp"
#include "covercollection.hpp"
//#include "knapsolver.hpp"

class Lift{

public:
	//problem data
    //KnapSolver knapslvr;

	const Data * data;
	int nnodes;
	int ndemands;
	int narcs;
    int dss;
    //lift data
    int *p;int *w;
    int * states;
    int max_n; //max n on (n*C) states
    int max_dbar; //max C on (n*C) states
    int max_up_dbar;
    int min_dbar;
    int max_u;
    int luc_tt_u;
    int sz_states;

    int dim; //cover current dimention
    
    //Main methods
    CoverL* lift_cover(std::deque<Pair2>& luc, std::deque<Pair2> & lift_down , std::deque<Pair2> & lift_up , int dbar, int rhs, int dss_);
    int solve( int u, int  dbar,  int strt, int end, bool down);
    int dpsolver(int W, int strt, int end, bool down);

    //Down Lifting methods
    void lift_vars_down(CoverL* vi,  std::deque<Pair2> & lift_down ,std::deque<Pair2>& luc, int & dbar, int & ttgamcul, int & ttgamc1, int rhs, bool down_first);
    int down_lifting(int strt, int end,int arc , CoverL* vi , std::deque<Pair2>& luc, int & dbar, int & ttgamcul, int & ttgamc1, int rhs);
    
    //Up Lifting methods
    void lift_vars_up(CoverL* vi,  std::deque<Pair2> & lift_up ,std::deque<Pair2>& luc, int& dbar, int& ttgamcul, int& ttgamc1, int rhs, bool down_first);
    int up_lifting(int strt, int end,int arc , CoverL* vi , std::deque<Pair2>& luc, int & dbar, int & ttgamcul, int & ttgamc1, int rhs);
    
    //setters and initializing methods
    int reset( const std::deque<Pair2> & lift_down, int dbar);
	void set_data(const Data * d);
    int build(const std::deque<Pair2>& luc,const std::deque<Pair2>& lift_up,const std::deque<Pair2>& lift_down);
    void add_lifted_var(bool down, CoverL* vi, std::deque<Pair2>& luc, int dbar, int u, int gam, int arc, int & ttgamcul, int & ttgamc1);

    //constructors/destructors and auxiliars
    inline Lift(): p(0), w(0), states(0){max_n = max_dbar = max_up_dbar = min_dbar = max_u = luc_tt_u = sz_states=0;
        dss=dim=0;
        nnodes = ndemands =narcs=0;
    }
    inline void clear(){if(sz_states) delete [] states; sz_states=0;}
    inline ~Lift(){delete [] p; delete [] w; if(sz_states) delete [] states;}
    void print(const std::deque<Pair2>& luc, int dbar);
};

#endif 

