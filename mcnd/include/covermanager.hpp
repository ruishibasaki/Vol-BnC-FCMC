//
//  covermanager.hpp
//  
//
//  Created by Rui Shibasaki on 26/07/2019.
//

#ifndef covermanager_hpp
#define covermanager_hpp

#include <stdio.h>
#include "covercollection.hpp"
#include "cutsetmanager.hpp"
#include "lift.hpp"
#include "MCND_cut.hpp"


class CoverManager: virtual public CutManager {

public:
    
    const Data * data;
    const double * colub; 
    const int * arc_map;
    double * capabl; 
    Lift coverlift;
    CoverCollection covers;
    int num_actv;
    int lim_to_remv;
    int gend;
    int ttgend_cover; //stat can be removed
    int ttgend_card; //stat can be removed
    std::deque<Cover*> purgbl;
    //-------------------------------------------------------------------------------------------
    //  initializing methods
    //-------------------------------------------------------------------------------------------
    
    inline CoverManager(): data(0), arc_map(0), colub(0){ num_actv = lim_to_remv =0; 
    													 ttgend_card = ttgend_cover =0; //stat can be removed
    													 }
    inline ~CoverManager(){ delete [] capabl;}

    inline void set_arc_map(const int * map){ arc_map = map;}
    void initialize(const Data * d, int lim);
    int reset_and_map_collection(int fsize, const double* topo, double * dual, int * actvS, int & csize, bool recheck_collct);
    void clean_collection();
    //-------------------------------------------------------------------------------------------
    //  main methods
    //-------------------------------------------------------------------------------------------
    
    int cover_generation_main(const double * ystar, const double * y,const CutSetCollection * sets, int curr_id, int max );
    int cover_generation(int ss_size, const int * SS_arcs, double uss, double dss,
                         const double * ystar, const double * y, int curr_id );
    Cover * make_cover(double& delta, const std::deque<Trio1> & ss_, const double * ystar, std::deque<Pair2>& lift_down, std::deque<Pair2>& cover, int id_vi  );
	
	
    //-------------------------------------------------------------------------------------------
    //  minimal Cover methods
    //-------------------------------------------------------------------------------------------

    void minimize_cover(double& delta, std::deque<Trio1>& candidates, std::deque<Pair2>& cover, const double * ystar, std::deque<Pair2>& lift_down);
    
    void form_c0(std::deque<Pair2> & lift_up, std::deque<Trio1> & ss_,
                 const double *ystar, double dss, double uss);
    
    void form_c1(std::deque<Pair2> & lift_down, std::deque<Trio1> & ss_,
                 const double *ystar,double & delta, double dss, double uss);
    
    double cutset_preprocess(int sz, double dss, const int * ss_,  std::deque<Trio1>& ss_deque, std::deque<Pair2> & lift_down,
                             const double *y, const double *ystar);
    
    void restrict_cutset(std::deque<Pair2> & lift_down, std::deque<Pair2> & lift_up, std::deque<Trio1> & ss_, const double *ystar,double & delta, double dss, double uss);
    
    
    //-------------------------------------------------------------------------------------------
    //  min cardinality methods
    //-------------------------------------------------------------------------------------------
    
    int mincard_generation_main(const double * ystar, const double * y, const CutSetCollection * sets, int curr_id, int max);
    int mincard_generation(int ss_size, const int * SS_arcs, int ss_ksize, const int * SS_comm, double uss, double dss,
								 const double * ystar, const double * y, int id_vi);
	int make_cardcs( std::deque<Pair2> & ss_, double delta );
	double card_cutset_preprocess(int sz, const int * ss_, int szk, const int * ss_k,  std::deque<Pair2>& ss_deque,
                                const double *y, const double *ystar, double dss);
    //-------------------------------------------------------------------------------------------
    //  auxiliary methods
    //-------------------------------------------------------------------------------------------
    
    double checkViol(const Cover * c, const double *y);
    int add_external_cover(const std::deque<Pair2>& c, int id, int maxNumrows_);
    void reposition_covers(int num_covers);
    //-------------------------------------------------------------------------------------------
    //  Multipliers methods
    //-------------------------------------------------------------------------------------------
    
    double set_new_mult_pos(double *rc, std::vector<Trio1>& ws, const double * dual,
                        const std::deque<int>& con_arcs, const std::vector<int> & con_arcs_map,
                        const std::vector<Pair2>&  con_arcs_wnid);
    double update_dual_pos( std::vector<Trio1>& ws, double * dual, double * h);
    double set_new_mult_neg(double *rc, std::vector<Trio1>& ws,  double * dual, std::vector<Pair2> & con_arcs_map,
                      std::deque<Pair2>& con_arcs, const std::vector<Pair2>& con_arcs_wnid,
                     const std::vector<const Cover *>& addrs ,const double *fk);
    double update_rc_neg(double dimsh, int nvi, const Cover * vi, std::vector<Trio1>& ws,
                         std::vector<Pair2> & con_arcs_map, double *rc);
    double update_con_arcs(std::deque<Pair2>& con_arcs, const double* fk, const double *rc);
    
    //-------------------------------------------------------------------------------------------
    //  Volume Integration methods
    //-------------------------------------------------------------------------------------------
    
    int compute_cover_rc(const double * dual, const int* actvS, int actvSSz, double * rc, double & B0);
    int compute_cover_sg(const double * x, const int * actvS, int actvSSz,  double * v);
    void add_cover_vi(int added, int * actvS, int & actvSSz,double * h, double * dstar, double * dual, double * dual_lb, double * dual_ub );
    void make_inactive(int index, const int* actvS, double* v);
    double recompute_mult_pos( double * dual, double * h,  double *rc, 
                               const double *xy, const int * actvS);
    double recompute_mult_neg(double * dual, double * rc, const double *fk,
                              const double *xy, const int * actvS, int actvSSz);
    double arc_dg_imp(int arc, const double * xy, const double * h, const int * actvS, int actvSSz);
};

#endif /* covermanager_hpp */
