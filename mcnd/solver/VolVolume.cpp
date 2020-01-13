// $Id: VolVolume.cpp 281 2011-01-05 23:58:31Z lou $
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>


#include "VolVolume.hpp"

//#############################################################################
/// Usage: v=w; where w is a VOL_dvector
VOL_dvector&
VOL_dvector::operator=(const VOL_dvector& w) {
    if (this == &w)
        return *this;
    delete[] v;
    const int wsz = w.size();
    if (wsz == 0) {
        v = 0;
        sz = 0;
    } else {
        v = new double[sz = wsz];
        for (int i = sz - 1; i >= 0; --i)
            v[i] = w[i];
    }
    return *this;
}
//-----------------------------------------------------------------------

/// Usage v=w; where w is a double. It copies w in every entry of v
VOL_dvector&
VOL_dvector::operator=(const double w) {
    for (int i = sz - 1; i >= 0; --i)
        v[i] = w;
    return *this;
}


//#############################################################################
/// Usage: v=w; where w is a VOL_ivector
VOL_ivector&
VOL_ivector::operator=(const VOL_ivector& w) {
    if (this == &w) return *this;
    delete[] v;
    const int wsz = w.size();
    if (wsz == 0) {
        v = 0;
        sz = 0;
    } else {
        v = new int[sz = wsz];
        for (int i = sz - 1; i >= 0; --i)
            v[i] = w[i];
    }
    return *this;
}

//-----------------------------------------------------------------------

/// Usage v=w; where w is an int. It copies w in every entry of v
VOL_ivector&
VOL_ivector::operator=(const int w) {
    for (int i = sz - 1; i >= 0; --i)
        v[i] = w;
    return *this;
}

//############################################################################
/// find maximum absolute value of the primal violations
void
VOL_primal::find_max_viol(const VOL_dvector& dual_lb, 
                          const VOL_dvector& dual_ub, int actvSSz)
{
    const int nc = actvSSz;
    viol = 0;
    for ( int i = 0; i < nc; ++i ) {
        if ( (v[i] > 0.0 && dual_ub[i] != 0.0) ||
            (v[i] < 0.0 && dual_lb[i] != 0.0) )
            viol = VolMax(viol, VolAbs(v[i]));
    }
}

//############################################################################
/// Dual step. It takes a step in the direction v
// lcost is a member of VOL_dual
double
VOL_dual::step(const double target, const double lambda,
               const VOL_dvector& dual_lb, const VOL_dvector& dual_ub,
               const VOL_dvector& v,  int actvSSz) {
    const int nc = actvSSz;
    int i;
    
    double viol = 0.0;
    for (i = 0; i < nc; ++i) {
        if (( v[i] > 0.0 && u[i] < dual_ub[i] ) ||
            ( v[i] < 0.0 && u[i] > dual_lb[i] )) {
            viol += v[i] * v[i];
        }
    }
    
    const double stp = viol == 0.0 ? 0.0 : (target - lcost) / viol * lambda;
    //std::cout<<"step "<<stp<<" trgt: "<<target<<" rho: "<<lambda<<" viol: "<<viol<<std::endl;
    for (i = 0; i < nc; ++i) {
        //if(index>=2000)std::cout<<"stp: "<<stp<<" "<<v[index]<<std::endl;
        if (( v[i] > 0.0 && u[i] < dual_ub[i] ) ||
            ( v[i] < 0.0 && u[i] > dual_lb[i] )) {
            u[i] += stp * v[i];
            if (u[i] < dual_lb[i])
                u[i] = dual_lb[i];
            else if (u[i] > dual_ub[i])
                u[i] = dual_ub[i];
        }
    }
    return stp;
}

//-----------------------------------------------------------------------

/// ascent = inner product(v, u - last_u)
double
VOL_dual::ascent(const VOL_dvector& v, const VOL_dvector& last_u, int actvSSz ) const
{
    const int nc = actvSSz;
    int i;
    double asc = 0.0;
    for (i = 0; i < nc; ++i){
        asc += v[i] * (u[i] - last_u[i]);
    }
    return asc;
}

//-----------------------------------------------------------------------

/** compute xrc. This is (c - u A) * ( xstar - x ). This is just
 miscellaneous information, it is not used in the algorithm. */
void 
VOL_dual::compute_xrc(const VOL_dvector& xstar, const VOL_dvector& x,
                      const VOL_dvector& rc)
{
    const int nc = x.size();
    xrc = 0;
    for (int i = 0; i < nc; ++i) {
        xrc += rc[i] * (xstar[i] - x[i]);
    }
}

//-----------------------------------------------------------------------

/** Copy only Atcive_Set values of <code>w</code> into the vector. */
void
VOL_dual::copy(const VOL_dual& w, int actvSSz){
    lcost = w.lcost;
    xrc = w.xrc;
    const int wsz = w.u.size();
    if (wsz != u.size()) {
        u.clear();
        u = w.u;
    }else{
        for (int i = actvSSz; i--; ){
            u[i] = w.u[i];
        }
    }
}

//############################################################################
/** Computing inner products. It computes v * ( alpha v + (1-alpha) h),
 v * h, v * v, h * h. Here v is the subgradient direction, and h is
 the conjugate direction. */
VOL_vh::VOL_vh(const double alpha,
               const VOL_dvector& dual_lb, const VOL_dvector& dual_ub,
               const VOL_dvector& v, const VOL_dvector& vstar,
               const VOL_dvector& u, int actvSSz ) :
hh(0), norm(0), vh(0), asc(0)
{
    int i;
    const int nc = actvSSz;
    double vv;
    
    for (i = 0; i < nc; ++i) {
        const double vi = v[i];
        const double vsi = vstar[i];
        vv = alpha * vi + (1.0 - alpha) * vsi;
        if (u[i] == 0.0  &&  dual_lb[i] == 0.0  &&  vv <= 0.0)
            continue;
        if (u[i] == 0.0  &&  dual_ub[i] == 0.0  &&  vv >= 0.0)
            continue;
        asc  += vi * vv;
        vh   += vi * vsi;
        norm += vi * vi;
        hh   += vsi * vsi;
    }
}

//############################################################################
/** Computes indicators for printing. They are v2=vstar * vstar, asc= v*v,
 vu= vstar * u, vabs = sum( abs(vstar[i]))/m, v2= sum( vstar[i]^2) / m.
 */
VOL_indc::VOL_indc(const VOL_dvector& dual_lb, const VOL_dvector& dual_ub,
                   const VOL_primal& primal, const VOL_primal& pstar,
                   const VOL_dual& dual, int actvSSz ) {
    
    v2 = vu = vabs = asc = 0.0;
    const VOL_dvector v = primal.v;
    const VOL_dvector vstar = pstar.v;
    const VOL_dvector u = dual.u;
    int i;
    const int nc = actvSSz;
    
    for (i = 0; i < nc; ++i) {
        if (u[i] == 0.0  &&  dual_lb[i] == 0.0  &&  vstar[i] <= 0.0)
            continue;
        if (u[i] == 0.0  &&  dual_ub[i] == 0.0  &&  vstar[i] >= 0.0)
            continue;
        v2   += vstar[i] * vstar[i];
        asc  += v[i] * v[i];
        vu   -= vstar[i] * u[i];
        vabs += VolAbs(vstar[i]);
    }
    
    v2 = sqrt(v2) / nc;
    vabs /= nc;
}

//############################################################################

// reading parameters that control the algorithm
void
VOL_problem::read_params(const char* filename) 
{
    char s[100];
    FILE* infile = fopen(filename, "r");
    if (!infile) {
        printf("Failure to open file: %s\n", filename);
        abort();
    }
    while (fgets(s, 100, infile)) {
        const int len = (int)strlen(s) - 1;
        if (s[len] == '\n')
            s[len] = 0;
        std::string ss(s);
        
        if (ss.find("printflag") == 0) {
            int i = (int)ss.find("=");
            parm.printflag = atoi(&s[i+1]);
            
        } else if (ss.find("printinvl") == 0) {
            int i = (int)ss.find("=");
            parm.printinvl = atoi(&s[i+1]);
            
        } else if (ss.find("maxsgriters") == 0) {
            int i = (int)ss.find("=");
            parm.maxsgriters = atoi(&s[i+1]);
            
        } else if (ss.find("heurinvl") == 0) {
            int i = (int)ss.find("=");
            parm.heurinvl = atoi(&s[i+1]);
            
        } else if (ss.find("greentestinvl") == 0) {
            int i =(int) ss.find("=");
            parm.greentestinvl = atoi(&s[i+1]);
            
        } else if (ss.find("yellowtestinvl") == 0) {
            int i =(int) ss.find("=");
            parm.yellowtestinvl = atoi(&s[i+1]);
            
        } else if (ss.find("redtestinvl") == 0) {
            int i = (int)ss.find("=");
            parm.redtestinvl = atoi(&s[i+1]);
            
        } else if (ss.find("lambdainit") == 0) {
            int i =(int) ss.find("=");
            parm.lambdainit = atof(&s[i+1]);
            
        } else if (ss.find("alphainit") == 0) {
            int i = (int)ss.find("=");
            parm.alphainit = atof(&s[i+1]);
            
        } else if (ss.find("alphamin") == 0) {
            int i =(int) ss.find("=");
            parm.alphamin = atof(&s[i+1]);
            
        } else if (ss.find("alphafactor") == 0) {
            int i = (int)ss.find("=");
            parm.alphafactor = atof(&s[i+1]);
            
        } else if (ss.find("alphaint") == 0) {
            int i = (int)ss.find("=");
            parm.alphaint = atoi(&s[i+1]);
            
        } else if (ss.find("primal_abs_precision") == 0) {
            int i = (int)ss.find("=");
            parm.primal_abs_precision = atof(&s[i+1]);
            
        } else if (ss.find("gap_abs_precision") == 0) {
            int i = (int)ss.find("=");
            parm.gap_abs_precision = atof(&s[i+1]);
            
        } else if (ss.find("gap_rel_precision") == 0) {
            int i =(int) ss.find("=");
            parm.gap_rel_precision = atof(&s[i+1]);
            
        } else if (ss.find("ascent_check_invl") == 0) {
            int i = (int)ss.find("=");
            parm.ascent_check_invl = atoi(&s[i+1]);
            
        } else if (ss.find("minimum_rel_ascent") == 0) {
            int i =(int) ss.find("=");
            parm.minimum_rel_ascent = atoi(&s[i+1]);
            
        } else if (ss.find("granularity") == 0) {
            int i = (int)ss.find("=");
            parm.granularity = atof(&s[i+1]);
            
        }else if (ss.find("maxtime") == 0) {
            int i = (int)ss.find("=");
            parm.maxtime = atof(&s[i+1]);
        }
    }
    fclose(infile);
}

//#############################################################################

void
VOL_problem::set_default_parm()
{
    parm.lambdainit = 0.1;
    parm.alphainit = 0.01;
    parm.alphamin = 0.001;
    parm.alphafactor = 0.5;
    parm.ubinit = DBL_MAX;
    parm.primal_abs_precision = 0.02;
    parm.gap_abs_precision = 0.0;
    parm.gap_rel_precision = 0.001;
    parm.granularity = 0.0;
    parm.minimum_rel_ascent = 0.0001;
    parm.ascent_first_check = 500;
    parm.ascent_check_invl = 100;
    parm.maxsgriters = 2000;
    parm.printflag = 3;
    parm.printinvl = 50;
    parm.heurinvl = 100000000;
    parm.greentestinvl = 1;
    parm.yellowtestinvl = 2;
    parm.redtestinvl = 10;
    parm.alphaint = 80;
    parm.temp_dualfile = 0;
    parm.maxtime=3600;
}

//#############################################################################

VOL_problem::VOL_problem() : 
alpha_(-1),
lambda_(-1),
iter_(0),
value(-1),
psize(-1),
dsize(-1)
{
    set_default_parm();
}

//-----------------------------------------------------------------------

VOL_problem::VOL_problem(const char *filename) : 
alpha_(-1),
lambda_(-1),
iter_(0),
value(-1),
psize(-1),
dsize(-1)
{
    set_default_parm();
    read_params(filename);
}

//######################################################################

VOL_problem::~VOL_problem()
{
    delete[] parm.temp_dualfile;
}

//######################################################################
/// print information about the current iteration
void
VOL_problem::print_info(const int iter,
                        const VOL_primal& primal, const VOL_primal& pstar,
                        const VOL_dual& dual)
{
    VOL_indc indc(dual_lb, dual_ub, primal, pstar, dual,active_size );
    printf("%i. L=%f P=%f vu=%f infeas=%f\n asc=%f vmax=%f P-vu=%f xrc =%f\n",
           iter, dual.lcost, pstar.value, indc.vu, indc.v2, indc.asc,
           pstar.viol, pstar.value - indc.vu, dual.xrc);
}

//######################################################################
/// this is the Volume Algorithm
int
VOL_problem::solve(VOL_user_hooks& hooks, const bool use_preset_dual) 
{
    clock_t t_u;
    clock_t t_tt;
    struct tms buff;
    times( &buff );
    t_u = buff.tms_utime;
    double timespent = 0;
    ttime_solve=0;
    ttime_search=0;
    if (initialize(use_preset_dual) < 0) // initialize several parameters
        return -1;

    double best_ub = parm.ubinit;      // upper bound
    int retval = 0;
    double step=0;
    original_sz = active_size;

    VOL_dvector rc(psize); // reduced costs
    VOL_dual dual(dsize); // dual vector
    dual.u = dsol;
    
    VOL_primal primal(psize, dsize);  // primal vector
    VOL_primal pstar(psize, dsize);
    VOL_dual dstar(dual);
    VOL_dual dlast(dual);

    retval = hooks.compute_rc(dual.u, rc, active_size); // compute reduced costs
    if (retval < 0)  return -1;
    // solve relaxed problem
    
    times( &buff );
    t_tt = buff.tms_utime;
    retval = hooks.solve_subproblem(pstar.x,dual.u, rc, dual.lcost, primal.x, primal.value);
    if (retval < 0)  return -1;
    retval = hooks.additional_settings(0, dual.lcost, dual.u, rc, pstar.v, primal.x, pstar.x, dstar.u, active_size);
    if (retval < 0)  return -1;
    retval = hooks.resolve_subproblem(dual.u, rc, dual.lcost, primal.x, primal.value);
    if (retval < 0)  return -1;
    times( &buff );
    ttime_solve += ( double( buff.tms_utime - t_tt ) ) / double( CLK_TCK );
    retval = hooks.addVI(iter_,dual.lcost, pstar.x, primal.x, dual_lb, dual_ub, primal.v, pstar.v, active_size );
    if (retval < 0)  return -1;
    retval = hooks.compute_sg(primal.x, active_size, primal.v);
    if (retval < 0)  return -1;
    std::cout<<"first vole iter: "<<dual.lcost<<std::endl;
    
    // set target for the lagrangian value
    double target = readjust_target(-DBL_MAX/2, dual.lcost);
    // find primal violation
    //primal.find_max_viol(dual_lb, dual_ub); // this may be left out for speed
    
    pstar = primal; // set pstar=primal
    //pstar.find_max_viol(dual_lb, dual_ub); // set violation of pstar
    
    //dual.compute_xrc(pstar.x, primal.x, rc); // compute xrc
    
    //
    dstar.copy(dual, active_size); // dstar is the best dual solution so far
    dlast.copy(dual, active_size); // set dlast=dual
    
    iter_ = 0;
    if (parm.printflag)
        print_info(iter_, primal, pstar, dual);
    
    VOL_swing swing;
    VOL_alpha_factor alpha_factor;
    
    double * lcost_sequence = new double[parm.ascent_check_invl];
    
    const int ascent_first_check = VolMax(parm.ascent_first_check,
                                          parm.ascent_check_invl);
    
    //std::cout<<"maxiteers: "<<parm.maxsgriters<<" acts:"<<active_size<<std::endl;
    for (iter_ = 1; iter_ <= parm.maxsgriters; ++iter_) {  // main iteration

        retval = hooks.removeVI( active_size, pstar.v, dstar.u, dual.u);
        if (retval < 0)  break;

        dlast.copy(dual, active_size);
        // take a dual step
        step = dual.step(target, lambda_, dual_lb, dual_ub, pstar.v, active_size);
        // compute reduced costs
        retval = hooks.compute_rc(dual.u, rc, active_size);
        if (retval < 0)  break;
        // solve relaxed problem

        times( &buff );
        t_tt = buff.tms_utime;
        retval = hooks.solve_subproblem(pstar.x,dual.u, rc, dual.lcost,primal.x, primal.value);
        if (retval < 0)  break;
        retval = hooks.additional_settings(iter_, dual.lcost, dual.u, rc, pstar.v, primal.x, pstar.x, dstar.u, active_size);
        if (retval < 0)  break;
        retval = hooks.resolve_subproblem(dual.u, rc, dual.lcost, primal.x, primal.value);
        if (retval < 0)  break;

        times( &buff );
        ttime_solve += ( double( buff.tms_utime - t_tt ) ) / double( CLK_TCK );
        
        times( &buff );
        t_tt = buff.tms_utime;
        retval = hooks.addVI(iter_,dual.lcost, pstar.x, primal.x, dual_lb, dual_ub, primal.v, pstar.v, active_size );
        if (retval < 0)  break;
        times( &buff );
        ttime_search += ( double( buff.tms_utime - t_tt ) ) / double( CLK_TCK );
        retval = hooks.compute_sg(primal.x, active_size, primal.v);
        if (retval < 0)  break;
        
        //primal.find_max_viol(dual_lb, dual_ub); // this may be left out for speed
        //dual.compute_xrc(pstar.x, primal.x, rc); // compute xrc
        if (dual.lcost > dstar.lcost) {
            //std::cout<<"aug: "<<dual.lcost-dstar.lcost<<std::endl;
            dstar.copy(dual, active_size);
            //std::cout<<"L*: "<<dstar.lcost<<std::endl;
        }

        // compute inner product between the new subgradient and the
        // last direction. This to decide among green, yellow, red
        const double ascent = dual.ascent(primal.v, dlast.u,active_size);
        // green, yellow, red
        swing.cond(dlast, dual.lcost, ascent, iter_);
        //std::cout<<iter_<<" swing: "<<swing.lastgreeniter<<" , "<<swing.lastyellowiter<<" , "<<swing.lastrediter<<std::endl;;
        // change lambda if needed
        lambda_ *= swing.lfactor(parm, lambda_, iter_);
        

        if (iter_ % parm.alphaint == 0) { // change alpha if needed
            const double fact = alpha_factor.factor(parm, dstar.lcost, alpha_);
            if (fact != 1.0 && (parm.printflag & 2)) {
                printf(" ------------decreasing alpha to %f\n", alpha_*fact);
            }
            alpha_ *= fact;
        }
        
        // convex combination with new primal vector
        pstar.cc(power_heur(primal, pstar, dual), primal,active_size);
        pstar.find_max_viol(dual_lb, dual_ub, active_size); // find maximum violation of pstar
        if (swing.rd)
            dual.copy(dstar, active_size);
        // if there is no improvement reset dual=dstar
        // check if target should be updated
        target = readjust_target(target, dstar.lcost);

        if ((iter_ % parm.printinvl == 0) && parm.printflag) { // printing iteration information
            print_info(iter_, primal, pstar, dual);
            swing.print();
        }

        if (iter_ % parm.heurinvl == 0) { // run primal heuristic
            double ub = DBL_MAX;
            retval = hooks.heuristics(*this, pstar.x, ub);
            if (retval < 0)  break;
            if (ub < best_ub)
                best_ub = ub;
        }

        // test terminating criteria
        const bool primal_feas = (pstar.viol < parm.primal_abs_precision);
        //const double gap = VolAbs(pstar.value - dstar.lcost);
        const double gap = pstar.value - dstar.lcost;
        const bool small_gap = VolAbs(dstar.lcost) < 0.0001 ?
        (gap < parm.gap_abs_precision) :
        ( (gap < parm.gap_abs_precision) ||
         (gap/VolAbs(dstar.lcost) < parm.gap_rel_precision) );

        
        // test optimality
        if (primal_feas && small_gap){
            if (parm.printflag) printf(" small lp gap \n");
            break;
        }
        
        // test proving integer optimality
        if (best_ub - dstar.lcost < parm.granularity){
            if (parm.printflag) printf(" small ip gap \n");
            break;
        }

        // test for non-improvement
        const int k = iter_ % parm.ascent_check_invl;
        if (iter_ > ascent_first_check) {
            if (dstar.lcost - lcost_sequence[k] <
                VolAbs(lcost_sequence[k]) * parm.minimum_rel_ascent){
                if (parm.printflag) printf(" small improvement \n");
                break;
            }
        }
        

        lcost_sequence[k] = dstar.lcost;
        

        times( &buff );
        timespent = ( double( buff.tms_utime - t_u ) ) / double( CLK_TCK );
        if (timespent > parm.maxtime) {
            retval=0;
            break;
        }


    }
    delete[] lcost_sequence;
    
    if (parm.printflag)
        print_info(iter_, primal, pstar, dual);
    // set solution to return
    value = dstar.lcost;
    psol = pstar.x;
    dsol = dstar.u;
    viol = pstar.v;
    
    return retval;
}

//-----------------------------------------------------------------------

/// A function to initialize a few variables
int
VOL_problem::initialize(const bool use_preset_dual) {
    // setting bounds for dual variables
    if (dual_lb.size() > 0) {
        if (dual_lb.size() != dsize) {
            printf("size inconsistent (dual_lb)\n");
            return -1;
        }
    } else {
        // fill it with -infinity
        dual_lb.allocate(dsize);
        dual_lb = - DBL_MAX;
    }
    
    if (dual_ub.size() > 0) {
        if (dual_ub.size() != dsize) {
            printf("size inconsistent (dual_ub)\n");
            return -1;
        }
    } else {
        // fill it with infinity
        dual_ub.allocate(dsize);
        dual_ub = DBL_MAX;
    }
    // setting initial values for parameters
    alpha_       = parm.alphainit;
    lambda_      = parm.lambdainit;
    // check if there is an initial dual solution
    if (use_preset_dual) {
        if (dsol.size() != dsize) {
            printf("size inconsistent (dsol)\n");
            return -1;
        }
    } else {
        dsol.clear();
        dsol.allocate(dsize);
        dsol = 0.0;
    }
    return 0;
}

//-----------------------------------------------------------------------

/// Here we increase the target once we get within 5% of it
double
VOL_problem::readjust_target(const double oldtarget, const double lcost) const
{
    double target = oldtarget;
    if (lcost >= target - VolAbs(target) * 0.05) {
        if (VolAbs(lcost) < 10.0) {
            target = 10.0;
        }else if(lcost >= target){
            target = VolMax(lcost + 0.05 * VolAbs(lcost), lcost + 0.75*(lcost - VolAbs(target)));
            //std::cout<<"  >=  **** readjusting target!!! new target = "<<target<<" , old: "<<oldtarget<<" lcost:  "<<lcost<<std::endl;
        } else {
            target += 0.025 * VolAbs(target);
            target = VolMax(target, lcost + 0.05 * VolAbs(lcost));
            //std::cout<<"  <  **** readjusting target!!! new target = "<<target<<" , old: "<<oldtarget<<std::endl;
        }
        if (target != oldtarget && (parm.printflag & 2)) {
            printf("     **** readjusting target!!! new target = %f *****\n",
                   target);
        }
    }
    return target;
}

//-----------------------------------------------------------------------

/** Here we decide the value of alpha_fb to be used in the convex
 combination. More details of this are in doc.ps
 IN:  alpha, primal, pstar, dual
 OUT: pstar = alpha_fb * pstar + (1 - alpha_fb) * primal
 */
double
VOL_problem::power_heur(const VOL_primal& primal, const VOL_primal& pstar,
                        const VOL_dual& dual) const
{

    const double alpha = alpha_;
    VOL_vh prod(alpha, dual_lb, dual_ub, primal.v, pstar.v, dual.u, active_size);
    
    double a_asc;
    double alpha_fb;
    
    if (prod.norm + prod.hh - 2.0 * prod.vh > 0.0)
        alpha_fb = (prod.hh - prod.vh) / (prod.norm + prod.hh - 2.0 * prod.vh);
    else
        alpha_fb = alpha;
    
    if((prod.norm - prod.vh)!=0) a_asc = (alpha * prod.norm - prod.vh) / (prod.norm - prod.vh);
    else a_asc = alpha_fb;
    if(a_asc<0) a_asc = alpha_fb;
    
    if (alpha_fb > alpha)
        alpha_fb = alpha;
    if (alpha_fb < a_asc)
        alpha_fb = a_asc;
    if (alpha_fb > 1.0)
        alpha_fb = alpha;
    if (alpha_fb < 0.0)
        alpha_fb = alpha / 10.0;
    
    //std::cout<<"a_asc: "<<a_asc<<": "<<(alpha * prod.norm - prod.vh)<<" ; "<<(prod.norm - prod.vh)<<" alpha: "<<alpha_fb<<" alph0: "<<alpha<<std::endl;
    
    return alpha_fb;
}

