// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


#ifndef OsiVolSolverInterface_H
#define OsiVolSolverInterface_H

#include <string>
#include <vector>
#include <deque>
#include <list>

#include "VolVolume.hpp"
#include "Structures.hpp"
#include "MCND_data.hpp"
#include "MCND_osidata.hpp"
#include "MCND_solution.hpp"

#include "cutsetmanager.hpp"
#include "covermanager.hpp"
#include "CoinPackedMatrix.hpp"

#include "OsiSolverInterface.hpp"
#include "CoinWarmStartDual.hpp"
#include "WarmStartDual.hpp"
static const double OsiVolInfinity = 1.0e31;

//#############################################################################


/** Vol(ume) Solver Interface
 
 Instantiation of OsiVolSolverInterface for the Volume Algorithm
 */

class OsiVolSolverInterface : virtual public OsiSolverInterface,
virtual public OsiBabSolver, public VOL_user_hooks {
    
    
public:
    //---------------------------------------------------------------------------
    /**@name Solve methods */
    
    /// Solve initial LP relaxation
    virtual void initialSolve();
    
    /// Resolve an LP relaxation after problem modification
    virtual void resolve();
    
    /// Invoke solver's built-in enumeration algorithm
    virtual void branchAndBound() {
        throw CoinError("Sorry, the Volume Algorithm doesn't implement B&B",
                        "branchAndBound", "OsiVolSolverInterface");
    }
    
    virtual bool isAbandoned() const{return (retval==-100)? true : false;}
    
    virtual bool isProvenOptimal() const{
        return ((retval==0) &&
                volprob_.iter() < volprob_.parm.maxsgriters);}
    
    virtual bool isProvenPrimalInfeasible()const{
        if(retval==-1 || volprob_.value<0 ){
            //std::cout<<"VOL: "<<retval<<" / "<< volprob_.value<<std::endl;
            return true;
        }else return false;
    }
    
    virtual bool isProvenDualInfeasible() const{return false;}
    
    virtual bool isIterationLimitReached() const{
        if(volprob_.iter()>= volprob_.parm.maxsgriters)
            return true;
        else return false;}
    
    /// Is the given primal objective limit reached?
    virtual bool isPrimalObjectiveLimitReached() const{
        if(!isProvenPrimalInfeasible())
            return (volprob_.value>volprob_.parm.ubinit)? true : false;
        else return false;}
    
    
    /// Is the given dual objective limit reached?
    virtual bool isDualObjectiveLimitReached() const{return false;}
    //---------------------------------------------------------------------------
    
    
    //---------------------------------------------------------------------------
    /**@name WarmStart related methods */
    //@{
    /*! \brief Get an empty warm start object
     
     This routine returns an empty warm start object. Its purpose is
     to provide a way to give a client a warm start object of the
     appropriate type, which can resized and modified as desired.
     */
    virtual CoinWarmStart *getEmptyWarmStart() const;
    
    /// Get warmstarting information
    virtual CoinWarmStart* getWarmStart() const;
    /** Set warmstarting information. Return true/false depending on whether
     the warmstart information was accepted or not. */
    virtual bool setWarmStart(const CoinWarmStart* warmstart);
    
    virtual void markHotStart();
    /// Optimize starting from the hot start snapshot.
    virtual void solveFromHotStart();
    /// Delete the hot start snapshot.
    virtual void unmarkHotStart();
    
    //@}
    //---------------------------------------------------------------------------
    /**@name Problem information methods
     
     These methods call the solver's query routines to return
     information about the problem referred to by the current object.
     Querying a problem that has no data associated with it result in
     zeros for the number of rows and columns, and NULL pointers from
     the methods that return vectors.
     
     Const pointers returned from any data-query method are valid as
     long as the data is unchanged and the solver is not called.
     */
    //@{
    /**@name Methods related to querying the input data */
    //@{
    virtual int getNumCols() const{return numcols_;};
    virtual int getNumRows() const{return numrows_;};
    virtual int getMaxNumCols() const {return maxNumcols_;};
    virtual int getMaxNumRows() const {return maxNumrows_;};
    virtual int getNumElements() const {return 0;}
    virtual const double * getColLower() const{return collb;};
    virtual const double * getColUpper() const{return colub;};
    virtual const char * getRowSense() const { return 0; }
    virtual const double * getRightHandSide() const { return 0; }
    virtual const double * getRowRange() const { return 0; }
    virtual const double * getRowLower() const {return rowlb;};
    virtual const double * getRowUpper() const {return rowub;};
    virtual const double * getObjCoefficients() const { return 0; }
    virtual double getObjSense() const { return 1; }
    virtual bool isContinuous(int colNumber) const {return true;}
    virtual const CoinPackedMatrix * getMatrixByRow() const{ return 0; }
    virtual const CoinPackedMatrix * getMatrixByCol() const{ return 0; }
    virtual double getInfinity() const { return OsiVolInfinity; }
    //@}
    
    /**@name Methods related to querying the solution */
    //@{
    /// Get pointer to array[getNumCols()] of primal solution vector
    virtual const double * getColSolution() const{//std::cout<<"get primal "<<std::endl;
        return solution;}
    
    virtual const double * getRowPrice() const {std::cout<<"get dual "<<std::endl;
        return dual; }
    
    virtual const double * getReducedCost() const {//std::cout<<"get rc "<<std::endl;
        return rc_; }
    
    virtual const double * getRowActivity() const {//std::cout<<"get ractiv"<<std::endl;
        return lhs_; }
    
    virtual double getObjValue() const {
        if(retval==-1 || volprob_.value<0 ){
            return getInfinity();
        }else return volprob_.value;
    }
    virtual int getIterationCount() const { return volprob_.iter(); }
    
    
    virtual std::vector<double*> getDualRays(int maxNumRays,
                                             bool fullRay = false) const{
        std::vector<double*> ret;
        ret.push_back(dual);
        return ret;
    }
    
    virtual std::vector<double*> getPrimalRays(int maxNumRays) const{
        std::vector<double*> ret;
        ret.push_back(solution);
        return ret;
    }
    
    //@}
    //@}
    
    //---------------------------------------------------------------------------
    
    /**@name Problem modifying methods */
    //@{
    //-------------------------------------------------------------------------
    /**@name Changing bounds on variables and constraints */
    //@{
    /** Set an objective function coefficient */
    virtual void setObjCoeff( int elementIndex, double elementValue ) {
        
    }
    
    using OsiSolverInterface::setColLower ;
    /** Set a single column lower bound<br>
     Use -DBL_MAX for -infinity. */
    virtual void setColLower( int elementIndex, double elementValue ) {
        
    }
    
    using OsiSolverInterface::setColUpper ;
    /** Set a single column upper bound<br>
     Use DBL_MAX for infinity. */
    virtual void setColUpper( int elementIndex, double elementValue ) {
    }
    
    /** Set a single row lower bound<br>
     Use -DBL_MAX for -infinity. */
    virtual void setRowLower( int elementIndex, double elementValue ) {
        
    }
    
    /** Set a single row upper bound<br>
     Use DBL_MAX for infinity. */
    virtual void setRowUpper( int elementIndex, double elementValue ) {
        
    }
    
    /** Set the type of a single row<br> */
    virtual void setRowType(int index, char sense, double rightHandSide,
                            double range) {
    }
    
    virtual void setRowSetBounds(const int* indexFirst,
                                 const int* indexLast,
                                 const double* boundList);
    
    virtual void setColSetBounds(const int* indexFirst,
                                 const int* indexLast,
                                 const double* boundList);
    //@}
    
    //-------------------------------------------------------------------------
    /**@name Integrality related changing methods */
    //@{
    /** Set the index-th variable to be a continuous variable */
    virtual void setContinuous(int index){}
    /** Set the index-th variable to be an integer variable */
    virtual void setInteger(int index){}
    
    //@}
    
    //-------------------------------------------------------------------------
    /// Set objective function sense (1 for min (default), -1 for max,)
    virtual void setObjSense(double s ) {  }
    
    
    virtual void setColSolution(const double * colsol){}
    virtual void setRowPrice(const double * rowprice){}
    
    //-------------------------------------------------------------------------
    /**@name Methods to expand a problem.<br>
     Note that if a column is added then by default it will correspond to a
     continuous variable. */
    //@{
    
    using OsiSolverInterface::addCol ;
    virtual void addCol(const CoinPackedVectorBase& vec,
                        const double collb, const double colub,
                        const double obj){}
    virtual void deleteCols(const int num, const int * colIndices){}
    using OsiSolverInterface::addRow ;
    virtual void addRow(const CoinPackedVectorBase& vec,
                        const double rowlb, const double rowub){}
    virtual void addRow(const CoinPackedVectorBase& vec,
                        const char rowsen, const double rowrhs,
                        const double rowrng){}
    virtual void deleteRows(const int num, const int * rowIndices);
    
    //@}
    //@}
    
    //---------------------------------------------------------------------------
public:
    
    
    virtual void loadProblem(const CoinPackedMatrix& matrix,
                             const double* collb, const double* colub,
                             const double* obj,
                             const double* rowlb, const double* rowub){}
    virtual void assignProblem(CoinPackedMatrix*& matrix,
                               double*& collb, double*& colub, double*& obj,
                               double*& rowlb, double*& rowub){}
    virtual void loadProblem(const CoinPackedMatrix& matrix,
                             const double* collb, const double* colub,
                             const double* obj,
                             const char* rowsen, const double* rowrhs,
                             const double* rowrng){}
    virtual void assignProblem(CoinPackedMatrix*& matrix,
                               double*& collb, double*& colub, double*& obj,
                               char*& rowsen, double*& rowrhs,
                               double*& rowrng){}
    virtual void loadProblem(const int numcols, const int numrows,
                             const int* start, const int* index,
                             const double* value,
                             const double* collb_, const double* colub_,
                             const double* obj,
                             const double* rowlb_, const double* rowub_);
    virtual void loadProblem(const int numcols, const int numrows,
                             const int* start, const int* index,
                             const double* value,
                             const double* collb, const double* colub,
                             const double* obj,
                             const char* rowsen, const double* rowrhs,
                             const double* rowrng){};
    virtual void writeMps(const char *filename,
                          const char *extension = "mps",
                          double objSense=0.0) const{}
    //@}
    
    //---------------------------------------------------------------------------
    
    /**@name Constructors and destructors */
    //@{
    /// Default Constructor
    OsiVolSolverInterface ();
    
    OsiVolSolverInterface(const char * volparfile);
    
    /// Clone
    virtual OsiSolverInterface * clone(bool copyData = true) const;
    
    /// Copy constructor
    OsiVolSolverInterface (const OsiVolSolverInterface &);
    
    /// Assignment operator
    OsiVolSolverInterface & operator=(const OsiVolSolverInterface& rhs);
    
    /// Destructor
    virtual ~OsiVolSolverInterface ();
    //@}
    
    //---------------------------------------------------------------------------
    
protected:
    ///@name Protected methods
    //@{
    /** Apply a row cut (append to constraint matrix). */
    virtual void applyRowCut(const OsiRowCut& rc){}
    
    /** Apply a column cut (adjust one or more bounds). */
    virtual void applyColCut(const OsiColCut& cc){}
    //@}
    
private:
    //---------------------------------------------------------------------------
    void gutsOfDestructor_();
    //---------------------------------------------------------------------------
    /** A method allocating sufficient space for the rim vectors corresponding
     to the rows. */
    void rowRimAllocator_();
    /** A method allocating sufficient space for the rim vectors corresponding
     to the columns. */
    void colRimAllocator_();
    //---------------------------------------------------------------------------
    /**@name The rim vectors */
    //@{
    /// Pointer to dense vector of structural variable upper bounds
    double  *colub;
    /// Pointer to dense vector of structural variable lower bounds
    double  *collb;
    /// Pointer to dense vector of bool to indicate if column is continuous
    //bool    *continuous_;
    /// Pointer to dense vector of slack variable upper bounds
    double  *rowub;
    /// Pointer to dense vector of slack variable lower bounds
    double  *rowlb;
    /// Pointer to dense vector of row sense indicators
    //char    *rowsense_;
    /// Pointer to dense vector of row right-hand side values
    //double  *rhs_;
    /** Pointer to dense vector of slack upper bounds for range
     constraints (undefined for non-range rows). */
    //double  *rowrange_;
    /// Pointer to dense vector of objective coefficients
    //double  *objcoeffs_;
    //@}
    
    //---------------------------------------------------------------------------
    /**@name The solution */
    //@{
    /// Pointer to dense vector of primal structural variable values
    double  *solution;
    /// Pointer to dense vector of dual row variable values
    double  *dual;
    /// Pointer to dense vector of reduced costs
    double  *rc_;
    /// Pointer to dense vector of left hand sides (row activity levels)
    double  *lhs_;
    /// The Lagrangean cost, a lower bound on the objective value
    double   lagrangeanCost_;
    //@}
    
    //---------------------------------------------------------------------------
    /** An array to store the hotstart information between solveHotStart()
     calls */
    const WarmStartDual *HotStart_;
    
    /// allocated size of the (dualized) row related rim vectors
    int maxNumrows_;
    /// allocated size of the column related rim vectors
    int maxNumcols_;
    /// number of the rows (dualized)
    int numrows_;
    /// number of the columns
    int numcols_;
    
    //---------------------------------------------------------------------------
    //  volume hook methods
    //---------------------------------------------------------------------------
private:
    int compute_rc(const VOL_dvector& dual, VOL_dvector& rc, int actvSSz);
    
    int compute_sg(const VOL_dvector& x, int actvSSz, VOL_dvector& v);
    
    int solve_subproblem(const VOL_dvector& xstar,
                         const VOL_dvector& dual,  VOL_dvector& rc,
                         double& lcost, VOL_dvector& x,
                         double& pcost);
    
    int resolve_subproblem(const VOL_dvector& dual, VOL_dvector& rc,
                           double& lcost,
                           VOL_dvector& x,double& pcost);
    
    int additional_settings(int iter, double& lcost, VOL_dvector& dual, VOL_dvector& rc, VOL_dvector& h, VOL_dvector& x, const VOL_dvector& xhist, const VOL_dvector& dstar,  int actvSSz);
    
    int heuristics(const VOL_problem& p,
                   const VOL_dvector& x, double& heur_val);
    
    int addVI(int iter,double lcost, const VOL_dvector& xstar,
              const VOL_dvector& x, VOL_dvector& dual_lb, VOL_dvector& dual, VOL_dvector& dual_ub,
              VOL_dvector& v, VOL_dvector& h, int & actvSSz);
    
    int removeVI( int & actvSSz,VOL_dvector& pstarv, VOL_dvector& dstaru,  VOL_dvector& dualu);
    
    //---------------------------------------------------------------------------
    //  solving methods
    //---------------------------------------------------------------------------
    
    double knapsack(int a, const double * rc, double*  x);
    
    //---------------------------------------------------------------------------
    //  auxiliary methods
    //---------------------------------------------------------------------------
    
    double arc_dg_imp(int arc, const double * xy, const double * h,  int actvSSz);
    int mark_topo( VOL_dvector& x, double lcost);
    void translate_primal(const VOL_dvector& xhist);
    void translate_hotstart();
    
public:
    void translate_sol();
    void set_start();
    
    void map_duals();
    void direct_solve(const std::deque<int>& topo, const CoinWarmStart* warmstart);
    
    const Data * data;
    int narcs, ndemands, nnodes;
    
    std::deque<int> nz_arcs;
    int szopnd, szunfxd, sznz;
    int fsize, csize;
    
    
    bool HotStartSet;
    int retval;
    
    //Volume attributes
    double VIub, VItt, B0;
    double * VItopo;
    double * addrc;
    double * yhit;
    
    int * arc_map;
    int * actv;
    int lim_to_remv, maxNumVI, intvlVI;
    
    CoverManager* cover_manager;
    CutSetManager* ss_manager;
    int mode;
    
    /// The volume solver
    VOL_problem volprob_;
    VOL_problem* volprob() { return &volprob_; }
    
};




#endif
