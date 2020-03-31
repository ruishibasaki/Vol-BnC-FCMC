//
//  MCND_cut.hpp
//  
//
//  Created by Rui Shibasaki on 09/01/2020.
//

#ifndef MCND_cut_h
#define MCND_cut_h

#include "BCP_cut.hpp"

class CutManager{
protected:
	static int ttgend;
	CutManager(){ }
	virtual ~CutManager(){};
};


class MCND_Cut : public BCP_cut_algo{
public:
	int type;
	MCND_Cut(): BCP_cut_algo(0, 1e40) {type =0; }
	virtual ~MCND_Cut(){};
	
	virtual void pack(BCP_buffer& buf) const =0;
	virtual void unpack(BCP_buffer& buf) =0;
	
	virtual bool check_viol(const BCP_vec<BCP_var*>& vars)=0;
	virtual double check_viol(const double* vars)=0;
	virtual bool check_logical_fix(const BCP_vec<BCP_var*>& vars, int* yarcs)=0;
	virtual bool purgbl()=0;
	virtual void mark_unpurgbl()=0;
	virtual int id_vi()=0;
	virtual int serial_nmbr()=0;

};

class MCND_CutUnit{
public:
	 
	virtual ~MCND_CutUnit(){};
};

#endif /* MCND_cut_h */
