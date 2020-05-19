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
	int cut_type;
	MCND_Cut(): BCP_cut_algo(0, 1e40) {cut_type =0; }
	virtual ~MCND_Cut(){};
	
	virtual void pack(BCP_buffer& buf) const =0;
	virtual void unpack(BCP_buffer& buf) =0;
	
	virtual bool check_viol_updt_fix(const BCP_vec<BCP_var*>& vars, BCP_vec<int>& var_changed_pos,
                                BCP_vec<double>& var_new_bd, bool & viol, bool & zrofx, int* fixd)=0;
	virtual bool purgbl()=0;
	virtual void mark_unpurgbl()=0;
	virtual void mark_purgbl()=0;
	virtual int id_vi()=0;
	virtual void print()=0;
	virtual int serial_nmbr()=0;

};

class MCND_CutUnit{
public:
	 
	virtual ~MCND_CutUnit(){};
};

#endif /* MCND_cut_h */
