//
//  MCND_cut.hpp
//  
//
//  Created by Rui Shibasaki on 09/01/2020.
//

#ifndef MCND_cut_h
#define MCND_cut_h

#include "BCP_cut.hpp"
#include "BCP_buffer.hpp"

#include "covercollection.hpp"


class CoverCut : public BCP_cut_algo{

	Cover* cover;

public:
	
	inline CoverCut(Cover* c): BCP_cut_algo(0, 1e40), cover(c){}
	inline CoverCut(): BCP_cut_algo(0, 1e40),cover(0){}
	inline ~CoverCut(){}

	
	inline void pack(BCP_buffer& buf) const{ buf.pack(cover); }
	inline void unpack(BCP_buffer& buf){ buf.unpack(cover);}
	
	inline Cover* get_cover(){ return cover;}
	inline void set_cover(Cover* c){ cover = c;}
	
	bool check_viol(const BCP_vec<BCP_var*>& vars);
	double check_viol(const double* vars);

};

#endif /* MCND_cut_h */
