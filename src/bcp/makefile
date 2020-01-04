##############################################################################
################################ makefile ####################################
##############################################################################
#									     #
#   makefile of BCP 							     #
#									     #
#
##############################################################################

libBcpcg_suffix = \
	BCP_cg_user \
	BCP_cg_param \
	BCP_cg \
	BCP_cg_main 

libBcplp_suffix = \
	BCP_lp_fathom \
	BCP_lp_branch \
	BCP_lp_branching \
	BCP_lp_misc \
	BCP_lp_main_loop \
	BCP_lp_msgproc \
	BCP_lp_user \
	BCP_lp_generate_vars \
	BCP_lp_msg_node_rec \
	BCP_lp_node \
	BCP_lp_create_lp \
	BCP_lp_msg_node_send \
	BCP_lp_convert_OsiWarmStart \
	BCP_lp_result \
	BCP_lp_generate_cuts \
	BCP_lp_colrow \
	BCP_lp \
	BCP_lp_main \
	BCP_lp_pool \
	BCP_lp_param 
	
libBcptm_suffix = \
	BCP_tm_trimming \
	BCP_tm_msg_node_send \
	BCP_tm_msg_node_rec \
	BCP_tm_param \
	BCP_tm_user \
	BCP_tm_commandline \
	BCP_tm_node \
	BCP_tm_nodes_to_storage \
	BCP_tm_msgproc \
	BCP_tm_functions \
	BCP_tm_main \
	BCP_tm_statistics \
	BCP_tm \
	BCP_tmstorage 

libBcpvg_suffix = \
	BCP_vg_param \
	BCP_vg_user \
	BCP_vg_main \
	BCP_vg 
	
libBcpm_suffix = \
	BCP_solution \
	BCP_node_change \
	BCP_indexed_pricing \
	BCP_USER \
	BCP_warmstart_basis \
	BCP_branch \
	BCP_vector_sanity \
	BCP_cut \
	BCP_obj_change \
	BCP_warmstart_dual \
	BCP_warmstart_primaldual \
	BCP_matrix \
	BCP_problem_core \
	BCP_process \
	BCP_warmstart_pack \
	BCP_matrix_pack \
	BCP_var \
	BCP_vector_bool \
	BCP_vector_char \
	BCP_vector_double \
	BCP_vector_int \
	BCP_vector_short \
	BCP_message_single \
	BCP_message_mpi \
	BCP_message_pvm
	
	
includebcp_HEADERS = \
	include/BcpConfig.h \
	include/BCP_branch.hpp \
	include/BCP_buffer.hpp \
	include/BCP_cg.hpp \
	include/BCP_cg_param.hpp \
	include/BCP_cg_user.hpp \
	include/BCP_cut.hpp \
	include/BCP_enum_branch.hpp \
	include/BCP_enum.hpp \
	include/BCP_enum_process_t.hpp \
	include/BCP_enum_tm.hpp \
	include/BCP_error.hpp \
	include/BCP_functions.hpp \
	include/BCP_indexed_pricing.hpp \
	include/BCP_lp_branch.hpp \
	include/BCP_lp_functions.hpp \
	include/BCP_lp.hpp \
	include/BCP_lp_main_loop.hpp \
	include/BCP_lp_node.hpp \
	include/BCP_lp_param.hpp \
	include/BCP_lp_pool.hpp \
	include/BCP_lp_result.hpp \
	include/BCP_lp_user.hpp \
	include/BCP_main_fun.hpp \
	include/BCP_math.hpp \
	include/BCP_matrix.hpp \
	include/BCP_mempool.hpp \
	include/BCP_message.hpp \
	include/BCP_message_mpi.hpp \
	include/BCP_message_pvm.hpp \
	include/BCP_message_single.hpp \
	include/BCP_message_tag.hpp \
	include/BCP_node_change.hpp \
	include/BCP_obj_change.hpp \
	include/BCP_os.hpp \
	include/BCP_parameters.hpp \
	include/BCP_problem_core.hpp \
	include/BCP_process.hpp \
	include/BCP_set_intersects.hpp \
	include/BCP_solution.hpp \
	include/BCP_string.hpp \
	include/BCP_tm_functions.hpp \
	include/BCP_tm.hpp \
	include/BCP_tm_node.hpp \
	include/BCP_tm_param.hpp \
	include/BCP_tmstorage.hpp \
	include/BCP_tm_user.hpp \
	include/BCP_USER.hpp \
	include/BCP_var.hpp \
	include/BCP_vector_bool.hpp \
	include/BCP_vector_bool.hpp \
	include/BCP_vector_change.hpp \
	include/BCP_vector_char.hpp \
	include/BCP_vector_char.hpp \
	include/BCP_vector_double.hpp \
	include/BCP_vector_double.hpp \
	include/BCP_vector_general.hpp \
	include/BCP_vector.hpp \
	include/BCP_vector_int.hpp \
	include/BCP_vector_int.hpp \
	include/BCP_vector_sanity.hpp \
	include/BCP_vector_short.hpp \
	include/BCP_vector_short.hpp \
	include/BCP_vg.hpp \
	include/BCP_vg_param.hpp \
	include/BCP_vg_user.hpp \
	include/BCP_warmstart_basis.hpp \
	include/BCP_warmstart_dual.hpp \
	include/BCP_warmstart_primaldual.hpp \
	include/BCP_warmstart.hpp
	

libBcp_suffix =  $(libBcpcg_suffix) $(libBcptm_suffix) $(libBcpvg_suffix) $(libBcpm_suffix) $(libBcplp_suffix)
LINK = -Wl,-rpath,../../lib -L../../lib -lCoinUtils -lOsi 

COININC = -I../CoinUtils/include/ -I../CoinUtils/config/
OSIINC = -I../Osi/include/

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

#STATIC LIBRARY
#COINLIB = -L$(COINDIR)build/ -lCoinUtils 
#libCoinUtils_ = $(COINDIR)build/libCoinUtils.a 

#SHARED LIBRARY
libBcp_ = ../../lib/libBcp.so

BCPOBJ = $(libBcp_suffix:%=build/%.o)
BCPINC = -Iinclude/ -Iconfig/
BCPH   = $(includebcp_HEADERS)

# dependencies: every .o from its .C + every recursively included .h- - - - -

all: $(libBcp_suffix)  $(libBcp_) 

$(libBcpcg_suffix):
	g++ -fPIC -c cg/$@.cpp -o build/$@.o $(COININC) $(OSIINC) $(BCPINC) -O3

$(libBcpvg_suffix):
	g++ -fPIC -c vg/$@.cpp -o build/$@.o $(COININC) $(OSIINC) $(BCPINC) -O3

$(libBcptm_suffix):
	g++ -fPIC -c tm/$@.cpp -o build/$@.o $(COININC) $(OSIINC) $(BCPINC) -O3
	
$(libBcpm_suffix):
	g++ -fPIC -c member/$@.cpp -o build/$@.o $(COININC) $(OSIINC) $(BCPINC) -O3
	
$(libBcplp_suffix):
	g++ -fPIC -c lp/$@.cpp -o build/$@.o $(COININC) $(OSIINC) $(BCPINC) -O3
	
#SHARED LIBRARY
$(libBcp_):
	g++ -shared  -Wl,-install_name,/Users/ruishibasaki/Documents/thesis/development/mcnd/bcp/lib/libBcp.so -o $@ $(LINK) $(BCPOBJ)

# STATIC LIBRARY
#$(libCoinUtils_): $(COINOBJ)
#	g++ -shared

clean::
	rm -f build/*.o build/*~ ../../lib/libBcp.so


############################ End of makefile #################################
