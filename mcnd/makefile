# compiler
CC = g++

# module name
NAME = mcnd_bcp

# basic directory
DIR = ./

# debug switches
#SW = -g -Wall -fPIC -fexceptions -DNDEBUG -DIL_STD
# production switches
SW = -O3

suffix_cg = \
	covercollection \
	covermanager \
	cutsetmanager \
	cutsetcollection \
	lift 

suffix_tm = \
	MCND_tm \
	#MCND_tm_param 
	
suffix_lp = \
	MCND_lp \
	MCND_lp_branch  

suffix_lp_cplex = \
	MCND_mincostflow \
	#MCND_pump \
	#MCND_checkfeas
	
suffix_member =\
	UtilsMethods	\
	MCND_data	\
	MCND_init	\
	MCND_solution \
	MCND_branch_score

suffix_solver = OsiVolSolverInterface VolVolume
	
OBJ = $(suffix_tm:%=build/%.o) $(suffix_lp_cplex:%=build/%.o) $(suffix_lp:%=build/%.o) $(suffix_member:%=build/%.o) $(suffix_solver:%=build/%.o) $(suffix_cg:%=build/%.o)

#SYSTEM = x86-64_linux
#LIBFORMAT = static_pic
#CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio128/concert
#CPLEXDIR = 	/opt/ibm/ILOG/CPLEX_Studio128/cplex
#CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

SYSTEM = x86-64_osx
LIBFORMAT = static_pic
CONCERTDIR = /Applications/CPLEX_Studio128/concert
CPLEXDIR = 	/Applications/CPLEX_Studio128/cplex
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXFLAGS =   -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -ldl -lpthread

	
CPLEXINC = -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include
COININC = -I../src/CoinUtils/include/ -I../src/CoinUtils/config/ 
OSIINC = -I../src/Osi/include/
VOLINC = -I../src/Vol/include/
BCPINC = -I../src/Bcp/include/ -I../src/Bcp/config/
INC= -I$(DIR)include

LINKPATH = /Users/ruishibasaki/Documents/thesis/development/mcnd/bcp/lib/
LINK = -Wl,-rpath,$(LINKPATH) -L$(LINKPATH) -lCoinUtils -lOsi -lBcp
CPLEXFLAGS =   -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -ldl -lpthread

tm_o= $(suffix_tm:%=build/%.o)
cg_o= $(suffix_cg:%=build/%.o)
member_o= $(suffix_member:%=build/%.o)
lp_o= $(suffix_lp:%=build/%.o)
lp_cplex_o= $(suffix_lp_cplex:%=build/%.o)
solver_o= $(suffix_solver:%=build/%.o)


all: $(NAME)

$(cg_o): build/%.o: cg/%.cpp
	g++ -c $< -o $@ $(COININC) $(OSIINC) $(BCPINC) $(INC) -O3
	
$(tm_o): build/%.o: tm/%.cpp
	g++ -c $< -o $@ $(COININC) $(OSIINC) $(BCPINC) $(INC) -O3
	
$(member_o): build/%.o: member/%.cpp
	g++ -c $< -o $@ $(COININC) $(OSIINC) $(BCPINC) $(INC) -O3
	
$(solver_o): build/%.o: solver/%.cpp
	g++ -c $< -o $@ $(COININC) $(OSIINC) $(VOLINC) $(BCPINC) $(INC) -O3

$(lp_cplex_o): build/%.o: lp/%.cpp
	g++ -c $< -o $@ $(CPLEXINC) $(COININC) $(OSIINC) $(BCPINC) $(INC) -O3
	
$(lp_o): build/%.o: lp/%.cpp
	g++ -c $< -o $@  $(CPLEXINC) $(COININC) $(OSIINC) $(VOLINC) $(BCPINC) $(INC) -O3
	

$(NAME): $(lp_cplex_o) $(cg_o) $(tm_o) $(member_o) $(solver_o)  $(lp_o) 
	$(CC) $(OBJ) -o $(NAME)  $(LINK) $(CPLEXFLAGS) $(INC) $(SW)

clean::
	rm -f $(DIR)build/*.o $(DIR)build/*~ $(NAME)
