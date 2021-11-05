# Copyright 2020, Gurobi Optimization, LLC

PLATFORM = linux64
INC      = /opt/gurobi912/linux64/include
CPP      = g++ -lemon
#CPP      = g++ -02
CARGS    = -m64 -g -std=c++11
CPPLIB   = /opt/gurobi912/linux64/lib 

CCFLAGS = $(CARGS)

all: VAP_Main 


VAP_Main: VAP_Main.cpp
	$(CPP) $(CARGS) -o $@ $< -I$(INC) -L$(CPPLIB) -lgurobi_c++ -lgurobi91 -lm
	

clean:
	rm -rf VAP_Main *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm; \

	
