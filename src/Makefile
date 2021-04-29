######## Use environment variables - set by user ###########
VIENNA=${VIENNA_INSTALL}
CDSHOME=${CDS_HOME}
############################################################

all: vienna

dummy: compile link

########## Default compiler, flags, and objects ###########
CXX = g++

CXXFLAGS = -O3 -std=c++17 -mfpmath=sse -msse -msse2 -msse3 \
	-Wall -Wextra -Wpedantic -pipe -Wundef --param large-function-growth=100000 -fPIC\

CDSFOLD_OBJS = fasta.o Util.o FileManager.o codon.o Problem.o AASeqConverter.o

LD_FLAGS = -fmessage-length=0 -fopenmp 

################### Vienna build  ##########################

VIENNA_FLAGS = -I$(VIENNA) -I$(VIENNA)/include/ -DUSE_VIENNA_ENERGY_MODEL 
VIENNA_OBJS = ViennaEnergyModel.o
VIENNA_LD_FLAGS = -L$(VIENNA) -lRNA

vienna: CXXFLAGS+=$(VIENNA_FLAGS)
vienna: LD_FLAGS+=$(VIENNA_LD_FLAGS)
vienna: CDSFOLD_OBJS+=$(VIENNA_OBJS)
vienna: compile $(CDSFOLD_OBJS) $(VIENNA_OBJS) link

################# Debug build ##############################
DEBUG_FLAGS = -DDEBUG -g -fverbose-asm

debug: CXXFLAGS+= $(DEBUG_FLAGS)
debug: vienna

################### General Build rules ####################

compile:
	$(CXX) ${CXXFLAGS} -c CDSfold.cpp 

link: $(CDSFOLD_OBJS)
	g++ -Wall -Wextra -Wpedantic -o CDSfold CDSfold.o $(CDSFOLD_OBJS) $(LD_FLAGS)

%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c $<

example:
	$(CDSHOME)/src/CDSfold $(CDSHOME)/example/test.faa

test:
	python3 $(CDSHOME)/test/run_tests.py

clean:
	-rm *.o
	-rm *.s
	-rm CDSfold

.PHONY: all clean dummy test link compile vienna
