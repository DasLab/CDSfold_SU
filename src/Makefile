######## Please set this valuable BY YOURSELF ##############
VIENNA=${VIENNA_INSTALL}
CDSHOME=${CDS_HOME}
############################################################

all: compile link

CXX = g++

CXXFLAGS = -O3 -std=c++17 -mfpmath=sse -msse -msse2 -msse3 \
	-DNDEBUG -Wall -Wextra -Wpedantic -pipe -Wundef --param large-function-growth=100000 -fPIC \
	-I$(VIENNA) -I$(VIENNA)/include/ 

CDSFOLD_OBJS = fasta.o Util.o FileManager.o codon.o Problem.o AASeqConverter.o ViennaEnergyModel.o

compile:
	$(CXX) ${CXXFLAGS} -c CDSfold.cpp 

link: $(CDSFOLD_OBJS)
	g++ -Wall -Wextra -Wpedantic -o CDSfold CDSfold.o $(CDSFOLD_OBJS) -fmessage-length=0 -L$(VIENNA) -fopenmp -lRNA

%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c $<

example:
	$(CDSHOME)/src/CDSfold $(CDSHOME)/example/test.faa

test:
	python3 $(CDSHOME)/test/run_tests.py

clean:
	-rm *.o
	-rm CDSfold
