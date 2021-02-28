######## Please set this valuable BY YOURSELF ##############
VIENNA=/home/andy/programs/ViennaRNA-2.4.9/src/ViennaRNA
############################################################

all: compile link

CXX = g++

CXXFLAGS = -O3 -std=c++17 -mfpmath=sse -msse -msse2 -msse3 \
	-DNDEBUG -Wall -Wextra -Wpedantic -pipe -Wundef --param large-function-growth=100000 -fPIC \
	-I$(VIENNA) -I$(VIENNA)/include/ViennaRNA 

CDSFOLD_OBJS = fasta.o Util.o FileManager.o codon.o backtracking.o

compile:
	$(CXX) ${CXXFLAGS} -c CDSfold.cpp 

link: $(CDSFOLD_OBJS)
	g++ -Wall -Wextra -Wpedantic -o CDSfold CDSfold.o $(CDSFOLD_OBJS) -fmessage-length=0 -L$(VIENNA) -fopenmp -lRNA

%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c $<

clean:
	-rm *.o
	-rm CDSfold
