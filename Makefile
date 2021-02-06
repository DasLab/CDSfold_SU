######## Please set this valuable BY YOURSELF ##############
VIENNA=/home/andy/programs/ViennaRNA-2.4.9/src/ViennaRNA
############################################################

all: compile link

compile:
	cd src; g++ -c -O3 -Wall -Wextra -Wpedantic CDSfold.cpp -I$(VIENNA) -I$(VIENNA)/include/ViennaRNA 

link:
	cd src; g++ -Wall -Wextra -Wpedantic -o CDSfold CDSfold.o -fmessage-length=0 -L$(VIENNA) -fopenmp -lRNA

clean:
	-rm src/CDSfold.o
	-rm src/CDSfold
