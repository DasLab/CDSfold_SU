#include <iostream>
#include <cstdlib>
//#include <string> // $B$J$/$F$bNI$$!#(B
#include <fstream>
#include <cstring> // to use strcpy for gcc-4.3 or later
#include <vector>
using namespace std;

#pragma once

class eachseq {
  public:
    char *desc;
    char *seq;
    int seqlen;
};

class fasta {
  public:
    fasta(const char *fname); //$B%3%s%9%H%i%/%?$K$OLa$jCM(B
                              //$B$,MW$i$J$$$H$$$&FCJL%k!<%k$,$"$k!#(B
    ~fasta();
    auto getDesc() -> char * { return data[p].desc; }
    auto getSeq() -> char * { return data[p].seq; }
    auto getSeqLen() -> int { return data[p].seqlen; }
    void initP() { p = 0; }
    void printP() { cout << p << endl; }

    auto next() -> int {
        p++;
        if (p > numSeq) {
            p = 0;
            return 0;
        }
        return 1;
    }

  private:
    vector<eachseq> data; // $B$$$/$DG[Ns$,$"$k$+J,$+$i$J$$$N$G!"(BVector$B7?$H$9$k!#(B
    int p;                // $B%]%$%s%?(B
    int numSeq;           // $BG[Ns?t(B
};
