// #include <cstdlib>
//#include <string> // $B$J$/$F$bNI$$!#(B
// #include <fstream>
#include <vector>
#include <string>
using namespace std;

#pragma once

class eachseq {
  public:
    string desc;
    string seq;
};

class fasta {
  public:
    fasta(const char *fname); //$B%3%s%9%H%i%/%?$K$OLa$jCM(B
                              //$B$,MW$i$J$$$H$$$&FCJL%k!<%k$,$"$k!#(B
    ~fasta();
    auto getDesc() -> string { return data[p].desc; }
    auto getSeq() -> string { return data[p].seq; }
    auto getSeqLen() -> int { return data[p].seq.size(); }
    void initP() { p = 0; }
    void printP();

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
