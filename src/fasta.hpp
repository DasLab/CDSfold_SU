#include <iostream>
#include <cstdlib>
//#include <string> // なくても良い。
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
    fasta(const char *fname); //コンストラクタには戻り値
                              //が要らないという特別ルールがある。
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
    vector<eachseq> data; // いくつ配列があるか分からないので、Vector型とする。
    int p;                // ポインタ
    int numSeq;           // 配列数
};
