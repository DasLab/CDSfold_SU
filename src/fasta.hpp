// #include <cstdlib>
//#include <string> // なくても良い。
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
    fasta(const char *fname); //コンストラクタには戻り値
                              //が要らないという特別ルールがある。
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
    vector<eachseq> data; // いくつ配列があるか分からないので、Vector型とする。
    int p;                // ポインタ
    int numSeq;           // 配列数
};
