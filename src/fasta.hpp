// #include <cstdlib>
//#include <string> // なくても良い。
// #include <fstream>
#include <vector>
#include <string>

using namespace std;

#pragma once

/* Class for holding a single sequence
 *
 */
class eachseq {
  public:
    string desc;
    string seq;
};

/* Class for holding a collection of amino acid sequences in the
 * FASTA format
 */
class fasta {
  public:
    
    /* constructor - pass in file name */ 
    fasta(const char *fname);
   
    /* destructor - TODO this isn't needed since it does nothing */
    ~fasta();
    auto getDesc() -> string { return data[p].desc; }
    auto getSeq() -> string { return data[p].seq; }
    auto getSeqLen() -> int { return data[p].seq.size(); }
    void initP() { p = 0; }
    void printP();

    /* iterator used to move onto next instance of eachseq in the data vector
     * input: none
     * returns: (int) if there is another sequence that in the vector */
    auto next() -> int {
        p++;
        if (p > numSeq) {
            p = 0;              // TODO - this should return bool not int
            return 0;
        }
        return 1;
    }

  private:
    vector<eachseq> data; // vector holding each individual sequence
    int p;                // index into the list of sequences 
    int numSeq;           // number of sequences
};
