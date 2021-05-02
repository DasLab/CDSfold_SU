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
    string desc_;
    string seq_;
};

/* Class for holding a collection of amino acid sequences in the
 * FASTA format
 */
class fasta {
  public:
    
    /* constructor - pass in file name */ 
    fasta(const char *fname);
   
    auto getDesc() -> string { return data_[p_].desc_; }
    auto getSeq() -> string  { return data_[p_].seq_; }
    auto getSeqLen() -> int  { return data_[p_].seq_.size(); }
    void initP() { p_ = 0; }
    void printP();

    /* iterator used to move onto next instance of eachseq in the data vector
     * input: none
     * returns: (int) if there is another sequence that in the vector */
    auto next() -> int {
        p_++;
        if (p_ > numSeq_) {
            p_ = 0;              // TODO - this should return bool not int
            return false;
        }
        return true;
    }

  private:
    vector<eachseq> data_; // vector holding each individual sequence
    int p_;                // index into the list of sequences 
    int numSeq_;           // number of sequences
};
