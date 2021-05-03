#include <map>
#include <iostream>
#include <vector>

#pragma once

using namespace std;

/* The table and extended table have differing L and R
 table	extendedTable
 L UUA <-> UVA
   UUG <-> UVG
   CUA <-> CVA
   CUC <-> CWC
   CUG <-> CVG
   CUU <-> CWU
 R AGA <-> AXA
   AGG <-> AXG
   CGU <-> CYU
   CGC <-> CYC
   CGA <-> CXA
   CGG <-> CXG

 V, X：The next base is A or G
 W, Y：The next base is C or U
 */

/* class containing information */
class codon {
  public:
    
    /* constructor */
    codon();

    /* Get codons that express a certain amino acid */
    auto getCodons(char a, string exceptedCodons) -> vector<string>;
    
    /* Get codons that express a certain amino acid using the extended table */
    auto getExtendedCodons(char a, string exceptedCodons) -> vector<string>;

    auto c2a(int p1, int p2, int p3) const -> char { return codonToAA[p1][p2][p3]; }

    /* print out the table mapping amino acids to codons */
    void showTable() {
        for (auto & it : table) {
            char aa = it.first;
            vector<string> codons = getCodons(aa, "");
            for (auto & codon : codons) {
                cout << aa << " " << codon << endl;
            }
        }
    }

  private:
    static const map<char, vector<string>> table;           /* amino acids -> all possible codons */
    static const map<char, vector<string>> extendedTable;   /* amino acids -> codons w/ extended table */
    static const map<string, string> expectedCodonOfCodon;  /* codon to extended table codon */
    vector<vector<vector<char>>> codonToAA;                 /* table mapping codons to amino acids  */

    auto split(string &str, char delim) -> vector<string>;  /* convert string into  a vector */
};
