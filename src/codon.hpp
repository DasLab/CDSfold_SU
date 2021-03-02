#include <map>
// #include <string>
#include <iostream>
// #include <iosfwd>
#include <vector>

#pragma once

using namespace std;

/* tableとextendedTableはL、Rが異なる
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

 V、X：次の塩基はAまたはG
 W、Y：次の塩基はCまたはU
 */

class codon {
  public:
    codon();

    auto getCodons(char c, string exceptedCodons) -> vector<string>;
    auto getExtendedCodons(char c, string exceptedCodons) -> vector<string>;

    auto c2a(int p1, int p2, int p3) const -> char { return table_rev[p1][p2][p3]; }

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
    static const map<char, vector<string>> table;
    static const map<char, vector<string>> extendedTable;
    static const map<string, string> expectedCodonOfCodon;
    char table_rev[5][5][5];

    auto split(string &str, char delim) -> vector<string>;
};
