#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
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

    vector<string> getCodons(char c, string exceptedCodons) {
        if (table.count(c) == 0) {
            cerr << "ERR: table doesn't have " << c << "." << endl;
            exit(1);
        }

        char delim = ',';
        vector<string> exceptVector = split(exceptedCodons, delim);
        map<string, int> exceptMap;
        for (auto const & tmp : exceptVector) {
            if (exceptMap.find(tmp) == exceptMap.end()) {
                exceptMap[tmp] = 1;
            } else {
                int count = exceptMap[tmp] + 1;
                exceptMap[tmp] = count;
            }
        }

        vector<string> filterTable;
        for (auto const & codon : table.at(c)) {
            if (exceptMap.find(codon) == exceptMap.end()) {
                filterTable.push_back(codon);
            }
        }

        return filterTable;
    }

    vector<string> getExtendedCodons(char c, string exceptedCodons) {
        if (extendedTable.count(c) == 0) {
            cerr << "ERR: extended table doesn't have " << c << "." << endl;
            exit(1);
        }
        char delim = ',';
        vector<string> exceptVector = split(exceptedCodons, delim);
        map<string, int> exceptMap;
        for (const auto& exceptCodon : exceptVector) {
            string convertedExceptCodon = exceptCodon;
            map<string, string>::const_iterator itr;
            if ((itr = codon::expectedCodonOfCodon.find(exceptCodon)) != codon::expectedCodonOfCodon.end()) {
                convertedExceptCodon = itr->second;
            }

            if (exceptMap.find(convertedExceptCodon) == exceptMap.end()) {
                exceptMap[convertedExceptCodon] = 1;
            } else {
                int count = exceptMap[convertedExceptCodon] + 1;
                exceptMap[convertedExceptCodon] = count;
            }
        }

        // 除外コドンを除いたコドンテーブルを作成
        vector<string> filterTable;
        for (auto const & codon : extendedTable.at(c)) {
            if (exceptMap.find(codon) == exceptMap.end()) {
                filterTable.push_back(codon);
            }
        }
        return filterTable;
    }

    char c2a(int p1, int p2, int p3) const { return table_rev[p1][p2][p3]; }

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

    vector<string> split(string &str, char delim) {
        istringstream iss(str);
        string tmp;
        vector<string> res;

        while (getline(iss, tmp, delim)) {
            res.push_back(tmp);
        }
        return res;
    }
};

const map<string, string> codon::expectedCodonOfCodon{
	{"UUA", "UVA"},
	{"UUG", "UVG"},
	{"CUU", "CWU"},
	{"CUC", "CWC"},
	{"CUA", "CVA"},
	{"CUG", "CVG"},
	{"CGU", "CYU"},
	{"CGC", "CYC"},
	{"CGA", "CXA"},
	{"CGG", "CXG"},
	{"AGA", "AXA"},
	{"AGG", "AXG"},
};

const map<char, vector<string>> codon::table{
	{'F', {
		"UUU", "UUC"
	}},
	{'L', {
		"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"
	}},
	{'I', {
		"AUU", "AUC", "AUA"
	}},
	{'M', {
		"AUG"
	}},
	{'V', {
		"GUU", "GUC", "GUA", "GUG"
	}},
	{'S', {
		"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"
	}},
	{'P', {
		"CCU", "CCC", "CCA", "CCG"
	}},
	{'T', {
		"ACU", "ACC", "ACA", "ACG"
	}},
	{'A', {
		"GCU", "GCC", "GCA", "GCG"
	}},
	{'Y', {
		"UAU", "UAC"
	}},
	{'*', {
		"UAA", "UAG", "UGA"
	}},
	{'H', {
		"CAU", "CAC"
	}},
	{'Q', {
		"CAA", "CAG"
	}},
	{'N', {
		"AAU", "AAC"
	}},
	{'K', {
		"AAA", "AAG"
	}},
	{'D', {
		"GAU", "GAC"
	}},
	{'E', {
		"GAA", "GAG"
	}},
	{'C', {
		"UGU", "UGC"
	}},
	{'W', {
		"UGG"
	}},
	{'R', {
		"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"
	}},
	{'G', {
		"GGU", "GGC", "GGA", "GGG"
	}}
};

const map<char, vector<string>> codon::extendedTable{
	{'F', {
		"UUU", "UUC"
	}},
	{'L', {
		"UVA", "UVG", "CWU", "CWC", "CVA", "CVG"
	}},
	{'I', {
		"AUU", "AUC", "AUA"
	}},
	{'M', {
		"AUG"
	}},
	{'V', {
		"GUU", "GUC", "GUA", "GUG"
	}},
	{'S', {
		"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"
	}},
	{'P', {
		"CCU", "CCC", "CCA", "CCG"
	}},
	{'T', {
		"ACU", "ACC", "ACA", "ACG"
	}},
	{'A', {
		"GCU", "GCC", "GCA", "GCG"
	}},
	{'Y', {
		"UAU", "UAC"
	}},
	{'*', {
		"UAA", "UAG", "UGA"
	}},
	{'H', {
		"CAU", "CAC"
	}},
	{'Q', {
		"CAA", "CAG"
	}},
	{'N', {
		"AAU", "AAC"
	}},
	{'K', {
		"AAA", "AAG"
	}},
	{'D', {
		"GAU", "GAC"
	}},
	{'E', {
		"GAA", "GAG"
	}},
	{'C', {
		"UGU", "UGC"
	}},
	{'W', {
		"UGG"
	}},
	{'R', {
		"CYU", "CYC", "CXA", "CXG", "AXA", "AXG"
	}},
	{'G', {
		"GGU", "GGC", "GGA", "GGG"
	}}
};

codon::codon() {
    table_rev[4][4][4] = 'F';
    table_rev[4][4][2] = 'F';
    table_rev[4][4][1] = 'L';
    table_rev[4][4][3] = 'L';

    table_rev[2][4][4] = 'L';
    table_rev[2][4][2] = 'L';
    table_rev[2][4][1] = 'L';
    table_rev[2][4][3] = 'L';

    table_rev[1][4][4] = 'I';
    table_rev[1][4][2] = 'I';
    table_rev[1][4][1] = 'I';
    table_rev[1][4][3] = 'M';

    table_rev[3][4][4] = 'V';
    table_rev[3][4][2] = 'V';
    table_rev[3][4][1] = 'V';
    table_rev[3][4][3] = 'V';

    table_rev[4][2][4] = 'S';
    table_rev[4][2][2] = 'S';
    table_rev[4][2][1] = 'S';
    table_rev[4][2][3] = 'S';

    table_rev[2][2][4] = 'P';
    table_rev[2][2][2] = 'P';
    table_rev[2][2][1] = 'P';
    table_rev[2][2][3] = 'P';

    table_rev[1][2][4] = 'T';
    table_rev[1][2][2] = 'T';
    table_rev[1][2][1] = 'T';
    table_rev[1][2][3] = 'T';

    table_rev[3][2][4] = 'A';
    table_rev[3][2][2] = 'A';
    table_rev[3][2][1] = 'A';
    table_rev[3][2][3] = 'A';

    table_rev[4][1][4] = 'Y';
    table_rev[4][1][2] = 'Y';
    table_rev[4][1][1] = '*';
    table_rev[4][1][3] = '*';

    table_rev[2][1][4] = 'H';
    table_rev[2][1][2] = 'H';
    table_rev[2][1][1] = 'Q';
    table_rev[2][1][3] = 'Q';

    table_rev[1][1][4] = 'N';
    table_rev[1][1][2] = 'N';
    table_rev[1][1][1] = 'K';
    table_rev[1][1][3] = 'K';

    table_rev[3][1][4] = 'D';
    table_rev[3][1][2] = 'D';
    table_rev[3][1][1] = 'E';
    table_rev[3][1][3] = 'E';

    table_rev[4][3][4] = 'C';
    table_rev[4][3][2] = 'C';
    table_rev[4][3][1] = '*';
    table_rev[4][3][3] = 'W';

    table_rev[2][3][4] = 'R';
    table_rev[2][3][2] = 'R';
    table_rev[2][3][1] = 'R';
    table_rev[2][3][3] = 'R';

    table_rev[1][3][4] = 'S';
    table_rev[1][3][2] = 'S';
    table_rev[1][3][1] = 'R';
    table_rev[1][3][3] = 'R';

    table_rev[3][3][4] = 'G';
    table_rev[3][3][2] = 'G';
    table_rev[3][3][1] = 'G';
    table_rev[3][3][3] = 'G';
}
