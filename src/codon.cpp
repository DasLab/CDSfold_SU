#include <sstream>
#include <iostream>
#include <string>
#include "codon.hpp"

/*  Given an amino acid, returns the list of codons that express it.
 *
 *  Inputs:
 *  a - amino acid
 *  exceptedCodons - comma delimited list of codons to ignore
 *
 *  Returns:
 *  Vector of codons that are expressed to the amino acid
 */

auto codon::getCodons(char a, string exceptedCodons) -> vector<string> {
   
    /* check that passed character is an amino acid */
    if (table.count(a) == 0) {
        cerr << "ERR: table doesn't have " << a << "." << endl;
        exit(1);
    }

    /* create vector of codons to skip over */
    char delim = ',';
    vector<string> exceptVector = split(exceptedCodons, delim);
   
    /* map between the codon to skip and the number of times it appears in exceptVector */
    map<string, int> exceptMap;
    for (auto const & tmp : exceptVector) {
        /* exception not found in the map - add it in */ 
        if (exceptMap.find(tmp) == exceptMap.end()) {
            exceptMap[tmp] = 1;
        }
    }

    vector<string> filterTable;
    /* iterate through the codons the express this amino acid */ 
    for (auto const & codon : table.at(a)) {
        
        /* if this codon shouldn't be filtered, then added to the table */
        if (exceptMap.find(codon) == exceptMap.end()) {
            filterTable.push_back(codon);
        }
    }

    return filterTable;
}

/*  Given an amino acid, returns the list of codons that express it using
 *  the extended codon table
 *
 *  Inputs:
 *  a - amino acid
 *  exceptedCodons - comma delimited list of codons to ignore
 *
 *  Returns:
 *  Vector of codons that are expressed to the amino acid
 */

auto codon::getExtendedCodons(char a, string exceptedCodons) -> vector<string> {
    
    /* Check if the amino acid is in the table */
    if (extendedTable.count(a) == 0) {
        cerr << "ERR: extended table doesn't have " << a << "." << endl;
        exit(1);
    }
    /* create a vector of codons to skip */ 
    char delim = ',';
    vector<string> exceptVector = split(exceptedCodons, delim);
    
    map<string, int> exceptMap;
    for (const auto& exceptCodon : exceptVector) {
        string convertedExceptCodon = exceptCodon;   // TODO this is already a string
        map<string, string>::const_iterator itr;
        
        /* if codon appears in the CodonofCodon table, use the table entry */
        if ((itr = codon::expectedCodonOfCodon.find(exceptCodon)) != codon::expectedCodonOfCodon.end()) {
            convertedExceptCodon = itr->second;
        }

        /* convert the list of exceptions to a map */
        if (exceptMap.find(convertedExceptCodon) == exceptMap.end()) {
            exceptMap[convertedExceptCodon] = 1;
        }
    }

    // create a codon table excluding codons in exceptedCodons
    vector<string> filterTable;
    for (auto const & codon : extendedTable.at(a)) {
        if (exceptMap.find(codon) == exceptMap.end()) {
            filterTable.push_back(codon);
        }
    }
    return filterTable;
}


/*  Parses a string of codons and converts it to a vector of strings
 *  Inputs:
 *  str - string holding a list of codons
 *  delim - character separating codons in the string
 *
 *  Return:
 *  vector of strings, each entry representing a codon
 */

auto codon::split(string &str, char delim) -> vector<string> {
    istringstream iss(str);
    string tmp;
    vector<string> res;

    while (getline(iss, tmp, delim)) {
        res.push_back(tmp);
    }
    return res;
}

/* map between codons and codons in the extended table */
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

/* map between amino acids and codons that express it */
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

/* map between amino acids and codons that express it using the extended table*/
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

/* Table mapping codons to amino acids. The dimensions represent the 
 * 1st, 2nd, and 3rd nucleotide in the codon respectively and the index
 * represents which nucleotide is at that position:
 * 1 - A
 * 2 - C
 * 3 - G
 * 4 - U
 */
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
