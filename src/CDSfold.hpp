#include <map>
#include <cstring>
#include <fstream>
#include <sstream>
#include "codon.hpp"
#include "backtracking.hpp" // just for bond
#include "energy.hpp"

using namespace std;

#pragma once



inline void showPos2Nuc(vector<vector<int>> &v) {

    for (unsigned int i = 1; i <= v.size(); i++) {
        cout << i;
        for (unsigned int j = 0; j < v[i].size(); j++) {
            cout << "\t" << v[i][j];
        }
        cout << endl;
    }
}


inline void showPos2Nuc(vector<vector<int>> &v, std::array<char, 20> const & i2n) {

    for (unsigned int i = 1; i < v.size(); i++) {
        //	  for(unsigned int i = 1; i <= v.size(); i++){ // This is segmentaion fault
        cout << i;
        for (unsigned int j = 0; j < v[i].size(); j++) {
            cout << "\t" << i2n[v[i][j]];
        }
        cout << endl;
    }
}


