/*
 * Util.h
 *
 *  Created on: 2014/09/22
 *      Author: kamegai
 */

#pragma once

#define MIN2(A, B) ((A) < (B) ? (A) : (B))
#define MAX2(A, B) ((A) > (B) ? (A) : (B))

#include <iosfwd>
#include <vector>
#include <random>
#include <memory>

extern "C" {
    #include <climits>
}

using namespace std;

class Util {
    public:
        static void baseReplace(string &base, const string& from, const string& to);
};

inline void shuffleStr(vector<string>(*ary), int size, std::unique_ptr< std::mt19937 > const & gen) {
    std::uniform_int_distribution<> distrib(0, size - 1);
    for (int i = 0; i < size; i++) {
        int j = distrib(*gen);
        string t = (*ary)[i];
        (*ary)[i] = (*ary)[j];
        (*ary)[j] = t;
    }
}

inline void shuffle(int ary[], int size, std::unique_ptr< std::mt19937 > const & gen) {
    std::uniform_int_distribution<> distrib(0, size - 1);
    for (int i = 0; i < size; i++) {
        int j = distrib(*gen);
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }
}

inline auto predict_memory(int len, int w, vector<vector<int>> &pos2nuc) -> float {
    // int size = getMatrixSize(len, w);
    long int total_bytes = 0;
    // int n_test;
    //	int limit = LONG_MAX;
    for (int i = 1; i <= len; i++) {
        for (int j = i; j <= MIN2(len, i + w - 1); j++) {
            // cout << i << "," << j << endl;
            total_bytes += sizeof(int *) * (pos2nuc[i].size() + 4) * 2; // +4 is an empirical value
            for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                total_bytes += sizeof(int) * (pos2nuc[j].size() + 4) * 2; // +4 is an empirical value
                // n_test++;
                if (total_bytes >= LONG_MAX - 10000)
                    return (float)LONG_MAX / (1024 * 1024);
            }
        }
    }
    // cout << n_test << endl;
    // cout << total_bytes << endl;
    total_bytes *= 1.2; // 1.2 is an empirical value
    if (total_bytes >= LONG_MAX - 10000)
        return (float)LONG_MAX / (1024 * 1024);

    return (float)total_bytes / (1024 * 1024);
}


