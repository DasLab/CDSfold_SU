#pragma once

#include <iosfwd>
#include <iostream>
#include <memory>
// maybe convert to unique_ptr to avoid need for whole headers
#include "codon.hpp"
#include "CDSfold_rev.hpp"
#include "AASeqConverter.hpp"
#include "Options.hpp"
#include "EnergyModel.hpp"
#include "ViennaEnergyModel.hpp"

//extern "C" {
//#include "params.h"       // not included through previous includes - keep
//#include "energy_const.h" // included through CDSfold_rev.hpp - can remove
//}

using namespace std;

struct Options;
          
using stack = struct stack {
    int i;
    int j;
    int Li;
    int Rj;
    int ml;
};

using bond = struct bond {
    int i;
    int j;
};

inline void clear_sec_bp(stack *s, bond *b, int len) {
    for (int i = 0; i < 500; i++) {
        s[i].i = -INF;
        s[i].j = -INF;
        s[i].Li = -INF;
        s[i].Rj = -INF;
        s[i].ml = -INF;
    }
    for (int i = 0; i < len / 2; i++) {
        b[i].i = -INF;
        b[i].j = -INF;
    }
}


static constexpr array<int, 20> make_i2r() {
    return {0, 1, 2, 3, 4, 4, 4, 3, 3, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // n[5] = 4; // 2nd position of L is converted to U
    // n[7] = 3; // 2nd position of R is converted to G
}


inline array<int, 100> make_ii2r() {
    array<int, 100> n;
    int s = 1;
    for (int i1 = 1; i1 <= 8; i1++) {
        for (int i2 = 1; i2 <= 8; i2++) {
            n[i1 * 10 + i2] = s++;
        }
    }
    return n;
}


inline map<char, int> make_n2i() {
    map<char, int> m;
    m['A'] = 1;
    m['C'] = 2;
    m['G'] = 3;
    m['U'] = 4;
    m['V'] = 5; // 2nd position of L, before A/G
    m['W'] = 6; // 2nd position of L, before U/C
    m['X'] = 7; // 2nd position of R, before A/G
    m['Y'] = 8; // 2nd position of R, before U/C

    return m;
}


static constexpr array<char, 20> make_i2n() {
    // array<char, 20> n = {};
    return {' ', 'A', 'C', 'G', 'U', 'V', 'W', 'X', 'Y', '\0',
            '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    // n[0] = ' ';
    // n[1] = 'A';
    // n[2] = 'C';
    // n[3] = 'G';
    // n[4] = 'U';
    // n[5] = 'V';
    // n[6] = 'W';
    // n[7] = 'X';
    // n[8] = 'Y';
    // return n;
}




vector<vector<int>> getPossibleNucleotide(std::string const & aaseq, codon &codon_table, map<char, int> const & n2i,
                                          string const & excludedCodons);

class Problem {
public:
    Problem(Options const & options, std::string const & aaseq);
    ~Problem() {
        // free(P_);
    }
    void calculate();
    
private:

    unique_ptr<EnergyModel> energyModel_;        /* energy model TODO change to unique ptr to access derived methods */
    Options options_;
    std::string aaseq_;
    unsigned int aalen_ = 0;
    unsigned int nuclen_ = 0;
    int max_bp_distance_final_ = 0;
    /* replace with a struct from the EnergyModel class; no longer directly access - must be through class */
    // shared_ptr<EnergyParams> P_ = nullptr;
    vector<vector<int>> Dep1_;
    vector<vector<int>> Dep2_;
    map<string, int> predefHPN_E_;
    vector<int> NucConst_;
    vector<vector<int>> pos2nuc_;
    vector<vector<vector<string>>> substr_;
    int n_inter_ = 0; // In the current implementation, n_inter = 1 or 2.
    int ofm_[5]; // I don't even think you need this much.
    int oto_[5]; // these have two different meanings for the max and min algs?
    vector<int> indx_;

    vector<vector<vector<int>>> C_, M_, F_, F2_;
    vector<array<array<int, 4>, 4>> DMl_, DMl1_, DMl2_;
    vector<int> ChkC_, ChkM_;
    vector<bond> base_pair_;

    // In theory it's a little silly to construct this once per sequence
    // but look, it's cheap and maybe it'll be helpful long term to have
    // this be a one to one relationship
    static const map<char, int> n2i;// = make_n2i();
    static const array<char, 20> i2n;// = make_i2n();
    static const array<int, 20> i2r;// = make_i2r();
    static const array<int, 100> ii2r;
    static constexpr char NucDef[] =
    "*AUGGAGGGGAUUGUCACGGGAGAUCGGCUUGCUUGCGUGGCGCUUCAUGGAAGCUCUUUGCUCCAUGAAGCGUCCGUAAGCAAGUAUACCGAUAUCCCGGGCAUUCUCC"
    "UCCAAUACAUCGAUGAAUUUCCCCUCACUGAUAUUGCCGCGCACGCGCCACGCGAGGCGUGGCAAAGCCUGUGCGAACAGGCGAUCUGUAUCGUCCAUCAUAUUAGCGAC"
    "CGGGGCAUCCUCAAUGAGGAUGUUAAAACCCGGUCGCUGACGAUACAGAUCAACAGUGAGGGGAUGUUCAAGAUGUUUAUG";
    static constexpr char dummy_str[10] = "XXXXXXXXX";

    // make_i2r(i2r);

    AASeqConverter conv;

    codon codon_table_;

    inline vector<int> set_ij_indx();


    auto rev_fold_step1() -> string;
    void rev_fold_step2(string & optseq_r);

    void fixed_fold(string & optseq);

    void allocate_arrays();
    void fill_F();
    void allocate_F2();
    void fill_F2();



    void backtrack(string *optseq, array<stack, 500> & sector, const int &initL, const int &initR);
    void backtrack2(string *optseq, array<stack, 500> & sector, const int &initL, const int &initR);
    void fixed_backtrack(string const & optseq, bond *base_pair, vector<int> const & c, vector<int> const & m, int *f, int nuclen, const int (&BP_pair)[5][5]);

};


inline auto getMatrixSize(int len, int w) -> int {
    int size = 0;
    for (int i = 1; i <= w; i++) {
        size += len - (i - 1); // Image of adding the number of diagonal elements of the matrix
                               // When i = 1, add diagonal elements (len).
    }

    cout << "The size of matrix is " << size << endl;
    return size;
}


inline auto getIndx(int const &i, int const &j, int const &w, vector<int> const &indx) -> int {
    return indx[j] + i - MAX2(0, j - w); // j-w is the number of unused elements.
                                         // If w is not specified (= length), j elements (1 <i <j) are prepared in the j column.
                                         // When w is specified, the number of elements used in column j is w, and the number of unused elements is j-w.
}

void allocate_F2(int len, vector<int> const & indx, int w, vector<vector<int>> &pos2nuc, vector<vector<vector<int>>> & f2);

inline void showChkMatrix(int *&m, vector<int > const &indx, int len, int w) {
    printf("REF:");
    for (int j = 1; j <= len; j++) {
        printf("\t%d", j);
    }
    printf("\n");

    for (int i = 1; i <= len; i++) {
        printf("REF:");
        printf("%d", i);
        for (int j = 1; j <= len; j++) {
            if (i >= j) {
                printf("\t-");
            } else {
                // printf("\t%d", m[indx[j]+i]);
                printf("\t%d", m[getIndx(i, j, w, indx)]);
            }
        }
        printf("\n");
    }
}


inline void showFixedMatrix(const int *m, vector<int> const & indx, const int len, const int w) {
    printf("REF:");
    for (int j = 1; j <= len; j++) {
        printf("\t%d", j);
    }
    printf("\n");

    for (int i = 1; i <= len; i++) {
        printf("REF:");
        printf("%d", i);
        for (int j = 1; j <= len; j++) {
            if (i >= j) {
                printf("\t-");
            } else {
                // printf("\t%d", m[indx[j]+i]);
                printf("\t%d", m[getIndx(i, j, w, indx)]);
            }
        }
        printf("\n");
    }
}

inline vector<int> createNucConstraint(const char *s, unsigned int const len, map<char, int> const & n2i) {
    vector<int> v(len + 1);
    for (int i = 1; i <= len; i++) {
        v[i] = n2i.at(s[i]);
    }
    return v;
}

inline vector<int> Problem::set_ij_indx() {
    vector<int> a( nuclen_ + 1, 0 );
    if (max_bp_distance_final_ <= 0) {
        cerr << "Invalid w:" << max_bp_distance_final_ << endl;
        exit(1);
    }
    max_bp_distance_final_ = MIN2(nuclen_, max_bp_distance_final_);
    int cum = 0;
    for (int n = 1; n <= nuclen_; n++) {
        a[n] = cum;
        // cout << n << ":" << a[n] << endl;
        if (n < max_bp_distance_final_) {
            cum += n;
        } else {
            cum += max_bp_distance_final_;
        }
    }
    return a;
}

auto getMemoryUsage(const string &fname) -> int;

inline void fixed_init_matrix(const int &nuclen, const int &size, vector<int> & C, vector<int> & M, vector<int> & F, vector<int> & DMl, vector<int> & DMl1, vector<int> & DMl2) {
    for (int i = 0; i <= nuclen; i++) {
        F[i] = 0;
        DMl[i] = INF;
        DMl1[i] = INF;
        DMl2[i] = INF;
    }

    for (int i = 0; i < size; i++) {
        C[i] = INF;
        M[i] = INF;
    }
}
