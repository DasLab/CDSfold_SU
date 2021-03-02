#include <map>
#include <cstring>
#include <fstream>
#include <sstream>
#include "codon.hpp"
#include "backtracking.hpp" // just for bond
#include "energy.hpp"

using namespace std;

#pragma once

vector<vector<int>> getPossibleNucleotide(std::string const & aaseq, codon &codon_table, map<char, int> const & n2i,
                                          string const & excludedCodons);

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


inline void allocate_arrays(
    int len,
    vector<int> const & indx,
    int w,
    vector<vector<int>> &pos2nuc,
    vector<vector<vector<int>>> & c,
    vector<vector<vector<int>>> & m,
    vector<vector<vector<int>>> & f,
    vector<array<array<int, 4>, 4>> & dml,
    vector<array<array<int, 4>, 4>> & dml1,
    vector<array<array<int, 4>, 4>> & dml2,
    vector<int> & chkc,
    vector<int> & chkm,
    vector<bond> & b
) {
    int size = getMatrixSize(len, w);
    c.resize(size + 1);
    m.resize(size + 1);
    for (int i = 1; i <= len; i++) {
        for (int j = i; j <= MIN2(len, i + w - 1); j++) {
            // cout << i << " " << j << endl;
            int ij = getIndx(i, j, w, indx);
            c[ij].resize(pos2nuc[i].size());
            m[ij].resize(pos2nuc[i].size());
            for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                c[ij][L].resize(pos2nuc[j].size());
                m[ij][L].resize(pos2nuc[j].size());
            }
        }
    }

    f.resize(len + 1);
    dml.resize(len + 1, {{{INF, INF, INF, INF},{INF, INF, INF, INF},{INF, INF, INF, INF},{INF, INF, INF, INF}}});
    dml1.resize(len + 1, {{{INF, INF, INF, INF},{INF, INF, INF, INF},{INF, INF, INF, INF},{INF, INF, INF, INF}}});
    dml2.resize(len + 1, {{{INF, INF, INF, INF},{INF, INF, INF, INF},{INF, INF, INF, INF},{INF, INF, INF, INF}}});
    for (int j = 1; j <= len; j++) {
        f[j].resize(pos2nuc[j].size());
        for (unsigned int L = 0; L < pos2nuc[1].size(); L++) { // The first position
            f[j][L].resize(pos2nuc[j].size());
        }
    }

    chkc.resize(size + 1, INF);
    chkm.resize(size + 1, INF);

    b.resize(len / 2);
}


inline void allocate_F2(int len, vector<int> const & indx, int w, vector<vector<int>> &pos2nuc, vector<vector<vector<int>>> & f2) {
    int size = getMatrixSize(len, w);
    f2.resize(size + 1);
    for (int i = 1; i <= len; i++) {
        for (int j = i; j <= MIN2(len, i + w - 1); j++) {
            int ij = getIndx(i, j, w, indx);
            f2[ij].resize(pos2nuc[i].size());
            for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                f2[ij][L].resize(pos2nuc[j].size());
            }
        }
    }
}


inline vector<int> set_ij_indx(int length) {
    vector<int> a( length + 1, 0 );
    for (int n = 1; n <= length; n++) {
        a[n] = (n * (n - 1)) / 2;
    }
    return a;
}


inline vector<int> set_ij_indx(int length, int w) {
    vector<int> a( length + 1, 0 );
    if (w <= 0) {
        cerr << "Invalid w:" << w << endl;
        exit(1);
    }
    w = MIN2(length, w);
    int cum = 0;
    for (int n = 1; n <= length; n++) {
        a[n] = cum;
        // cout << n << ":" << a[n] << endl;
        if (n < w) {
            cum += n;
        } else {
            cum += w;
        }
    }
    return a;
}


inline void make_i2r(array<int, 20> & n) {
    n[1] = 1;
    n[2] = 2;
    n[3] = 3;
    n[4] = 4;
    n[5] = 4; // 2nd position of L is converted to U
    n[6] = 4;
    n[7] = 3; // 2nd position of R is converted to G
    n[8] = 3;
}


inline void make_ii2r(array<int, 100> & n) {
    int s = 1;
    for (int i1 = 1; i1 <= 8; i1++) {
        for (int i2 = 1; i2 <= 8; i2++) {
            n[i1 * 10 + i2] = s++;
        }
    }
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

// void view_n2i(map<char, int> n2i, char const &c){

inline auto view_n2i(map<char, int> n2i, char c) -> int {
    // cout << c << endl;
    return n2i[c];
}


inline void make_i2n(array<char, 20> & n) {
    n[0] = ' ';
    n[1] = 'A';
    n[2] = 'C';
    n[3] = 'G';
    n[4] = 'U';
    n[5] = 'V';
    n[6] = 'W';
    n[7] = 'X';
    n[8] = 'Y';
}


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


inline vector<int> createNucConstraint(const char *s, unsigned int const len, map<char, int> &n2i) {
    vector<int> v(len + 1);
    for (int i = 1; i <= len; i++) {
        v[i] = n2i[s[i]];
    }
    return v;
}


inline void InitRand() { srand((unsigned int)time(nullptr)); }


inline void shuffleStr(vector<string>(*ary), int size) {
    for (int i = 0; i < size; i++) {
        int j = rand() % size;
        string t = (*ary)[i];
        (*ary)[i] = (*ary)[j];
        (*ary)[j] = t;
    }
}

inline void shuffle(int ary[], int size) {
    for (int i = 0; i < size; i++) {
        int j = rand() % size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }
}

inline vector<pair<int, int>> shufflePair(vector<pair<int, int>> ary, int size) {
    for (int i = 0; i < size; i++) {
        int j = rand() % size;
        // cout << rand() << ":" << j << endl;
        pair<int, int> t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }

    return ary;
}

inline auto getMemoryUsage(const string &fname) -> int {
    // cout << fname << endl;
    ifstream ifs(fname.c_str());
    if (ifs) {
        string line;
        while (getline(ifs, line)) { // 行の読み込み
            string index = line.substr(0, 6);
            if (index == "VmRSS:") {
                // cout << line << endl;
                string mem_str;
                for (unsigned int i = 6; i < line.size(); i++) {
                    if (line[i] == ' ' || line[i] == '\t')
                        continue;
                    if (isdigit(line[i])) {
                        mem_str.append(1, line[i]);
                    } else if (line[i] == 'k' || line[i] == 'B') {
                        continue;
                    } else {
                        cerr << "Unexpected letter found in " << fname << "(" << line[i] << ")" << endl;
                        return -1;
                    }
                }
                stringstream ss(mem_str);
                int val;
                ss >> val;
                return val;
            }
        }

    } else {
        cerr << "Error: cannot open file(" << fname << ")" << endl;
        return -1;
    }

    cerr << "VmRSS line was not found in " << fname << endl;
    return -1;
}


inline void fixed_init_matrix(const int &nuclen, const int &size, vector<int> & C, vector<int> & M, int *F, int *DMl, int *DMl1, int *DMl2) {
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


void fixed_fold(string optseq, vector<int> const & indx, const int &w, map<string, int> &predefE, const int (&BP_pair)[5][5],
                paramT *P, char const *aaseq, const codon& codon_table);