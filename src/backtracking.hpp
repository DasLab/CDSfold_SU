
#include <map>
// #include <string>
// #include <iostream>
#include <iosfwd>
#include <vector>

#pragma once

using namespace std;

           
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

void backtrack(string *optseq, array<stack, 500> & sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f2,
               vector<int> const & indx, const int &initL, const int &initR, paramT * const P, const vector<int> &NucConst,
               const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
               int const (&BP_pair)[5][5], array<char, 20> const &i2n, int *const &rtype, array<int, 100> const &ii2r,
               vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, int &DEPflg, map<string, int> &predefE,
               vector<vector<vector<string>>> &substr, map<char, int> &n2i, const char *nucdef);

void backtrack2(string *optseq, array<stack, 500> & sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f2,
                vector<int> const &, const int &initL, const int &initR, paramT * const P, const vector<int> &NucConst,
                const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
                int const (&BP_pair)[5][5], array<char, 20> const &i2n, int *const &rtype, array<int, 100> const &ii2r,
                vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, int &DEPflg, map<string, int> &predefE,
                vector<vector<vector<string>>> &substr, map<char, int> &n2i, const char *nucdef);
                               
void fixed_backtrack(string optseq, bond *base_pair, vector<int> const & c, vector<int> const & m, int *f, vector<int> const & indx, paramT *P, int nuclen, int w,
                     const int (&BP_pair)[5][5], map<string, int> predefE);