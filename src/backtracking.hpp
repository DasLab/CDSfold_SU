
#include <map>
// #include <string>
// #include <iostream>
#include <iosfwd>
#include <vector>

extern "C" {
#include <cctype>
#include "fold.h"
// #include "fold_vars.h"
#include <climits>
#include <cmath>
#include "params.h"
// #include "part_func.h"
#include <cstdio>
#include <cstdlib>
#include "utils.h"
}

#include "Problem.hpp" // for getIndx

#pragma once

using namespace std;


void backtrack(string *optseq, array<stack, 500> & sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f2,
               vector<int> const & indx, const int &initL, const int &initR, paramT * const P, const vector<int> &NucConst,
               const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
               int const (&BP_pair)[5][5], array<char, 20> const &i2n, int const *const &rtype, array<int, 100> const &ii2r,
               vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, bool const DEPflg, map<string, int> &predefE,
               vector<vector<vector<string>>> &substr, map<char, int> const & n2i, const char *nucdef);

void backtrack2(string *optseq, array<stack, 500> & sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f2,
                vector<int> const &, const int &initL, const int &initR, paramT * const P, const vector<int> &NucConst,
                const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
                int const (&BP_pair)[5][5], array<char, 20> const &i2n, int const *const &rtype, array<int, 100> const &ii2r,
                vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, bool const DEPflg, map<string, int> &predefE,
                vector<vector<vector<string>>> &substr, map<char, int> const & n2i, const char *nucdef);
                               
void fixed_backtrack(string optseq, bond *base_pair, vector<int> const & c, vector<int> const & m, int *f, vector<int> const & indx, paramT *P, int nuclen, int w,
                     const int (&BP_pair)[5][5], map<string, int> predefE);