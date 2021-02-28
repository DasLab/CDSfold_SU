
#include <iostream>


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

#include "backtracking.hpp"
#include "CDSfold.hpp" // for getIndx, TermAU...

void backtrack(string *optseq, array<stack, 500> & sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f,
               vector<int> const & indx, const int &initL, const int &initR, paramT *const P, const vector<int> &NucConst,
               const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
               int const (&BP_pair)[5][5], array<char, 20> const &i2n, int *const &rtype, array<int, 100> const &ii2r,
               vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, int &DEPflg, map<string, int> &predefE,
               vector<vector<vector<string>>> &substr, map<char, int> &n2i, const char *nucdef) {

    int s = 0;
    int b = 0;
    sector[++s].i = 1;
    sector[s].j = length;
    sector[s].Li = initL;
    sector[s].Rj = initR;
    sector[s].ml = 0;

OUTLOOP:
    while (s > 0) {
        int fij, fi, ij, cij, traced, traced_Lk, i1, j1, k, p, q;
        //	    int canonical = 1;     /* (i,j) closes a canonical structure */

        //The value should be reflected in optseq here？
        int i = sector[s].i;
        int j = sector[s].j;
        int Li = sector[s].Li;
        int Rj = sector[s].Rj;
        int ml = sector[s--].ml; /* ml is a flag indicating if backtracking is to
                                    occur in the M- (1) or in the F-array (0) */
        int Li_nuc = pos2nuc[i][Li];
        int Rj_nuc = pos2nuc[j][Rj];

        if (i + 1 == j && Dep1[ii2r[Li_nuc * 10 + Rj_nuc]][i] == 0) {
            continue;
        }
        if (i + 2 == j && Dep2[ii2r[Li_nuc * 10 + Rj_nuc]][i] == 0) {
            continue;
        }

        int type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];

        (*optseq)[i] = i2n[Li_nuc];
        (*optseq)[j] = i2n[Rj_nuc];

        if (ml == 2) {
            base_pair[++b].i = i;
            base_pair[b].j = j;
            goto repeat1;
        }

        //	    if (j < i+TURN+1) continue; /* no more pairs in this interval */
        if (j == i)
            break;

        //	    fij = (ml == 1)? m[indx[j]+i][Li][Rj] : f[j][Li][Rj];
        fij = (ml == 1) ? m[getIndx(i, j, w, indx)][Li][Rj] : f[j][Li][Rj];
        //	    if(TB_CHK_flg == 1)
        cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << fij << ")"
             << ":" << *optseq << endl;

        for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
            int Rj1_nuc = pos2nuc[j - 1][Rj1];
            if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                continue;
            }

            if (Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                continue;
            }

            //	    	fi  = (ml == 1)? m[indx[j-1]+i][Li][Rj1] + P->MLbase: f[j-1][Li][Rj1];
            fi = (ml == 1) ? m[getIndx(i, j - 1, w, indx)][Li][Rj1] + P->MLbase : f[j - 1][Li][Rj1];

            if (fij == fi) { /* 3' end is unpaired */
                sector[++s].i = i;
                sector[s].j = j - 1;
                sector[s].Li = Li;
                sector[s].Rj = Rj1;
                sector[s].ml = ml;
                // continue;
                goto OUTLOOP;
            }
        }

        if (ml == 0) { /* backtrack in f */

            if (i != 1) {
                cerr << "Traceback failure: i must be 1 during bachtrack in f" << endl;
            }

            // f[j]とC[1][j]が一致している時の処理。Vieenaでは、次for文に統合されている。
            if (type_LiRj && j <= w) {
                // note i == 1
                // int en_c = TermAU(type_LiRj, P) +  c[indx[j]+i][Li][Rj];
                int en_c = TermAU(type_LiRj, P) + c[getIndx(i, j, w, indx)][Li][Rj];
                int en_f = f[j][Li][Rj];
                if (en_c == en_f) {
                    k = i;
                    traced = j;
                    traced_Lk = Li;
                    goto LABEL1;
                }
            }

            // for(k=j-TURN-1,traced=0; k>=1; k--){
            for (k = j - TURN - 1, traced = 0; k >= MAX2(2, j - w + 1); k--) {

                for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                    int Rk1_nuc = pos2nuc[k - 1][Rk1];
                    if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                        continue;
                    }
                    if (DEPflg && k == 3 && Dep1[ii2r[Li_nuc * 10 + Rk1_nuc]][i] == 0) {
                        continue;
                    } // dependency between 1(i) and 2(k-1)
                    if (DEPflg && k == 4 && Dep2[ii2r[Li_nuc * 10 + Rk1_nuc]][i] == 0) {
                        continue;
                    } // dependency between 1(i) and 3(k-1)

                    for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                        int Lk_nuc = pos2nuc[k][Lk];
                        if (NCflg == 1 && i2r[Lk_nuc] != NucConst[k]) {
                            continue;
                        }
                        if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                            continue;
                        } // dependency between k-1 and k

                        int type_LkRj = BP_pair[i2r[Lk_nuc]][i2r[Rj_nuc]];
                        if (type_LkRj) {
                            //	int en_c = TermAU(type_LkRj, P) +  c[indx[j]+k][Lk][Rj];
                            int en_c = TermAU(type_LkRj, P) + c[getIndx(k, j, w, indx)][Lk][Rj];
                            int en_f = f[k - 1][Li][Rk1];
                            //	                    	int test = en_c + en_f;
                            //	                    	cout << fij << "=" << test << "(" << en_c << "+" << en_f << ")"
                            //<< i << ":" << k << ":" << j << endl;
                            if (fij == en_c + en_f) {
                                traced = j;
                                traced_Lk = Lk;
                                /* push back the remaining f portion */
                                sector[++s].i = i;
                                sector[s].j = k - 1;
                                sector[s].Li = Li;
                                sector[s].Rj = Rk1;
                                sector[s].ml = 0;

                                goto LABEL1;
                            }
                        }
                    }
                }
            }
        LABEL1:

            if (!traced) {
                fprintf(stderr, "backtrack failed in f\n");
                fprintf(stderr, "cannot trace f[%d][%d][%d] Lnuc=%c Rnuc=%c \n", j, Li, Rj, i2n[Li_nuc], i2n[Rj_nuc]);
                exit(0);
            }

            /* trace back the base pair found */
            // [1]
            i = k;          // iを更新
            j = traced;     // この代入は多分必要ない。jに関しては不変だから
            Li = traced_Lk; // Liを更新
            // Rjは不変
            base_pair[++b].i = i;
            base_pair[b].j = j;
            goto repeat1;
        } else { /* trace back in fML array */

            for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                int Li1_nuc = pos2nuc[i + 1][Li1];
                if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                    continue;
                }
                if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                    continue;
                } // dependency between k-1 and k

                //	  	if (m[indx[j]+i+1][Li1][Rj]+P->MLbase == fij) { /* 5' end is unpaired */
                if (m[getIndx(i + 1, j, w, indx)][Li1][Rj] + P->MLbase == fij) { /* 5' end is unpaired */
                    sector[++s].i = i + 1;
                    sector[s].j = j;
                    sector[s].Li = Li1;
                    sector[s].Rj = Rj;
                    sector[s].ml = ml;
                    goto OUTLOOP;
                }

                //	ij  = indx[j]+i;
                ij = getIndx(i, j, w, indx);

                if (fij == c[ij][Li][Rj] + TermAU(type_LiRj, P) + P->MLintern[type_LiRj]) {
                    base_pair[++b].i = i;
                    base_pair[b].j = j;
                    goto repeat1;
                }
            }

            //		    for(k = i + 1 + TURN; k <= j - 2 - TURN; k++){
            for (k = i + 2 + TURN; k <= j - 1 - TURN; k++) {

                for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                    int Rk1_nuc = pos2nuc[k - 1][Rk1];
                    if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                        continue;
                    }
                    // check dependency is not needed because i+2<k,k+2<J
                    for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                        int Lk_nuc = pos2nuc[k][Lk];
                        if (NCflg == 1 && i2r[Lk_nuc] != NucConst[k]) {
                            continue;
                        }
                        if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                            continue;
                        } // dependency between k-1 and k

                        // if(fij == (m[indx[k-1]+i][Li][Rk1]+m[indx[j]+k][Lk][Rj])){
                        if (fij == (m[getIndx(i, k - 1, w, indx)][Li][Rk1] + m[getIndx(k, j, w, indx)][Lk][Rj])) {
                            sector[++s].i = i;
                            sector[s].j = k - 1;
                            sector[s].Li = Li;
                            sector[s].Rj = Rk1;
                            sector[s].ml = ml;
                            sector[++s].i = k;
                            sector[s].j = j;
                            sector[s].Li = Lk;
                            sector[s].Rj = Rj;
                            sector[s].ml = ml;
                            goto OUTLOOP;
                        }
                    }
                }
            }
            if (k > j - 1 - TURN) {
                fprintf(stderr, "backtrack failed in fML\n");
                exit(1);
            }
        }

    repeat1: // いちいちスタックに積まずに、ここで部分的なトレースバックをしてしまう。
        // continue;
        /*----- begin of "repeat:" -----*/
        if (j - i + 1 > w) {
            cerr << "backtrack failed at " << i << "," << j << " : the length must at most << w << endl";
        }
        //	    ij = indx[j]+i; // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
        ij = getIndx(i, j, w, indx); // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
        Li_nuc = pos2nuc[i][Li]; // Liは更新されている。
        Rj_nuc = pos2nuc[j][Rj]; // Rj_nucは更新されていない場合もある[1]。
        type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];
        cij = c[ij][Li][Rj];
        (*optseq)[i] = i2n[Li_nuc]; //塩基対部分を記録
        (*optseq)[j] = i2n[Rj_nuc];

        //		if(TB_CHK_flg == 1)
        //			cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << cij << ")" << ":" << *optseq
        //<< endl;

        // predefinedなヘアピンのトレースバック
        //		if ((j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8) &&
        //   			cij == predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].first){
        //			string s1 = predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].second;
        //			for(unsigned int k = 0; k < s1.size(); k++){
        //				(*optseq)[i+k] = s1[k]; //塩基を記録
        //			}
        //			cout << "Predefined Hairpin " << s1 << " at " << i << ":" << j << endl;
        //			goto OUTLOOP;
        //   	}
        if (j - i + 1 == 5 || j - i + 1 == 6 || j - i + 1 == 8) {
            int l = j - i + 1;
            for (unsigned int s = 0; s < substr[i][l].size(); s++) {
                string hpn = substr[i][l][s];
                int hL_nuc = n2i[hpn[0]];
                int hL2_nuc = n2i[hpn[1]];
                int hR2_nuc = n2i[hpn[l - 2]];
                int hR_nuc = n2i[hpn[l - 1]];
                if (hL_nuc != i2r[Li_nuc])
                    continue;
                if (hR_nuc != i2r[Rj_nuc])
                    continue;

                if (NCflg == 1) {
                    string s1 = string(nucdef).substr(i, l);
                    if (hpn != s1)
                        continue;
                }

                if (DEPflg && Li_nuc > 4 && Dep1[ii2r[Li_nuc * 10 + hL2_nuc]][i] == 0) {
                    continue;
                } // Dependency is already checked.
                if (DEPflg && Rj_nuc > 4 && Dep1[ii2r[hR2_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                    continue;
                } // ただし、Li_nuc、Rj_nucがVWXYのときだけは、一つ内側との依存関係をチェックする必要がある。

                // predefinedなヘアピンとの比較
                if (predefE.count(hpn) > 0) {
                    if (c[ij][Li][Rj] == predefE[hpn]) {
                        cout << "Predefined Hairpin at " << i << "," << j << endl;
                        for (unsigned int k = 0; k < hpn.size(); k++) {
                            (*optseq)[i + k] = hpn[k]; //塩基を記録
                        }
                        goto OUTLOOP;
                    }

                } else { // 普通のヘアピンとの比較
                    //					int energy = HairpinE(j - i - 1, type_LiRj,
                    //							i2r[hL2_nuc], i2r[hR2_nuc],
                    //							"NNNNNNNNN");
                    int energy = E_hairpin(j - i - 1, type_LiRj, i2r[hL2_nuc], i2r[hR2_nuc], "NNNNNNNNN", P);

                    if (c[ij][Li][Rj] == energy) {
                        for (unsigned int k = 0; k < hpn.size(); k++) {
                            (*optseq)[i + k] = hpn[k]; //塩基を記録
                        }
                        goto OUTLOOP;
                    }
                }
            }

        } else {
            // 普通のヘアピンのトレースバック
            for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                int Li1_nuc = pos2nuc[i + 1][Li1];
                if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                    continue;
                }
                if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                    continue;
                } // dependency between i and i+1

                for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                    int Rj1_nuc = pos2nuc[j - 1][Rj1];
                    if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                        continue;
                    }
                    if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                        continue;
                    } // dependency between j-1 and j

                    // if (cij == HairpinE(j-i-1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN")){
                    if (cij == E_hairpin(j - i - 1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN", P)) {
                        (*optseq)[i + 1] = i2n[Li1_nuc]; //塩基対の内側のミスマッチ塩基を記録
                        (*optseq)[j - 1] = i2n[Rj1_nuc];
                        goto OUTLOOP;
                    } else {
                        continue;
                    }
                }
            }
        }
        // Hairpinに該当がなければ、Internal loopのトレースバック。もっとも手強い。
        for (p = i + 1; p <= MIN2(j - 2 - TURN, i + MAXLOOP + 1); p++) {
            for (unsigned int Lp = 0; Lp < pos2nuc[p].size(); Lp++) {
                int Lp_nuc = pos2nuc[p][Lp];
                if (NCflg == 1 && i2r[Lp_nuc] != NucConst[p]) {
                    continue;
                }
                if (DEPflg && p == i + 1 && Dep1[ii2r[Li_nuc * 10 + Lp_nuc]][i] == 0) {
                    continue;
                } // dependency between i and q
                if (DEPflg && p == i + 2 && Dep2[ii2r[Li_nuc * 10 + Lp_nuc]][i] == 0) {
                    continue;
                } // dependency between i and q

                int minq = j - i + p - MAXLOOP - 2;
                if (minq < p + 1 + TURN)
                    minq = p + 1 + TURN;
                for (q = j - 1; q >= minq; q--) {
                    for (unsigned int Rq = 0; Rq < pos2nuc[q].size(); Rq++) {
                        int Rq_nuc = pos2nuc[q][Rq];
                        if (NCflg == 1 && i2r[Rq_nuc] != NucConst[q]) {
                            continue;
                        }
                        if (DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc * 10 + Rj_nuc]][q] == 0) {
                            continue;
                        } // dependency between q and j
                        if (DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc * 10 + Rj_nuc]][q] == 0) {
                            continue;
                        } // dependency between q and j

                        int type_LpRq = BP_pair[i2r[Lp_nuc]][i2r[Rq_nuc]];
                        if (type_LpRq == 0)
                            continue;
                        type_LpRq = rtype[type_LpRq];

                        for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                            int Li1_nuc = pos2nuc[i + 1][Li1];
                            if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                                continue;
                            }
                            if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                                continue;
                            } // dependency between i and i+1
                            if (i + 1 == p && Li1_nuc != Lp_nuc) {
                                continue;
                            } // i,pの時は、i+1の塩基とpの塩基は一致していないといけない。(1)

                            for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                                int Rj1_nuc = pos2nuc[j - 1][Rj1];
                                if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                    continue;
                                }
                                if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                                    continue;
                                } // dependency between j-1 and j
                                if (q == j - 1 && Rj1_nuc != Rq_nuc) {
                                    continue;
                                } // q,jの時は、qの塩基とj-1の塩基は一致していないといけない。(2)

                                for (unsigned int Lp1 = 0; Lp1 < pos2nuc[p - 1].size(); Lp1++) {
                                    int Lp1_nuc = pos2nuc[p - 1][Lp1];
                                    if (NCflg == 1 && i2r[Lp1_nuc] != NucConst[p - 1]) {
                                        continue;
                                    }

                                    if (DEPflg && Dep1[ii2r[Lp1_nuc * 10 + Lp_nuc]][p - 1] == 0) {
                                        continue;
                                    } // dependency between p-1 and p
                                    if (DEPflg && i == p - 2 && Dep1[ii2r[Li_nuc * 10 + Lp1_nuc]][i] == 0) {
                                        continue;
                                    } // i,X,p: dependency between i and p-1
                                    if (DEPflg && i == p - 3 && Dep1[ii2r[Li1_nuc * 10 + Lp1_nuc]][i + 1] == 0) {
                                        continue;
                                    } // i,X,X,p: dependency between i+1 and p-1

                                    if (i == p - 1 && Li_nuc != Lp1_nuc) {
                                        continue;
                                    } // i,pの時は、iの塩基とp-1の塩基は一致していないといけない。(1)の逆
                                    if (i == p - 2 && Li1_nuc != Lp1_nuc) {
                                        continue;
                                    } // i,X,pの時は、i+1の塩基とp-1の塩基(X)は一致していないといけない。

                                    for (unsigned int Rq1 = 0; Rq1 < pos2nuc[q + 1].size(); Rq1++) {
                                        int Rq1_nuc = pos2nuc[q + 1][Rq1];
                                        if (NCflg == 1 && i2r[Rq1_nuc] != NucConst[q + 1]) {
                                            continue;
                                        }
                                        if (DEPflg && Dep1[ii2r[Rq_nuc * 10 + Rq1_nuc]][q] == 0) {
                                            continue;
                                        } // dependency between q and q+1

                                        if (DEPflg && j == q + 2 && Dep1[ii2r[Rq1_nuc * 10 + Rj_nuc]][q + 1] == 0) {
                                            continue;
                                        } // q,X,j: dependency between j and q-1
                                        if (DEPflg && j == q + 3 && Dep1[ii2r[Rq1_nuc * 10 + Rj1_nuc]][q + 1] == 0) {
                                            continue;
                                        } // q,X,X,j: dependency between j+1 and q-1

                                        if (q + 1 == j && Rq1_nuc != Rj_nuc) {
                                            continue;
                                        } // q,jの時は、q+1の塩基とjの塩基は一致していないといけない。(2)の逆
                                        if (q + 2 == j && Rq1_nuc != Rj1_nuc) {
                                            continue;
                                        } // q,X,jの時は、q+1の塩基とj-1の塩基は一致していないといけない。

                                        // int energy = LoopEnergy(p-i-1, j-q-1, type_LiRj, type_LpRq,
                                        // 			i2r[Li1_nuc], i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc]);
                                        int energy = E_intloop(p - i - 1, j - q - 1, type_LiRj, type_LpRq, i2r[Li1_nuc],
                                                               i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc], P);

                                        //	int energy_new = energy+c[indx[q]+p][Lp][Rq];
                                        int energy_new = energy + c[getIndx(p, q, w, indx)][Lp][Rq];
                                        traced = (cij == energy_new);
                                        if (traced) {
                                            base_pair[++b].i = p;
                                            base_pair[b].j = q;

                                            (*optseq)[p] = i2n[Lp_nuc];
                                            (*optseq)[q] = i2n[Rq_nuc];

                                            (*optseq)[i + 1] = i2n[Li1_nuc];
                                            (*optseq)[p - 1] = i2n[Lp1_nuc];
                                            (*optseq)[j - 1] = i2n[Rj1_nuc];
                                            (*optseq)[q + 1] = i2n[Rq1_nuc];

                                            i = p, j = q; // i,jの更新
                                            Li = Lp;
                                            Rj = Rq;

                                            goto repeat1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /* end of repeat: --------------------------------------------------*/

        /* (i.j) must close a multi-loop */

        int rtype_LiRj = rtype[type_LiRj];
        i1 = i + 1;
        j1 = j - 1;

        sector[s + 1].ml = sector[s + 2].ml = 1;

        int en = cij - TermAU(rtype_LiRj, P) - P->MLintern[rtype_LiRj] - P->MLclosing;
        //	    for(k = i+2+TURN; k < j-2-TURN; k++){
        int Li1_save, Rk1_save, Lk_save, Rj1_save;
        Li1_save = Rk1_save = Lk_save = Rj1_save = -1;
        for (k = i + 3 + TURN; k < j - 1 - TURN; k++) {
            for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                int Rk1_nuc = pos2nuc[k - 1][Rk1];
                if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                    continue;
                }
                // i,k,jでの塩基の矛盾はチェックする必要がない。なぜなら、i+4<k, k+4<jだから。

                for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                    int Lk_nuc = pos2nuc[k][Lk];
                    if (NCflg == 1 && i2r[Lk_nuc] != NucConst[k]) {
                        continue;
                    }
                    if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                        continue;
                    } // dependency between k-1 and k

                    for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                        int Li1_nuc = pos2nuc[i + 1][Li1];
                        if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                            continue;
                        }
                        if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                            continue;
                        } // dependency between i and i+1

                        for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                            int Rj1_nuc = pos2nuc[j - 1][Rj1];
                            if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                continue;
                            }
                            if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                                continue;
                            } // dependency between j-1 and j

                            //マルチループを閉じるところと、bifucationを同時に探している。
                            // if(en == m[indx[k-1]+i+1][Li1][Rk1] + m[indx[j-1]+k][Lk][Rj1]){
                            if (en ==
                                m[getIndx(i + 1, k - 1, w, indx)][Li1][Rk1] + m[getIndx(k, j - 1, w, indx)][Lk][Rj1]) {
                                Li1_save = Li1;
                                Rk1_save = Rk1;
                                Lk_save = Lk;
                                Rj1_save = Rj1;
                                goto LABEL2;
                            }
                        }
                    }
                }
            }
        }
    LABEL2:

        if (k <= j - 2 - TURN) { /* found the decomposition successfully*/
            sector[++s].i = i1;
            sector[s].j = k - 1;
            sector[s].Li = Li1_save;
            sector[s].Rj = Rk1_save;

            sector[++s].i = k;
            sector[s].j = j1;
            sector[s].Li = Lk_save;
            sector[s].Rj = Rj1_save;

        } else {
            fprintf(stderr, "backtracking failed in repeat %d %d\n", i, j);
            //	    	exit(1);
        }
    }

    base_pair[0].i = b; /* save the total number of base pairs */
    //	cout << base_pair[0].i << endl;
}

void backtrack2(string *optseq, array<stack, 500> & sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f2,
                vector<int> const & indx, const int &initL, const int &initR, paramT *const P, const vector<int> &NucConst,
                const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
                int const (&BP_pair)[5][5], array<char, 20> const &i2n, int *const &rtype, array<int, 100> const &ii2r,
                vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, int &DEPflg, map<string, int> &predefE,
                vector<vector<vector<string>>> &substr, map<char, int> &n2i, const char *nucdef) {

    InitRand();

    int s = 0;
    int b = 0;
    sector[++s].i = 1;
    sector[s].j = length;
    sector[s].Li = initL;
    sector[s].Rj = initR;
    sector[s].ml = 0;

OUTLOOP:
    while (s > 0) {
        //	    int fij, fi, ij, cij, traced, traced_Lk, i1, j1, k, p , q;
        int fij, fi, ij, cij, traced, i1, j1, k, p, q;

        //ここでoptseqに値を反映させるべき？
        int i = sector[s].i;
        int j = sector[s].j;
        int Li = sector[s].Li;
        int Rj = sector[s].Rj;
        int ml = sector[s--].ml; /* ml is a flag indicating if backtracking is to
                                    occur in the M- (1) or in the F-array (0) */
        ij = getIndx(i, j, w, indx);
        int Li_nuc = pos2nuc[i][Li];
        int Rj_nuc = pos2nuc[j][Rj];

        if (i + 1 == j && Dep1[ii2r[Li_nuc * 10 + Rj_nuc]][i] == 0) {
            continue;
        }
        if (i + 2 == j && Dep2[ii2r[Li_nuc * 10 + Rj_nuc]][i] == 0) {
            continue;
        }

        int type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];

        (*optseq)[i] = i2n[Li_nuc];
        (*optseq)[j] = i2n[Rj_nuc];

        if (ml == 2) {
            base_pair[++b].i = i;
            base_pair[b].j = j;
            goto repeat1;
        }

        //	    if (j < i+TURN+1) continue; /* no more pairs in this interval */
        if (j == i + 1)
            continue;

        //	    fij = (ml == 1)? m[indx[j]+i][Li][Rj] : f[j][Li][Rj];
        fij = (ml == 1) ? m[getIndx(i, j, w, indx)][Li][Rj] : f2[ij][Li][Rj];
        cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << fij << ")"
             << ":" << *optseq << ":" << s << endl;

        // trace i,j from i,j-1 for multi-loop
        if (ml == 1) {
            for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                int Rj1_nuc = pos2nuc[j - 1][Rj1];
                if (Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                    continue;
                }
                int mi = m[getIndx(i, j - 1, w, indx)][Li][Rj1] + P->MLbase;

                if (fij == mi) { /* 3' end is unpaired */
                    sector[++s].i = i;
                    sector[s].j = j - 1;
                    sector[s].Li = Li;
                    sector[s].Rj = Rj1;
                    sector[s].ml = ml;
                    // continue;
                    goto OUTLOOP;
                }
            }
        }

        if (ml == 0) { /* backtrack in f */
            // traced = 0;
            int label[4];
            for (int l = 0; l < 4; l++) {
                label[l] = l;
            }
            shuffle(label, 4);

            for (int l : label) {
                cout << "go to label F" << label[l] + 1 << endl;
                if (l == 0) {
                    goto F1;
                } else if (l == 1) {
                    goto F2;
                } else if (l == 2) {
                    goto F3;
                } else {
                    goto F4;
                }

            F1:
                // trace i,j from i,j-1
                for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                    int Rj1_nuc = pos2nuc[j - 1][Rj1];
                    if (Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                        continue;
                    }

                    fi = f2[getIndx(i, j - 1, w, indx)][Li][Rj1];

                    if (fij == fi) { /* 3' end is unpaired */
                        sector[++s].i = i;
                        sector[s].j = j - 1;
                        sector[s].Li = Li;
                        sector[s].Rj = Rj1;
                        sector[s].ml = ml;
                        // continue;
                        cout << "Traceback path found." << endl;
                        goto OUTLOOP;
                    }
                }
                continue;

            F2:
                // trace i,j from i+1,j
                for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                    int Li1_nuc = pos2nuc[i + 1][Li1];
                    if (Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i + 1] == 0) {
                        continue;
                    }

                    fi = f2[getIndx(i + 1, j, w, indx)][Li1][Rj];

                    if (fij == fi) { /* 5' end is unpaired */
                        sector[++s].i = i + 1;
                        sector[s].j = j;
                        sector[s].Li = Li1;
                        sector[s].Rj = Rj;
                        sector[s].ml = ml;
                        // continue;
                        cout << "Traceback path found." << endl;
                        goto OUTLOOP;
                    }
                }
                continue;

            F3:
                // trace i,j from C(i,j)
                if (type_LiRj && j - i + 1 <= w) {
                    int en_c = TermAU(type_LiRj, P) + c[getIndx(i, j, w, indx)][Li][Rj];
                    int en_f = f2[ij][Li][Rj];
                    cout << en_c << "," << en_f << endl;
                    if (en_c == en_f) {
                        //	    			k = i;
                        //	    			traced = j;
                        base_pair[++b].i = i;
                        base_pair[b].j = j;
                        //	    			traced_Lk = Li;
                        //	    			goto LABEL1;
                        cout << "Traceback path found." << endl;
                        goto repeat1;
                    }
                }
                continue;

            F4:
                //	    		for(k=j-TURN-1; k>= i + TURN + 2; k--){
                //	    		for(k=j-1; k>= i + 2; k--){
                for (k = j - 1; k >= i + 2; k--) {

                    for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                        int Rk1_nuc = pos2nuc[k - 1][Rk1];
                        if (DEPflg && k == 3 && Dep1[ii2r[Li_nuc * 10 + Rk1_nuc]][i] == 0) {
                            continue;
                        } // dependency between 1(i) and 2(k-1)
                        if (DEPflg && k == 4 && Dep2[ii2r[Li_nuc * 10 + Rk1_nuc]][i] == 0) {
                            continue;
                        } // dependency between 1(i) and 3(k-1)

                        //	    				if((k - 1) - i + 1 > w ||
                        //	    					j - k + 1 > w)
                        //	    						continue;

                        for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                            int Lk_nuc = pos2nuc[k][Lk];
                            // if(NCflg == 1 && i2r[Lk_nuc] != NucConst[k]){continue;}
                            if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                                continue;
                            } // dependency between k-1 and k

                            int en_f1 = f2[getIndx(i, k - 1, w, indx)][Li][Rk1];
                            int en_f2 = f2[getIndx(k, j, w, indx)][Lk][Rj];
                            if (fij == en_f1 + en_f2) {
                                //	                    		traced = j;
                                //	                    		traced_Lk = Lk;
                                /* push back the remaining f portion */
                                sector[++s].i = i;
                                sector[s].j = k - 1;
                                sector[s].Li = Li;
                                sector[s].Rj = Rk1;
                                sector[s].ml = ml;
                                sector[++s].i = k;
                                sector[s].j = j;
                                sector[s].Li = Lk;
                                sector[s].Rj = Rj;
                                sector[s].ml = ml;
                                cout << "Traceback path found in " << i << "," << k - 1 << " and " << k << "," << j
                                     << endl;
                                goto OUTLOOP;
                            }
                        }
                    }
                }
                continue;
            }
            //	    	LABEL1:
            //	    	if (!traced){
            fprintf(stderr, "backtrack failed in f2\n");
            fprintf(stderr, "cannot trace f2[%d][%d][%d][%d] Lnuc=%c Rnuc=%c \n", i, j, Li, Rj, i2n[Li_nuc],
                    i2n[Rj_nuc]);
            exit(0);
            //	    	}

            /* trace back the base pair found */
            /*
            // [1]
            i=k;            // iを更新
            j=traced;       // この代入は多分必要ない。jに関しては不変だから
            Li = traced_Lk; // Liを更新
            //Rjは不変
            base_pair[++b].i = i;
            base_pair[b].j   = j;
            goto repeat1;
            */
        } else { /* trace back in fML array */

            for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                int Li1_nuc = pos2nuc[i + 1][Li1];
                if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                    continue;
                }
                if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                    continue;
                } // dependency between k-1 and k

                //	  	if (m[indx[j]+i+1][Li1][Rj]+P->MLbase == fij) { /* 5' end is unpaired */
                if (m[getIndx(i + 1, j, w, indx)][Li1][Rj] + P->MLbase == fij) { /* 5' end is unpaired */
                    sector[++s].i = i + 1;
                    sector[s].j = j;
                    sector[s].Li = Li1;
                    sector[s].Rj = Rj;
                    sector[s].ml = ml;
                    goto OUTLOOP;
                }

                //	ij  = indx[j]+i;
                ij = getIndx(i, j, w, indx);

                if (fij == c[ij][Li][Rj] + TermAU(type_LiRj, P) + P->MLintern[type_LiRj]) {
                    base_pair[++b].i = i;
                    base_pair[b].j = j;
                    goto repeat1;
                }
            }

            //		    for(k = i + 1 + TURN; k <= j - 2 - TURN; k++){
            for (k = i + 2 + TURN; k <= j - 1 - TURN; k++) {

                for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                    int Rk1_nuc = pos2nuc[k - 1][Rk1];
                    if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                        continue;
                    }
                    // check dependency is not needed because i+2<k,k+2<J
                    for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                        int Lk_nuc = pos2nuc[k][Lk];
                        if (NCflg == 1 && i2r[Lk_nuc] != NucConst[k]) {
                            continue;
                        }
                        if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                            continue;
                        } // dependency between k-1 and k

                        // if(fij == (m[indx[k-1]+i][Li][Rk1]+m[indx[j]+k][Lk][Rj])){
                        if (fij == (m[getIndx(i, k - 1, w, indx)][Li][Rk1] + m[getIndx(k, j, w, indx)][Lk][Rj])) {
                            sector[++s].i = i;
                            sector[s].j = k - 1;
                            sector[s].Li = Li;
                            sector[s].Rj = Rk1;
                            sector[s].ml = ml;
                            sector[++s].i = k;
                            sector[s].j = j;
                            sector[s].Li = Lk;
                            sector[s].Rj = Rj;
                            sector[s].ml = ml;
                            goto OUTLOOP;
                        }
                    }
                }
            }
            if (k > j - 1 - TURN) {
                fprintf(stderr, "backtrack failed in fML\n");
                exit(1);
            }
        }

    repeat1: // いちいちスタックに積まずに、ここで部分的なトレースバックをしてしまう。
        // continue;
        /*----- begin of "repeat:" -----*/
        if (j - i + 1 > w) {
            cerr << "backtrack failed at " << i << "," << j << " : the length must at most << w << endl";
        }
        //	    ij = indx[j]+i; // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
        ij = getIndx(i, j, w, indx); // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
        Li_nuc = pos2nuc[i][Li]; // Liは更新されている。
        Rj_nuc = pos2nuc[j][Rj]; // Rj_nucは更新されていない場合もある[1]。
        type_LiRj = BP_pair[i2r[Li_nuc]][i2r[Rj_nuc]];
        cij = c[ij][Li][Rj];
        (*optseq)[i] = i2n[Li_nuc]; //塩基対部分を記録
        (*optseq)[j] = i2n[Rj_nuc];

        //		if(TB_CHK_flg == 1)
        //			cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << cij << ")" << ":" << *optseq
        //<< endl;

        // predefinedなヘアピンのトレースバック
        //		if ((j-i+1 == 5 ||j-i+1 == 6 ||j-i+1 == 8) &&
        //   			cij == predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].first){
        //			string s1 = predefH[i][j-i+1][i2r[Li_nuc]][i2r[Rj_nuc]].second;
        //			for(unsigned int k = 0; k < s1.size(); k++){
        //				(*optseq)[i+k] = s1[k]; //塩基を記録
        //			}
        //			cout << "Predefined Hairpin " << s1 << " at " << i << ":" << j << endl;
        //			goto OUTLOOP;
        //   	}
        if (j - i + 1 == 5 || j - i + 1 == 6 || j - i + 1 == 8) {
            int l = j - i + 1;
            for (unsigned int s = 0; s < substr[i][l].size(); s++) {
                string hpn = substr[i][l][s];
                int hL_nuc = n2i[hpn[0]];
                int hL2_nuc = n2i[hpn[1]];
                int hR2_nuc = n2i[hpn[l - 2]];
                int hR_nuc = n2i[hpn[l - 1]];
                if (hL_nuc != i2r[Li_nuc])
                    continue;
                if (hR_nuc != i2r[Rj_nuc])
                    continue;

                if (NCflg == 1) {
                    string s1 = string(nucdef).substr(i, l);
                    if (hpn != s1)
                        continue;
                }

                if (DEPflg && Li_nuc > 4 && Dep1[ii2r[Li_nuc * 10 + hL2_nuc]][i] == 0) {
                    continue;
                } // Dependency is already checked.
                if (DEPflg && Rj_nuc > 4 && Dep1[ii2r[hR2_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                    continue;
                } // ただし、Li_nuc、Rj_nucがVWXYのときだけは、一つ内側との依存関係をチェックする必要がある。

                // predefinedなヘアピンとの比較
                if (predefE.count(hpn) > 0) {
                    if (c[ij][Li][Rj] == predefE[hpn]) {
                        cout << "Predefined Hairpin at " << i << "," << j << endl;
                        for (unsigned int k = 0; k < hpn.size(); k++) {
                            (*optseq)[i + k] = hpn[k]; //塩基を記録
                        }
                        goto OUTLOOP;
                    }

                } else { // Comparison with ordinary hairpins
                    //					int energy = HairpinE(j - i - 1, type_LiRj,
                    //							i2r[hL2_nuc], i2r[hR2_nuc],
                    //							"NNNNNNNNN");
                    int energy = E_hairpin(j - i - 1, type_LiRj, i2r[hL2_nuc], i2r[hR2_nuc], "NNNNNNNNN", P);

                    if (c[ij][Li][Rj] == energy) {
                        for (unsigned int k = 0; k < hpn.size(); k++) {
                            (*optseq)[i + k] = hpn[k]; // Record base
                        }
                        goto OUTLOOP;
                    }
                }
            }

        } else {
            // 普通のヘアピンのトレースバック
            for (int const Li1_nuc : pos2nuc[i + 1]) {
                if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                    continue;
                }
                if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                    continue;
                } // dependency between i and i+1

                for (int const Rj1_nuc : pos2nuc[j - 1]) {
                    if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                        continue;
                    }
                    if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                        continue;
                    } // dependency between j-1 and j

                    // if (cij == HairpinE(j-i-1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN")){
                    if (cij == E_hairpin(j - i - 1, type_LiRj, i2r[Li1_nuc], i2r[Rj1_nuc], "NNNNNNNNN", P)) {
                        (*optseq)[i + 1] = i2n[Li1_nuc]; // Record mismatched bases inside base pairs
                        (*optseq)[j - 1] = i2n[Rj1_nuc];
                        goto OUTLOOP;
                    } else {
                        continue;
                    }
                }
            }
        }
        // Hairpinに該当がなければ、Internal loopのトレースバック。もっとも手強い。
        for (p = i + 1; p <= MIN2(j - 2 - TURN, i + MAXLOOP + 1); p++) {
            for (unsigned int Lp = 0; Lp < pos2nuc[p].size(); Lp++) {
                int const Lp_nuc = pos2nuc[p][Lp];
                if (NCflg == 1 && i2r[Lp_nuc] != NucConst[p]) {
                    continue;
                }
                if (DEPflg && p == i + 1 && Dep1[ii2r[Li_nuc * 10 + Lp_nuc]][i] == 0) {
                    continue;
                } // dependency between i and q
                if (DEPflg && p == i + 2 && Dep2[ii2r[Li_nuc * 10 + Lp_nuc]][i] == 0) {
                    continue;
                } // dependency between i and q

                int minq = j - i + p - MAXLOOP - 2;
                if (minq < p + 1 + TURN)
                    minq = p + 1 + TURN;
                for (q = j - 1; q >= minq; q--) {
                    for (unsigned int Rq = 0; Rq < pos2nuc[q].size(); Rq++) {
                        int const Rq_nuc = pos2nuc[q][Rq];
                        if (NCflg == 1 && i2r[Rq_nuc] != NucConst[q]) {
                            continue;
                        }
                        if (DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc * 10 + Rj_nuc]][q] == 0) {
                            continue;
                        } // dependency between q and j
                        if (DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc * 10 + Rj_nuc]][q] == 0) {
                            continue;
                        } // dependency between q and j

                        int type_LpRq = BP_pair[i2r[Lp_nuc]][i2r[Rq_nuc]];
                        if (type_LpRq == 0)
                            continue;
                        type_LpRq = rtype[type_LpRq];

                        for (int const Li1_nuc : pos2nuc[i + 1]) {
                            if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                                continue;
                            }
                            if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                                continue;
                            } // dependency between i and i+1
                            if (i + 1 == p && Li1_nuc != Lp_nuc) {
                                continue;
                            } // In the case of i and p, the base of i + 1 and the base of p must match. (1)

                            for (int const Rj1_nuc : pos2nuc[j - 1]) {
                                if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                    continue;
                                }
                                if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                                    continue;
                                } // dependency between j-1 and j
                                if (q == j - 1 && Rj1_nuc != Rq_nuc) {
                                    continue;
                                } // In the case of q and j, the base of q and the base of j-1 must match. (2)

                                for (int const Lp1_nuc : pos2nuc[p - 1]) {
                                    if (NCflg == 1 && i2r[Lp1_nuc] != NucConst[p - 1]) {
                                        continue;
                                    }

                                    if (DEPflg && Dep1[ii2r[Lp1_nuc * 10 + Lp_nuc]][p - 1] == 0) {
                                        continue;
                                    } // dependency between p-1 and p
                                    if (DEPflg && i == p - 2 && Dep1[ii2r[Li_nuc * 10 + Lp1_nuc]][i] == 0) {
                                        continue;
                                    } // i,X,p: dependency between i and p-1
                                    if (DEPflg && i == p - 3 && Dep1[ii2r[Li1_nuc * 10 + Lp1_nuc]][i + 1] == 0) {
                                        continue;
                                    } // i,X,X,p: dependency between i+1 and p-1

                                    if (i == p - 1 && Li_nuc != Lp1_nuc) {
                                        continue;
                                    } // In the case of i and p, the base of i and the base of p-1 must match. The
                                      // reverse of (1)
                                    if (i == p - 2 && Li1_nuc != Lp1_nuc) {
                                        continue;
                                    } // In the case of i, X, p, the base of i + 1 and the base of p-1 (X) must match.

                                    for (int const Rq1_nuc : pos2nuc[q + 1]) {
                                        if (NCflg == 1 && i2r[Rq1_nuc] != NucConst[q + 1]) {
                                            continue;
                                        }
                                        if (DEPflg && Dep1[ii2r[Rq_nuc * 10 + Rq1_nuc]][q] == 0) {
                                            continue;
                                        } // dependency between q and q+1

                                        if (DEPflg && j == q + 2 && Dep1[ii2r[Rq1_nuc * 10 + Rj_nuc]][q + 1] == 0) {
                                            continue;
                                        } // q,X,j: dependency between j and q-1
                                        if (DEPflg && j == q + 3 && Dep1[ii2r[Rq1_nuc * 10 + Rj1_nuc]][q + 1] == 0) {
                                            continue;
                                        } // q,X,X,j: dependency between j+1 and q-1

                                        if (q + 1 == j && Rq1_nuc != Rj_nuc) {
                                            continue;
                                        } // q,jの時は、q+1の塩基とjの塩基は一致していないといけない。(2)の逆
                                        if (q + 2 == j && Rq1_nuc != Rj1_nuc) {
                                            continue;
                                        } // q,X,jの時は、q+1の塩基とj-1の塩基は一致していないといけない。

                                        // int energy = LoopEnergy(p-i-1, j-q-1, type_LiRj, type_LpRq,
                                        // 			i2r[Li1_nuc], i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc]);
                                        int energy = E_intloop(p - i - 1, j - q - 1, type_LiRj, type_LpRq, i2r[Li1_nuc],
                                                               i2r[Rj1_nuc], i2r[Lp1_nuc], i2r[Rq1_nuc], P);

                                        //	int energy_new = energy+c[indx[q]+p][Lp][Rq];
                                        int energy_new = energy + c[getIndx(p, q, w, indx)][Lp][Rq];
                                        traced = (cij == energy_new);
                                        if (traced) {
                                            base_pair[++b].i = p;
                                            base_pair[b].j = q;

                                            (*optseq)[p] = i2n[Lp_nuc];
                                            (*optseq)[q] = i2n[Rq_nuc];

                                            (*optseq)[i + 1] = i2n[Li1_nuc];
                                            (*optseq)[p - 1] = i2n[Lp1_nuc];
                                            (*optseq)[j - 1] = i2n[Rj1_nuc];
                                            (*optseq)[q + 1] = i2n[Rq1_nuc];

                                            i = p, j = q; // i,jの更新
                                            Li = Lp;
                                            Rj = Rq;

                                            goto repeat1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /* end of repeat: --------------------------------------------------*/

        /* (i.j) must close a multi-loop */

        int rtype_LiRj = rtype[type_LiRj];
        i1 = i + 1;
        j1 = j - 1;

        sector[s + 1].ml = sector[s + 2].ml = 1;

        int en = cij - TermAU(rtype_LiRj, P) - P->MLintern[rtype_LiRj] - P->MLclosing;
        //	    for(k = i+2+TURN; k < j-2-TURN; k++){
        int Li1_save, Rk1_save, Lk_save, Rj1_save;
        Li1_save = Rk1_save = Lk_save = Rj1_save = -1;
        for (k = i + 3 + TURN; k < j - 1 - TURN; k++) {
            for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                int Rk1_nuc = pos2nuc[k - 1][Rk1];
                if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                    continue;
                }
                // There is no need to check for base inconsistencies at i, k, j. Because i + 4 <k, k + 4 <j.

                for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                    int Lk_nuc = pos2nuc[k][Lk];
                    if (NCflg == 1 && i2r[Lk_nuc] != NucConst[k]) {
                        continue;
                    }
                    if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                        continue;
                    } // dependency between k-1 and k

                    for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                        int Li1_nuc = pos2nuc[i + 1][Li1];
                        if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                            continue;
                        }
                        if (DEPflg && Dep1[ii2r[Li_nuc * 10 + Li1_nuc]][i] == 0) {
                            continue;
                        } // dependency between i and i+1

                        for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                            int Rj1_nuc = pos2nuc[j - 1][Rj1];
                            if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                continue;
                            }
                            if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                                continue;
                            } // dependency between j-1 and j

                            // I'm looking for bifucation at the same time as closing the multi-loop.
                            // if(en == m[indx[k-1]+i+1][Li1][Rk1] + m[indx[j-1]+k][Lk][Rj1]){
                            if (en ==
                                m[getIndx(i + 1, k - 1, w, indx)][Li1][Rk1] + m[getIndx(k, j - 1, w, indx)][Lk][Rj1]) {
                                Li1_save = Li1;
                                Rk1_save = Rk1;
                                Lk_save = Lk;
                                Rj1_save = Rj1;
                                goto LABEL2;
                            }
                        }
                    }
                }
            }
        }
    LABEL2:

        if (k <= j - 2 - TURN) { /* found the decomposition successfully*/
            sector[++s].i = i1;
            sector[s].j = k - 1;
            sector[s].Li = Li1_save;
            sector[s].Rj = Rk1_save;

            sector[++s].i = k;
            sector[s].j = j1;
            sector[s].Li = Lk_save;
            sector[s].Rj = Rj1_save;

        } else {
            fprintf(stderr, "backtracking failed in repeat %d %d\n", i, j);
            //	    	exit(1);
        }
    }

    base_pair[0].i = b; /* save the total number of base pairs */
    //	cout << base_pair[0].i << endl;
}


void fixed_backtrack(string optseq, bond *base_pair, vector<int> const & c, vector<int> const & m, int *f, vector<int> const & indx, paramT * const P, int nuclen, int w,
                     const int (&BP_pair)[5][5], map<string, int> predefE) {
    int rtype[7] = {0, 2, 1, 4, 3, 6, 5};
    int s = 0;
    int b = 0;
    array<stack, 500> sector;
    sector[++s].i = 1;
    sector[s].j = nuclen;
    sector[s].ml = 0;

    map<char, int> n2i = make_n2i();
    vector<int> ioptseq;
	ioptseq.reserve(nuclen + 1);
    ioptseq[0] = 0;
    for (int i = 1; i <= nuclen; i++) {
        ioptseq[i] = n2i[optseq[i]];
    }

OUTLOOP:
    while (s > 0) {
        int fij, fi, ij, cij, traced, i1, j1, k, p, q;

        int i = sector[s].i;
        int j = sector[s].j;
        int ml = sector[s--].ml; /* ml is a flag indicating if backtracking is to
                                    occur in the M- (1) or in the F-array (0) */

        int type = BP_pair[ioptseq[i]][ioptseq[j]];

        if (ml == 2) {
            base_pair[++b].i = i;
            base_pair[b].j = j;
            goto repeat1;
        }

        if (j == i)
            break;

        fij = (ml == 1) ? m[getIndx(i, j, w, indx)] : f[j];
        cout << "TB_CHK:" << i << ":" << j << " " << ml << "(" << fij << ")" << endl;

        fi = (ml == 1) ? m[getIndx(i, j - 1, w, indx)] + P->MLbase : f[j - 1];

        if (fij == fi) { /* 3' end is unpaired */
            sector[++s].i = i;
            sector[s].j = j - 1;
            sector[s].ml = ml;
            // continue;
            goto OUTLOOP;
        }

        if (ml == 0) { /* backtrack in f */

            if (i != 1) {
                cerr << "Traceback failure: i must be 1 during bachtrack in f" << endl;
            }

            // f[j]とC[1][j]が一致している時の処理。Vieenaでは、次for文に統合されている。
            if (type && j <= w) {
                int en_c = TermAU(type, P) + c[getIndx(i, j, w, indx)];
                int en_f = f[j];
                if (en_c == en_f) {
                    k = i;
                    traced = j;
                    goto LABEL1;
                }
            }

            for (k = j - TURN - 1, traced = 0; k >= MAX2(2, j - w + 1); k--) {

                int type_kj = BP_pair[ioptseq[k]][ioptseq[j]];
                if (type_kj) {
                    int en_c = TermAU(type_kj, P) + c[getIndx(k, j, w, indx)];
                    int en_f = f[k - 1];
                    if (fij == en_c + en_f) {
                        traced = j;
                        /* push back the remaining f portion */
                        sector[++s].i = i;
                        sector[s].j = k - 1;
                        sector[s].ml = 0;

                        goto LABEL1;
                    }
                }
            }
        LABEL1:

            if (!traced) {
                fprintf(stderr, "backtrack failed in f\n");
                fprintf(stderr, "cannot trace f[%d] \n", j);
                exit(0);
            }

            /* trace back the base pair found */
            // [1]
            i = k;      // iを更新
            j = traced; // この代入は多分必要ない。jに関しては不変だから
            base_pair[++b].i = i;
            base_pair[b].j = j;
            goto repeat1;
        } else { /* trace back in fML array */

            if (m[getIndx(i + 1, j, w, indx)] + P->MLbase == fij) { /* 5' end is unpaired */
                sector[++s].i = i + 1;
                sector[s].j = j;
                sector[s].ml = ml;
                goto OUTLOOP;
            }

            ij = getIndx(i, j, w, indx);

            if (fij == c[ij] + TermAU(type, P) + P->MLintern[type]) {
                base_pair[++b].i = i;
                base_pair[b].j = j;
                goto repeat1;
            }

            for (k = i + 2 + TURN; k <= j - 1 - TURN; k++) {
                if (fij == (m[getIndx(i, k - 1, w, indx)] + m[getIndx(k, j, w, indx)])) {
                    sector[++s].i = i;
                    sector[s].j = k - 1;
                    sector[s].ml = ml;
                    sector[++s].i = k;
                    sector[s].j = j;
                    sector[s].ml = ml;
                    goto OUTLOOP;
                }
            }

            if (k > j - 1 - TURN) {
                fprintf(stderr, "backtrack failed in fML\n");
                exit(1);
            }
        }

    repeat1: // いちいちスタックに積まずに、ここで部分的なトレースバックをしてしまう。
        // continue;
        /*----- begin of "repeat:" -----*/
        if (j - i + 1 > w) {
            cerr << "backtrack failed at " << i << "," << j << " : the length must at most << w << endl";
        }
        ij = getIndx(i, j, w, indx); // ここでは元々のi,jから変換していることに注意。jは更新されないこともある[1]
        type = BP_pair[ioptseq[i]][ioptseq[j]];
        cij = c[ij];

        if (j - i + 1 == 5 || j - i + 1 == 6 || j - i + 1 == 8) {
            string hpn = "";
            for (int k = i; k <= j; k++) {
                hpn += optseq[k];
            }

            // predefinedなヘアピンとの比較
            if (predefE.count(hpn) > 0) {
                if (c[ij] == predefE[hpn]) {
                    cout << "Predefined Hairpin at " << i << "," << j << endl;
                    goto OUTLOOP;
                }
            } else { // 普通のヘアピンとの比較
                int energy = E_hairpin(j - i - 1, type, ioptseq[i + 1], ioptseq[j - 1], "NNNNNNNNN", P);

                if (c[ij] == energy) {
                    goto OUTLOOP;
                }
            }

        } else {
            // 普通のヘアピンのトレースバック
            if (cij == E_hairpin(j - i - 1, type, ioptseq[i + 1], ioptseq[j - 1], "NNNNNNNNN", P)) {
                goto OUTLOOP;
            }
        }
        // Hairpinに該当がなければ、Internal loopのトレースバック。もっとも手強い。
        for (p = i + 1; p <= MIN2(j - 2 - TURN, i + MAXLOOP + 1); p++) {
            int minq = j - i + p - MAXLOOP - 2;
            if (minq < p + 1 + TURN)
                minq = p + 1 + TURN;
            for (q = j - 1; q >= minq; q--) {

                int type_pq = BP_pair[ioptseq[p]][ioptseq[q]];
                if (type_pq == 0)
                    continue;
                type_pq = rtype[type_pq];

                int energy = E_intloop(p - i - 1, j - q - 1, type, type_pq, ioptseq[i + 1], ioptseq[j - 1],
                                       ioptseq[p - 1], ioptseq[q + 1], P);

                int energy_new = energy + c[getIndx(p, q, w, indx)];
                traced = (cij == energy_new);
                if (traced) {
                    base_pair[++b].i = p;
                    base_pair[b].j = q;

                    i = p, j = q; // i,jの更新
                    goto repeat1;
                }
            }
        }
        /* end of repeat: --------------------------------------------------*/

        /* (i.j) must close a multi-loop */

        int type_rev = rtype[type];
        i1 = i + 1;
        j1 = j - 1;

        sector[s + 1].ml = sector[s + 2].ml = 1;

        int en = cij - TermAU(type_rev, P) - P->MLintern[type_rev] - P->MLclosing;
        for (k = i + 3 + TURN; k < j - 1 - TURN; k++) {
            //マルチループを閉じるところと、bifucationを同時に探している。
            if (en == m[getIndx(i + 1, k - 1, w, indx)] + m[getIndx(k, j - 1, w, indx)]) {
                goto LABEL2;
            }
        }
    LABEL2:

        if (k <= j - 2 - TURN) { /* found the decomposition successfully*/
            sector[++s].i = i1;
            sector[s].j = k - 1;

            sector[++s].i = k;
            sector[s].j = j1;

        } else {
            fprintf(stderr, "fixed_backtracking failed in repeat %d %d\n", i, j);
            //	    	exit(1);
        }
    }

    base_pair[0].i = b; /* save the total number of base pairs */
    //	cout << base_pair[0].i << endl;
}