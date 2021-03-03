#include <map>
#include <array>
#include <vector>
#include <string>
#include <unistd.h>

#include "Problem.hpp"
#include "Options.hpp"
#include "AASeqConverter.hpp"
#include "codon.hpp"
#include "energy.hpp"
#include "backtracking.hpp"

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

using namespace std;

static const int BP_pair[5][5] =
    /* _  A  C  G  U  */
    {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 5}, {0, 0, 0, 1, 0}, {0, 0, 2, 0, 3}, {0, 6, 0, 4, 0}};

static const int rtype[7] = {0, 2, 1, 4, 3, 6, 5};

const map<char, int> Problem::n2i = make_n2i();
const array<char, 20> Problem::i2n = make_i2n();
const array<int, 20> Problem::i2r = make_i2r();
const array<int, 100> Problem::ii2r = make_ii2r();

Problem::Problem(Options const & options, string const & aaseq):
    options_( options ),
    aaseq_( aaseq )
{
    unsigned int aalen = aaseq.size();

    if (aalen <= 2) {
        cerr << "The amino acid sequence is too short.\n";
        exit(1);
    }

    int n_inter = 0; // In the current implementation, n_inter = 1 or 2.
    int ofm[5]; // I don't even think you need this much.
    int oto[5]; // these have two different meanings for the max and min algs?
    if (options_.partial_opt) {
        if (options_.opt_to > aalen) {
            options_.opt_to = aalen;
        }

        // Creating partial inverse optimization information
        if (options_.maximize_mfe) { // Structural removal of the specified area
            ofm[0] = (options_.opt_from - 1) * 3 + 1;
            oto[0] = options_.opt_to * 3;
            n_inter = 1;
        } else { // Structural stabilization of the specified area
            int l = 0;
            if (options_.opt_from != 1) {
                ofm[l] = 1;
                oto[l++] = (options_.opt_from - 1) * 3;
            }
            if (options_.opt_to != aalen) {
                ofm[l] = options_.opt_to * 3 + 1;
                oto[l++] = aalen * 3;
            }
            n_inter = l;
        }
        // Checking the optimization area
        // for(int I = 0; I < n_inter; I++){
        //	cout << ofm[I] << "-" << oto[I] << endl;
        //}
        // exit(0);
    }

    unsigned int nuclen = aalen * 3;
    int max_bp_distance_final = 0;

    if (options_.max_bp_distance == 0) {
        max_bp_distance_final = nuclen;
    } else if (options_.max_bp_distance < 10) {
        cerr << "W must be more than 10"
                << "(you used " << options_.max_bp_distance << ")" << endl;
        exit(1);
    } else if (options_.max_bp_distance > nuclen) {
        max_bp_distance_final = nuclen;
    } else {
        max_bp_distance_final = options_.max_bp_distance;
    }

    //		w_tmp = 50;// test!
    //		vector<vector<vector<string> > >  substr = conv.getBases(string(aaseq),8, exc);
    vector<vector<vector<string>>> substr = conv.getOriginalBases(aaseq, options_.codons_excluded);
    vector<vector<int>> Dep1;
    vector<vector<int>> Dep2;
    float ptotal_Mb_base = 0;

    if (options_.estimate_memory_use) {
        ptotal_Mb_base = 2 + nuclen * 0.006956;
    } else {
        Dep1 = conv.countNeighborTwoBase(aaseq, options_.codons_excluded);
        Dep2 = conv.countEveryOtherTwoBase(aaseq, options_.codons_excluded);
    }

    map<string, int> predefHPN_E = conv.getBaseEnergy();

    // vector<int> NucConst = createNucConstraint(NucDef, nuclen, n2i);

    vector<int> NucConst;
    if (options_.nucleotide_constraints) {
        NucConst = createNucConstraint(NucDef, nuclen, n2i);
    }
    //		for(int i = 1; i <= nuclen; i++)
    //		printf("%d %d\n", i, NucConst[i]);

    // createNucConstraint

    cout << aaseq << endl;
    //		cout << aalen << endl;

    vector<vector<int>> pos2nuc = getPossibleNucleotide(aaseq, codon_table, n2i, options_.codons_excluded);
    //		vector<vector<int> > pos2nuc = getPossibleNucleotide(aaseq, aalen, codon_table, n2i, 'R');
    //		showPos2Nuc(pos2nuc, i2n);
    //		exit(0);

    if (options_.estimate_memory_use) {
        float ptotal_Mb_alloc = predict_memory(nuclen, max_bp_distance_final, pos2nuc);
        float ptotal_Mb = ptotal_Mb_alloc + ptotal_Mb_base;
        cout << "Estimated memory usage: " << ptotal_Mb << " Mb" << endl;
        exit(0);
    }

    paramT *P = scale_parameters();
    update_fold_params();

    //		rev_flg = 0;
    //		if(rev_flg && num_interval == 0){
    vector<int> const indx = set_ij_indx(nuclen, max_bp_distance_final);
    if (options_.maximize_mfe && !options_.partial_opt) {
        // reverse mode
        string optseq_rev = rev_fold_step1(&aaseq[0], aalen, codon_table, options_.codons_excluded);
        //			rev_fold_step2(&optseq_rev, aaseq, aalen, codon_table, exc, ofm, oto, 1);
        rev_fold_step2(&optseq_rev, &aaseq[0], aalen, codon_table, options_.codons_excluded);
        fixed_fold(optseq_rev, indx, max_bp_distance_final, predefHPN_E, BP_pair, P, &aaseq[0], codon_table);
        free(P);
        return; // returnすると、実行時間が表示されなくなるためbreakすること。
    }

    vector<vector<vector<int>>> C, M, F, F2;
    vector<array<array<int, 4>, 4>> DMl, DMl1, DMl2;
    vector<int> chkC, chkM;
    vector<bond> base_pair;
    allocate_arrays(nuclen, indx, max_bp_distance_final, pos2nuc, C, M, F, DMl, DMl1, DMl2, chkC, chkM, base_pair);
    if (options_.random_backtrack) {
        allocate_F2(nuclen, indx, max_bp_distance_final, pos2nuc, F2);
    }
    // float ptotal_Mb = ptotal_Mb_alloc + ptotal_Mb_base;

    //		pid_t pid1 = getpid();
    //		stringstream ss1;
    //		ss1 << "/proc/" << pid1 << "/status";
    //		int m1 = getMemoryUsage(ss1.str());
    //		cout << "Memory(VmRSS): "  << float(m1)/1024 << " Mb" << endl;
    //		exit(0);

    // main routine
    for (int l = 2; l <= 4; l++) {
        for (int i = 1; i <= nuclen - l + 1; i++) {
            //				test = 1;
            int j = i + l - 1;
            // int ij = indx[j] + i;
            int ij = getIndx(i, j, max_bp_distance_final, indx);

            chkC[ij] = INF;
            chkM[ij] = INF;

            for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                int L_nuc = pos2nuc[i][L];
                if (options_.nucleotide_constraints == 1 && i2r[L_nuc] != NucConst[i]) {
                    continue;
                }
                for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {
                    int R_nuc = pos2nuc[j][R];
                    if (options_.nucleotide_constraints == 1 && i2r[R_nuc] != NucConst[j]) {
                        continue;
                    }
                    // L-R pair must be filtered
                    //						if(j-1==1){
                    //							cout << i << ":" << j << " " << L_nuc << "-"
                    //<< R_nuc << " " << Dep1[ii2r[L_nuc*10+R_nuc]][i] <<endl;
                    //						}
                    if (options_.DEPflg && j - i == 1 && i <= nuclen - 1 && Dep1[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
                        continue;
                    } // nuclen - 1はいらないのでは？
                    if (options_.DEPflg && j - i == 2 && i <= nuclen - 2 && Dep2[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
                        continue;
                    }

                    C[ij][L][R] = INF;
                    M[ij][L][R] = INF;
                    if (options_.random_backtrack)
                        F2[ij][L][R] = 0;
                }
            }
        }
    }

    //		cout << "TEST" << M[13][0][0] << endl;
    // main routine
    for (int l = 5; l <= nuclen; l++) {
        if (l > max_bp_distance_final)
            break;
        cout << "process:" << l << endl;

        //	  for(int l = 5; l <= 5; l++){
        for (int i = 1; i <= nuclen - l + 1; i++) {
            int j = i + l - 1;

            int opt_flg_ij = 1;
            if (options_.partial_opt) {
                for (int I = 0; I < n_inter; I++) {
                    if ((ofm[I] <= i && oto[I] >= i) || (ofm[I] <= j && oto[I] >= j)) {
                        opt_flg_ij = 0;
                        break;
                    }
                }
            }

            for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                int L_nuc = pos2nuc[i][L];
                //					cout << NCflg << endl;
                if (options_.nucleotide_constraints == 1 && i2r[L_nuc] != NucConst[i]) {
                    continue;
                }
                //					cout << "ok" << endl;

                for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {

                    int R_nuc = pos2nuc[j][R];

                    if (options_.nucleotide_constraints == 1 && i2r[R_nuc] != NucConst[j]) {
                        continue;
                    }

                    // int ij = indx[j] + i;
                    int ij = getIndx(i, j, max_bp_distance_final, indx);

                    C[ij][L][R] = INF;
                    M[ij][L][R] = INF;
                    //						cout << i << " " << j << ":" << M[ij][L][R] <<
                    // endl;

                    int type = BP_pair[i2r[L_nuc]][i2r[R_nuc]];

                    if (!type || !opt_flg_ij) {
                        C[ij][L][R] = INF;
                    } else {
                        // hairpin
                        if ((l == 5 || l == 6 || l == 8)) {
                            for (string const & hpn : substr[i][l]) {
                                int hL_nuc = n2i.at(hpn[0]);
                                int hL2_nuc = n2i.at(hpn[1]);
                                int hR2_nuc = n2i.at(hpn[l - 2]);
                                int hR_nuc = n2i.at(hpn[l - 1]);
                                if (hL_nuc != i2r[L_nuc])
                                    continue;
                                if (hR_nuc != i2r[R_nuc])
                                    continue;

                                if (options_.nucleotide_constraints) {
                                    string s1 = string(NucDef).substr(i, l);
                                    if (hpn != s1)
                                        continue;
                                }

                                //									cout << hpn <<
                                // endl;
                                //
                                if (options_.DEPflg && L_nuc > 4 && Dep1[ii2r[L_nuc * 10 + hL2_nuc]][i] == 0) {
                                    continue;
                                } // Dependencyをチェックした上でsubstringを求めているので,
                                    // hpnの内部についてはチェックする必要はない。
                                if (options_.DEPflg && R_nuc > 4 && Dep1[ii2r[hR2_nuc * 10 + R_nuc]][j - 1] == 0) {
                                    continue;
                                } // ただし、L_nuc、R_nucがVWXYのときだけは、一つ内側との依存関係をチェックする必要がある。
                                    // その逆に、一つ内側がVWXYのときはチェックの必要はない。既にチェックされているので。
                                if (predefHPN_E.count(hpn) > 0) {
                                    C[ij][L][R] = MIN2(predefHPN_E[hpn], C[ij][L][R]);

                                } else {
                                    //										int energy
                                    //= HairpinE(j - i - 1, type,
                                    // i2r[hL2_nuc], i2r[hR2_nuc],
                                    // dummy_str);
                                    int energy =
                                        E_hairpin(j - i - 1, type, i2r[hL2_nuc], i2r[hR2_nuc], dummy_str, P);
                                    C[ij][L][R] = MIN2(energy, C[ij][L][R]);
                                }
                            }
                            // exit(0);
                        } else {
                            for (unsigned int L2 = 0; L2 < pos2nuc[i + 1].size(); L2++) {
                                int L2_nuc = pos2nuc[i + 1][L2];
                                if (options_.nucleotide_constraints == 1 && i2r[L2_nuc] != NucConst[i + 1]) {
                                    continue;
                                }
                                // if(chkDep2){continue:}
                                for (unsigned int R2 = 0; R2 < pos2nuc[j - 1].size(); R2++) {
                                    int R2_nuc = pos2nuc[j - 1][R2];
                                    if (options_.nucleotide_constraints == 1 && i2r[R2_nuc] != NucConst[j - 1]) {
                                        continue;
                                    }

                                    if (options_.DEPflg && Dep1[ii2r[L_nuc * 10 + L2_nuc]][i] == 0) {
                                        continue;
                                    }
                                    if (options_.DEPflg && Dep1[ii2r[R2_nuc * 10 + R_nuc]][j - 1] == 0) {
                                        continue;
                                    }

                                    int energy;
                                    // cout << j-i-1 << ":" << type << ":" << i2r[L2_nuc] << ":" << i2r[R2_nuc] <<
                                    // ":" << dummy_str << endl;
                                    // energy = HairpinE(j - i - 1, type,
                                    // i2r[L2_nuc],
                                    // i2r[R2_nuc],
                                    // dummy_str);
                                    energy = E_hairpin(j - i - 1, type, i2r[L2_nuc], i2r[R2_nuc], dummy_str, P);
                                    // cout << "HairpinE(" << j-i-1 << "," << type << "," << i2r[L2_nuc] << "," <<
                                    // i2r[R2_nuc] << ")" << " at " << i << "," << j << ":" << energy << endl; cout
                                    // << i << " " << j  << " " << energy << ":" << i2n[L_nuc] << "-" << i2n[R_nuc]
                                    // << "<-" << i2n[L2_nuc] << "-" << i2n[R2_nuc] << endl;
                                    C[ij][L][R] = MIN2(energy, C[ij][L][R]);
                                }
                            }
                        }

                        // interior loop
                        // cout << i+1 << " " <<  MIN2(j-2-TURN,i+MAXLOOP+1) << endl;
                        for (int p = i + 1; p <= MIN2(j - 2 - TURN, i + MAXLOOP + 1);
                                p++) { // loop for position q, p
                            int minq = j - i + p - MAXLOOP - 2;
                            if (minq < p + 1 + TURN)
                                minq = p + 1 + TURN;
                            for (int q = minq; q < j; q++) {

                                int pq = getIndx(p, q, max_bp_distance_final, indx);

                                for (unsigned int Lp = 0; Lp < pos2nuc[p].size(); Lp++) {
                                    int Lp_nuc = pos2nuc[p][Lp];
                                    if (options_.nucleotide_constraints == 1 && i2r[Lp_nuc] != NucConst[p]) {
                                        continue;
                                    }

                                    if (options_.DEPflg && p == i + 1 && Dep1[ii2r[L_nuc * 10 + Lp_nuc]][i] == 0) {
                                        continue;
                                    }
                                    if (options_.DEPflg && p == i + 2 && Dep2[ii2r[L_nuc * 10 + Lp_nuc]][i] == 0) {
                                        continue;
                                    }

                                    for (unsigned int Rq = 0; Rq < pos2nuc[q].size(); Rq++) { // nucleotide for p, q
                                        int Rq_nuc = pos2nuc[q][Rq];
                                        if (options_.nucleotide_constraints == 1 && i2r[Rq_nuc] != NucConst[q]) {
                                            continue;
                                        }

                                        if (options_.DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc * 10 + R_nuc]][q] == 0) {
                                            continue;
                                        }
                                        if (options_.DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc * 10 + R_nuc]][q] == 0) {
                                            continue;
                                        }

                                        int type_2 = BP_pair[i2r[Lp_nuc]][i2r[Rq_nuc]];

                                        if (type_2 == 0)
                                            continue;
                                        type_2 = rtype[type_2];

                                        //											if
                                        //(noGUclosure) 												if ((type_2 ==
                                        //3)
                                        //														|| (type_2
                                        //== 4))
                                        // if ((p > i + 1)
                                        //															|| (q < j
                                        //- 1))
                                        // continue; /* continue unless stack *//* no_close is removed. It is related
                                        // with BONUS */

                                        //											if(i==8&&j==19){
                                        //												cout << "test:" << p << "-" << q
                                        //<< endl;
                                        //											}

                                        // for each intloops
                                        for (int const L2_nuc : pos2nuc[i + 1]) {
                                            if (options_.nucleotide_constraints == 1 && i2r[L2_nuc] != NucConst[i + 1]) {
                                                continue;
                                            }

                                            if (options_.DEPflg && Dep1[ii2r[L_nuc * 10 + L2_nuc]][i] == 0) {
                                                continue;
                                            }

                                            for (int const R2_nuc : pos2nuc[j - 1]) {
                                                if (options_.nucleotide_constraints == 1 && i2r[R2_nuc] != NucConst[j - 1]) {
                                                    continue;
                                                }

                                                if (options_.DEPflg && Dep1[ii2r[R2_nuc * 10 + R_nuc]][j - 1] == 0) {
                                                    continue;
                                                }

                                                for (int const Lp2_nuc : pos2nuc[p - 1]) {
                                                    if (options_.nucleotide_constraints == 1 && i2r[Lp2_nuc] != NucConst[p - 1]) {
                                                        continue;
                                                    }

                                                    if (options_.DEPflg && Dep1[ii2r[Lp2_nuc * 10 + Lp_nuc]][p - 1] == 0) {
                                                        continue;
                                                    }
                                                    if (p == i + 2 && L2_nuc != Lp2_nuc) {
                                                        continue;
                                                    } // check when a single nucleotide between i and p, this
                                                        // sentence confirm the dependency between Li_nuc and Lp2_nuc
                                                    if (options_.DEPflg && i + 3 == p &&
                                                        Dep1[ii2r[L2_nuc * 10 + Lp2_nuc]][i + 1] == 0) {
                                                        continue;
                                                    } // check dependency between i+1, p-1 (i,X,X,p)

                                                    for (int const Rq2_nuc : pos2nuc[q + 1]) {
                                                        if (q == j - 2 && R2_nuc != Rq2_nuc) {
                                                            continue;
                                                        } // check when a single nucleotide between q and j,this
                                                            // sentence confirm the dependency between Rj_nuc and
                                                            // Rq2_nuc

                                                        if (options_.nucleotide_constraints == 1 && i2r[Rq2_nuc] != NucConst[q + 1]) {
                                                            continue;
                                                        }

                                                        if (options_.DEPflg && Dep1[ii2r[Rq_nuc * 10 + Rq2_nuc]][q] == 0) {
                                                            continue;
                                                        }
                                                        if (options_.DEPflg && q + 3 == j &&
                                                            Dep1[ii2r[Rq2_nuc * 10 + R2_nuc]][q + 1] == 0) {
                                                            continue;
                                                        } // check dependency between q+1, j-1 (q,X,X,j)

                                                        int int_energy = E_intloop(p - i - 1, j - q - 1, type,
                                                                                    type_2, i2r[L2_nuc], i2r[R2_nuc],
                                                                                    i2r[Lp2_nuc], i2r[Rq2_nuc], P);
                                                        // LoopEnergy(p- i- 1,j- q-
                                                        // 1,type,type_2,i2r[L2_nuc],i2r[R2_nuc],i2r[Lp2_nuc],i2r[Rq2_nuc]);

                                                        // int energy =
                                                        //		int_energy
                                                        //		+ C[indx[q]
                                                        //			+ p][Lp][Rq];

                                                        int energy = int_energy + C[pq][Lp][Rq];
                                                        C[ij][L][R] = MIN2(energy, C[ij][L][R]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            } /* end q-loop */
                        }     /* end p-loop */

                        // multi-loop
                        for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                            int Li1_nuc = pos2nuc[i + 1][Li1];
                            if (options_.nucleotide_constraints == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                                continue;
                            }

                            if (options_.DEPflg && Dep1[ii2r[L_nuc * 10 + Li1_nuc]][i] == 0) {
                                continue;
                            }

                            for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                                int Rj1_nuc = pos2nuc[j - 1][Rj1];
                                if (options_.nucleotide_constraints == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                    continue;
                                }

                                if (options_.DEPflg && Dep1[ii2r[Rj1_nuc * 10 + R_nuc]][j - 1] == 0) {
                                    continue;
                                }
                                // if(DEPflg && j-i == 2 && i <= nuclen - 2 && Dep2[ii2r[L_nuc*10+R_nuc]][i] ==
                                // 0){continue;}
                                if (options_.DEPflg && (j - 1) - (i + 1) == 2 &&
                                    Dep2[ii2r[Li1_nuc * 10 + Rj1_nuc]][i + 1] == 0) {
                                    continue;
                                } // 2014/10/8
                                    // i-jが近いときは、MLclosingする必要はないのでは。少なくとも3つのステムが含まれなければならない。それには、５＋５＋２（ヘアピン2個分＋2塩基）の長さが必要。

                                int energy = DMl2
                                    [i + 1][Li1]
                                    [Rj1]; // 長さが2個短いときの、複合マルチループ。i'=i+1を選ぶと、j'=(i+1)+(l-2)-1=i+l-2=j-1(because:j=i+l-1)
                                int tt = rtype[type];

                                energy += P->MLintern[tt];
                                if (tt > 2)
                                    energy += P->TerminalAU;

                                energy += P->MLclosing;
                                // cout << "TEST:" << i << " " << j << " " << energy << endl;
                                C[ij][L][R] = MIN2(energy, C[ij][L][R]);

                                //									if(C[ij][L][R] == -1130 && ij
                                //== 10091){
                                // exit(0);
                                //									}
                            }
                        }

                        //							cout << "ok" << endl;
                    }


                    // fill M
                    // create M[ij] from C[ij]
                    if (type) {
                        int energy_M = C[ij][L][R];
                        if (type > 2)
                            energy_M += P->TerminalAU;

                        energy_M += P->MLintern[type];
                        M[ij][L][R] = energy_M;
                    }

                    // create M[ij] from M[i+1][j]
                    for (unsigned int Li1 = 0; Li1 < pos2nuc[i + 1].size(); Li1++) {
                        int Li1_nuc = pos2nuc[i + 1][Li1];
                        if (options_.nucleotide_constraints == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                            continue;
                        }
                        if (options_.DEPflg && Dep1[ii2r[L_nuc * 10 + Li1_nuc]][i] == 0) {
                            continue;
                        }

                        // int energy_M = M[indx[j]+i+1][Li1][R]+P->MLbase;
                        int energy_M = M[getIndx(i + 1, j, max_bp_distance_final, indx)][Li1][R] + P->MLbase;
                        M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
                    }

                    // create M[ij] from M[i][j-1]
                    for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                        int Rj1_nuc = pos2nuc[j - 1][Rj1];
                        if (options_.nucleotide_constraints == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                            continue;
                        }
                        if (options_.DEPflg && Dep1[ii2r[Rj1_nuc * 10 + R_nuc]][j - 1] == 0) {
                            continue;
                        }

                        // int energy_M = M[indx[j-1]+i][L][Rj1]+P->MLbase;
                        int energy_M = M[getIndx(i, j - 1, max_bp_distance_final, indx)][L][Rj1] + P->MLbase;
                        M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
                    }

                    /* modular decomposition -------------------------------*/
                    for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
                        // cout << k << endl;
                        for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                            int Rk1_nuc = pos2nuc[k - 1][Rk1];
                            if (options_.nucleotide_constraints == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                                continue;
                            }
                            // if(DEPflg && k == i + 2 && Dep1[ii2r[L_nuc*10+Rk1_nuc]][k-1] == 0){ continue;} //
                            // dependency between i and k - 1(=i+1) if(DEPflg && k == i + 3 &&
                            // Dep2[ii2r[L_nuc*10+Rk1_nuc]][k-1] == 0){ continue;} // dependency between i and k -
                            // 1(=i+2)

                            for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                                int Lk_nuc = pos2nuc[k][Lk];
                                if (options_.nucleotide_constraints == 1 && i2r[Lk_nuc] != NucConst[k]) {
                                    continue;
                                }
                                if (options_.DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                                    continue;
                                } // dependency between k - 1 and k
                                // if(DEPflg && (k-1) - i + 1 == 2 && Dep2[ii2r[Rk1_nuc*10+L_nuc]][k-1] == 0){
                                // continue;} // dependency between i and k - 1

                                // cout << i << " " << k-1 << ":" << M[indx[k-1]+i][L][Rk1] << "," << k << " " << j
                                // << ":" << M[indx[j]+k][Lk][R] << endl; int energy_M =
                                // M[indx[k-1]+i][L][Rk1]+M[indx[j]+k][Lk][R];
                                int energy_M = M[getIndx(i, k - 1, max_bp_distance_final, indx)][L][Rk1] +
                                                M[getIndx(k, j, max_bp_distance_final, indx)][Lk][R];
                                DMl[i][L][R] = MIN2(energy_M, DMl[i][L][R]);
                                M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
                            }
                        }
                    }

                    //						if(i == 3 && j == 7)
                    // cout << i << " " << j << ":" << C[ij][L][R] << " " << L << "-" << R << endl;
                    // vwxyがあるので、ここを複数回訪れることがある。
                    //なので、MIN2を取っておく。
                    // if(i2r[L_nuc] == NucConst[i] && i2r[R_nuc] == NucConst[j]){
                    chkC[ij] = MIN2(chkC[ij], C[ij][L][R]);
                    chkM[ij] = MIN2(chkM[ij], M[ij][L][R]);
                    //} このループは多分意味がない。
                }
            }
        }
        // rotate DMl arrays
        vector<array<array<int, 4>, 4>> FF;
        FF = DMl2;
        DMl2 = DMl1;
        DMl1 = DMl;
        DMl = FF;
        for (int j = 1; j <= nuclen; j++) {
            DMl[j].fill({{INF, INF, INF, INF}});
        }
    }

    // Fill F matrix
    // F[1], as well as the rest of F, is default-initialized to 0.

    for (unsigned int L1 = 0; L1 < pos2nuc[1].size(); L1++) {
        int L1_nuc = pos2nuc[1][L1];
        if (options_.nucleotide_constraints == 1 && i2r[L1_nuc] != NucConst[1]) {
            continue;
        }

        for (int j = 2; j <= nuclen; j++) {
            for (unsigned int Rj = 0; Rj < pos2nuc[j].size(); Rj++) {
                int Rj_nuc = pos2nuc[j][Rj];
                if (options_.nucleotide_constraints == 1 && i2r[Rj_nuc] != NucConst[j]) {
                    continue;
                }

                if (options_.DEPflg && j == 2 && Dep1[ii2r[L1_nuc * 10 + Rj_nuc]][1] == 0) {
                    continue;
                }
                if (options_.DEPflg && j == 3 && Dep2[ii2r[L1_nuc * 10 + Rj_nuc]][1] == 0) {
                    continue;
                }

                F[j][L1][Rj] = INF;

                int type_L1Rj = BP_pair[i2r[L1_nuc]][i2r[Rj_nuc]];
                if (type_L1Rj) {
                    //						if(opt_flg_1 && opt_flg_j){
                    int au_penalty = 0;
                    if (type_L1Rj > 2)
                        au_penalty = P->TerminalAU;
                    if (j <= max_bp_distance_final)
                        F[j][L1][Rj] =
                            MIN2(F[j][L1][Rj], C[getIndx(1, j, max_bp_distance_final, indx)][L1][Rj] + au_penalty); // recc 1
                    // F[j][L1][Rj] = MIN2(F[j][L1][Rj], C[indx[j] + 1][L1][Rj] + au_penalty); // recc 1
                    //						}
                }

                // create F[j] from F[j-1]
                for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                    int Rj1_nuc = pos2nuc[j - 1][Rj1];
                    if (options_.nucleotide_constraints == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                        continue;
                    }
                    if (options_.DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                        continue;
                    }

                    F[j][L1][Rj] = MIN2(F[j][L1][Rj], F[j - 1][L1][Rj1]); // recc 2
                }

                // create F[j] from F[k-1] and C[k][j]
                // for (int k = 2; k <= j - TURN - 1; k++) { // Is this correct?
                for (int k = MAX2(2, j - max_bp_distance_final + 1); k <= j - TURN - 1; k++) { // Is this correct?
                    for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                        int Rk1_nuc = pos2nuc[k - 1][Rk1];
                        if (options_.nucleotide_constraints == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                            continue;
                        }
                        if (options_.DEPflg && k == 3 && Dep1[ii2r[L1_nuc * 10 + Rk1_nuc]][1] == 0) {
                            continue;
                        } // dependency between 1(i) and 2(k-1)
                        if (options_.DEPflg && k == 4 && Dep2[ii2r[L1_nuc * 10 + Rk1_nuc]][1] == 0) {
                            continue;
                        } // dependency between 1(i) and 3(k-1)

                        for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                            int Lk_nuc = pos2nuc[k][Lk];
                            if (options_.nucleotide_constraints == 1 && i2r[Lk_nuc] != NucConst[k]) {
                                continue;
                            }

                            if (options_.DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                                continue;
                            } // dependency between k-1 and k

                            int type_LkRj = BP_pair[i2r[Lk_nuc]][i2r[Rj_nuc]];

                            int au_penalty = 0;
                            if (type_LkRj > 2)
                                au_penalty = P->TerminalAU;
                            // int kj = indx[j] + k;
                            int kj = getIndx(k, j, max_bp_distance_final, indx);

                            int energy = F[k - 1][L1][Rk1] + C[kj][Lk][Rj] + au_penalty; // recc 4

                            F[j][L1][Rj] = MIN2(F[j][L1][Rj], energy);
                        }
                    }
                }

                // cout << j << ":" << F[j][L1][Rj] << " " << i2n[L1_nuc] << "-" << i2n[Rj_nuc] << endl;

                // test
                if (j == nuclen) {
                    cout << i2n[L1_nuc] << "-" << i2n[Rj_nuc] << ":" << F[j][L1][Rj] << endl;
                }
            }
        }
    }

    int minL, minR, MFE;
    MFE = INF;
    for (unsigned int L = 0; L < pos2nuc[1].size(); L++) {
        int L_nuc = pos2nuc[1][L];
        if (options_.nucleotide_constraints == 1 && i2r[L_nuc] != NucConst[1]) {
            continue;
        }
        for (unsigned int R = 0; R < pos2nuc[nuclen].size(); R++) {
            int R_nuc = pos2nuc[nuclen][R];
            if (options_.nucleotide_constraints == 1 && i2r[R_nuc] != NucConst[nuclen]) {
                continue;
            }

            if (F[nuclen][L][R] < MFE) {
                MFE = F[nuclen][L][R];
                minL = L;
                minR = R;
            }
        }
    }

    if (MFE == INF) {
        printf("Mininum free energy is not defined.\n");
        exit(1);
    }

    if (options_.random_backtrack) {
        // Fill F2 matrix
        for (int l = 5; l <= nuclen; l++) {
            if (l > max_bp_distance_final)
                break;
            cout << "process F2:" << l << endl;

            for (int i = 1; i <= nuclen - l + 1; i++) {
                int j = i + l - 1;

                for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                    int L_nuc = pos2nuc[i][L];

                    for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {
                        int R_nuc = pos2nuc[j][R];
                        int ij = getIndx(i, j, max_bp_distance_final, indx);

                        F2[ij][L][R] = 0;

                        int type = BP_pair[i2r[L_nuc]][i2r[R_nuc]];

                        // from i, j-1 -> i, j
                        for (unsigned int R1 = 0; R1 < pos2nuc[j - 1].size(); R1++) {
                            int R1_nuc = pos2nuc[j - 1][R1];
                            if (options_.DEPflg && Dep1[ii2r[R1_nuc * 10 + R_nuc]][j - 1] == 0) {
                                continue;
                            }
                            int ij1 = getIndx(i, j - 1, max_bp_distance_final, indx);
                            F2[ij][L][R] = MIN2(F2[ij][L][R], F2[ij1][L][R1]);
                        }
                        // from i-1, j -> i, j
                        for (unsigned int L1 = 0; L1 < pos2nuc[i + 1].size(); L1++) {
                            int L1_nuc = pos2nuc[i + 1][L1];
                            if (options_.DEPflg && Dep1[ii2r[L_nuc * 10 + L1_nuc]][i] == 0) {
                                continue;
                            }
                            int i1j = getIndx(i + 1, j, max_bp_distance_final, indx);
                            F2[ij][L][R] = MIN2(F2[ij][L][R], F2[i1j][L1][R]);
                        }

                        // from C
                        int au_penalty = 0;
                        if (type > 2)
                            au_penalty = P->TerminalAU;
                        if (j - i + 1 <= max_bp_distance_final) {
                            F2[ij][L][R] = MIN2(F2[ij][L][R], C[ij][L][R] + au_penalty);
                            // cout << "test:" << F2[ij][L][R] << endl;
                        }

                        // Bifucation
                        /* modular decomposition -------------------------------*/
                        for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
                            // cout << k << endl;
                            //			    				if((k - 1) - i + 1 > w_tmp ||
                            //			    					j - k + 1 > w_tmp)
                            //			    						continue;

                            for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                                int Rk1_nuc = pos2nuc[k - 1][Rk1];

                                for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                                    int Lk_nuc = pos2nuc[k][Lk];
                                    if (options_.DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                                        continue;
                                    } // dependency between k - 1 and k

                                    int energy = F2[getIndx(i, k - 1, max_bp_distance_final, indx)][L][Rk1] +
                                                    F2[getIndx(k, j, max_bp_distance_final, indx)][Lk][R];
                                    F2[ij][L][R] = MIN2(F2[ij][L][R], energy);
                                }
                            }
                        }

                        //							cout << i << "," << j << "," << L << "," << R << "," << F2[ij][L][R]
                        //<< endl;
                    }
                }
            }
        }
    }

    //		string optseq;
    //		optseq.resize(nuclen+1, 'N');
    //		optseq[0] = ' ';
    //		backtrackR(&optseq, &*sector, &*base_pair, C, M, F,
    //					indx, minL, minR, P, NucConst, pos2nuc, NCflg, i2r, nuclen, w_tmp, BP_pair,
    //i2n, rtype, ii2r, Dep1, Dep2, DEPflg, predefHPN, predefHPN_E, substr, n2i, NucDef);

    array<stack, 500> sector;
    string optseq(nuclen + 1, 'N');
    optseq[0] = ' ';

    string optseq_org(nuclen + 1, 'N');
    optseq_org[0] = ' ';

    if (options_.random_backtrack) {
        backtrack2(&optseq, sector, base_pair, C, M, F2, indx, minL, minR, P, NucConst, pos2nuc, options_.nucleotide_constraints, i2r,
                    nuclen, max_bp_distance_final, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, options_.DEPflg, predefHPN_E, substr,
                    n2i, NucDef);
    } else {
        backtrack(&optseq, sector, base_pair, C, M, F, indx, minL, minR, P, NucConst, pos2nuc, options_.nucleotide_constraints, i2r,
                    nuclen, max_bp_distance_final, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, options_.DEPflg, predefHPN_E, substr, n2i,
                    NucDef);
    }

    //塩基Nの修正
    for (int i = 1; i <= nuclen; i++) {
        if (optseq[i] == 'N') {
            for (unsigned int R = 0; R < pos2nuc[i].size();
                    R++) { // check denendency with the previous and next nucleotide
                int R_nuc = pos2nuc[i][R];
                if (options_.nucleotide_constraints == 1 && i2r[R_nuc] != NucConst[i]) {
                    continue;
                }

                if (i != 1 && optseq[i - 1] != 'N') { // check consistensy with the previous nucleotide
                    int R_prev_nuc = n2i.at(optseq[i - 1]);
                    if (options_.DEPflg && Dep1[ii2r[R_prev_nuc * 10 + R_nuc]][i - 1] == 0) {
                        continue;
                    }
                }
                if (i != nuclen && optseq[i + 1] != 'N') { // check consistensy with the next nucleotide
                    int R_next_nuc = n2i.at(optseq[i + 1]);
                    if (options_.DEPflg && Dep1[ii2r[R_nuc * 10 + R_next_nuc]][i] == 0) {
                        continue;
                    }
                }
                if (i < nuclen - 1 && optseq[i + 2] != 'N') { // check consistensy with the next nucleotide
                    int R_next_nuc = n2i.at(optseq[i + 2]);
                    if (options_.DEPflg && Dep2[ii2r[R_nuc * 10 + R_next_nuc]][i] == 0) {
                        continue;
                    }
                }

                optseq[i] = i2n[R_nuc];
                break;
            }
            // cout << i << ":" << optseq[i] << endl;
        }
    }
    //塩基V,W,X,Yの修正
    for (int i = 1; i <= nuclen; i++) {

        if (optseq[i] == 'V' || optseq[i] == 'W') {
            optseq[i] = 'U';
        } else if (optseq[i] == 'X' || optseq[i] == 'Y') {
            optseq[i] = 'G';
        }
    }

    // 2次構造情報の表示
    string optstr;
    optstr.resize(nuclen + 1, '.');
    optstr[0] = ' ';
    for (int i = 1; i <= base_pair[0].i; i++) {
        optstr[base_pair[i].i] = '(';
        optstr[base_pair[i].j] = ')';
    }

    // show original amino acids
    for (int i = 0; i < aalen; i++) {
        cout << aaseq[i] << "  ";
    }
    cout << endl;
    // check amino acids of desinged DNA
    int j = 0;
    for (unsigned int i = 1; i < optseq.size(); i = i + 3) {
        char aa = codon_table.c2a(n2i.at(optseq[i]), n2i.at(optseq[i + 1]), n2i.at(optseq[i + 2]));
        cout << aa << "  ";
        if (aaseq[j] != aa) {
            cerr << j + 1 << "-th amino acid differs:" << aaseq[j] << ":" << aa << endl;
        }
        j++;
    }
    cout << endl;
    // Display of optimal bases and secondary structure
    string optseq_disp = optseq;
    optseq_disp.erase(0, 1);
    optstr.erase(0, 1);
    cout << optseq_disp << endl;
    cout << optstr << endl;
    cout << "MFE:" << float(MFE) / 100 << " kcal/mol" << endl;

    if (options_.partial_opt == true) {
        // Creation of partial amino acid sequence
        for (int I = 0; I < n_inter; I++) {
            int aa_fm = (ofm[I] - 1) / 3 + 1; // 1-based
            int aa_to = oto[I] / 3;           // 1-based
            int part_aalen = aa_to - aa_fm + 1;

            vector<char> part_aaseq;
            part_aaseq.reserve(part_aalen);
            part_aaseq[part_aalen] = '\0';
            int j = 0;
            for (int i = aa_fm; i <= aa_to; i++) {
                part_aaseq[j++] = aaseq[i - 1]; // convert to 0-based
            }

            cout << aa_fm << ":" << aa_to << endl;
            cout << part_aalen << endl;
            cout << &part_aaseq << endl;

            string part_optseq = rev_fold_step1(&part_aaseq[0], part_aalen, codon_table, options_.codons_excluded);
            rev_fold_step2(&part_optseq, &part_aaseq[0], part_aalen, codon_table, options_.codons_excluded);
            // combine optseq_rev and optseq
            j = 1;
            for (int i = ofm[I]; i <= oto[I]; i++) {
                optseq[i] = part_optseq[j++];
            }
        }
        fixed_fold(optseq, indx, max_bp_distance_final, predefHPN_E, BP_pair, P, &aaseq[0], codon_table);
        // fixed_fold(optseq, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
    }

    if (options_.show_memory_use) {
        // get process ID

        pid_t pid = getpid();
        stringstream ss;
        ss << "/proc/" << pid << "/status";
        int m = getMemoryUsage(ss.str());
        if (m == -1) {
            cerr << "Cannot get memory usage information from " << ss.str() << endl;
        } else {
            cout << "Memory(VmRSS): " << float(m) / 1024 << " Mb" << endl;
        }
    }

    free(P);
}


void fixed_fold(string optseq, vector<int> const & indx, const int &w, map<string, int> &predefE, const int (&BP_pair)[5][5],
                paramT *P, char const *aaseq, const codon& codon_table) {
    int nuclen = optseq.size() - 1;
    int aalen = (optseq.size() - 1) / 3;
    int size = getMatrixSize(nuclen, w);
    vector<int> C;
	C.reserve(size);
    vector<int> M;
	M.reserve(size);
    vector<int> F;
	F.reserve(nuclen + 1);
    vector<int> DMl;
	DMl.reserve(nuclen + 1);
    vector<int> DMl1;
	DMl1.reserve(nuclen + 1);
    vector<int> DMl2;
	DMl2.reserve(nuclen + 1);
    vector<bond> base_pair;
	base_pair.reserve(nuclen / 2);

    map<char, int> n2i = make_n2i();
    vector<int> ioptseq;
	ioptseq.reserve(nuclen + 1);
    ioptseq[0] = 0;
    for (int i = 1; i <= nuclen; i++) {
        ioptseq[i] = n2i[optseq[i]];
        //		cout << optseq[i] << endl;
    }
    //	exit(0);

    fixed_init_matrix(nuclen, size, C, M, &F[0], &DMl[0], &DMl1[0], &DMl2[0]);
    int rtype[7] = {0, 2, 1, 4, 3, 6, 5};

    // investigate mfe
    const char dummy_str[10] = "XXXXXXXXX";
    for (int l = 5; l <= nuclen; l++) {
        if (l > w)
            break;
        // cout << "process:" << l << endl;

        //	  for(int l = 5; l <= 5; l++){
        for (int i = 1; i <= nuclen - l + 1; i++) {
            int j = i + l - 1;
            int ij = getIndx(i, j, w, indx);
            C[ij] = INF;
            M[ij] = INF;
            //			cout << "test:" << j << endl;
            int type = BP_pair[ioptseq[i]][ioptseq[j]];

            if (type) {
                // hairpin
                int energy = E_hairpin(j - i - 1, type, ioptseq[i + 1], ioptseq[j - 1], dummy_str, P);
                C[ij] = MIN2(energy, C[ij]);

                if (l == 5 || l == 6 || l == 8) {
                    string hpn = "";
                    for (int k = i; k <= j; k++) {
                        hpn += optseq[k];
                    }
                    // cout << i << ":" << j << "=" << hpn << endl;
                    if (predefE.count(hpn) > 0) {
                        C[ij] = predefE[hpn];
                    }
                }

                // interior loop
                for (int p = i + 1; p <= MIN2(j - 2 - TURN, i + MAXLOOP + 1); p++) { // loop for position q, p
                    int minq = j - i + p - MAXLOOP - 2;
                    if (minq < p + 1 + TURN) {
                        minq = p + 1 + TURN;
                    }
                    for (int q = minq; q < j; q++) {

                        int pq = getIndx(p, q, w, indx);

                        int type_2 = BP_pair[ioptseq[p]][ioptseq[q]];

                        if (type_2 == 0) {
                            continue;
                        }
                        type_2 = rtype[type_2];

                        int int_energy = E_intloop(p - i - 1, j - q - 1, type, type_2, ioptseq[i + 1], ioptseq[j - 1],
                                                   ioptseq[p - 1], ioptseq[q + 1], P);

                        int energy = int_energy + C[pq];
                        C[ij] = MIN2(energy, C[ij]);

                    } /* end q-loop */
                }     /* end p-loop */
                // cout << i << "," << j << ":" << C[ij] << endl;

                // multi-loop
                energy = DMl2[i + 1];
                int tt = rtype[type];

                energy += P->MLintern[tt];
                if (tt > 2)
                    energy += P->TerminalAU;

                energy += P->MLclosing;
                C[ij] = MIN2(energy, C[ij]);
            } else
                C[ij] = INF;

            // fill M
            // create M[ij] from C[ij]
            if (type) {
                int energy_M = C[ij];
                if (type > 2) {
                    energy_M += P->TerminalAU;
                }

                energy_M += P->MLintern[type];
                M[ij] = energy_M;
            }

            // create M[ij] from M[i+1][j]
            int energy_M = M[getIndx(i + 1, j, w, indx)] + P->MLbase;
            M[ij] = MIN2(energy_M, M[ij]);

            // create M[ij] from M[i][j-1]
            energy_M = M[getIndx(i, j - 1, w, indx)] + P->MLbase;
            M[ij] = MIN2(energy_M, M[ij]);

            /* modular decomposition -------------------------------*/
            for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
                int energy_M = M[getIndx(i, k - 1, w, indx)] + M[getIndx(k, j, w, indx)];
                DMl[i] = MIN2(energy_M, DMl[i]);
                M[ij] = MIN2(energy_M, M[ij]);
            }
        }
        // rotate DMl arrays
        vector<int> FF;
        FF.reserve(nuclen + 1);
        for (int j = 1; j <= nuclen; j++) {
            FF[j] = DMl2[j];
            DMl2[j] = DMl1[j];
            DMl1[j] = DMl[j];
            DMl[j] = FF[j];
        }
        for (int j = 1; j <= nuclen; j++) {
            DMl[j] = INF;
        }
    }

    // Fill F matrix
    // Initialize F[1]
    F[1] = 0;
    for (int j = 2; j <= nuclen; j++) {
        F[j] = INF;
        int type = BP_pair[ioptseq[1]][ioptseq[j]];
        if (type) {
            int au_penalty = 0;
            if (type > 2)
                au_penalty = P->TerminalAU;
            if (j <= w)
                F[j] = MIN2(F[j], C[getIndx(1, j, w, indx)] + au_penalty); // recc 1
        }

        // create F[j] from F[j-1]
        F[j] = MIN2(F[j], F[j - 1]); // recc 2

        for (int k = MAX2(2, j - w + 1); k <= j - TURN - 1; k++) { // Is this correct?
            int type_k = BP_pair[ioptseq[k]][ioptseq[j]];

            int au_penalty = 0;
            if (type_k > 2)
                au_penalty = P->TerminalAU;
            int kj = getIndx(k, j, w, indx);

            int energy = F[k - 1] + C[kj] + au_penalty; // recc 4
            F[j] = MIN2(F[j], energy);
        }
    }

    int MFE = F[nuclen];
    //	showFixedMatrix(C, indx, nuclen, w);
    //	showFixedMatrix(M, indx, nuclen, w);
    //	for(int i = 1; i <= nuclen; i++){
    //		cout << "FMAT: " << i << " " << F[i] << "  " <<  w << endl;
    //	}
	// AMW TODO
    fixed_backtrack(optseq, &base_pair[0], C, M, &F[0], indx, P, nuclen, w, BP_pair, predefE);
    //}

    //	cout << "end" << endl;
    //	exit(0);
    // 2次構造情報の表示
    string optstr;
    optstr.resize(nuclen + 1, '.');
    optstr[0] = ' ';
    for (int i = 1; i <= base_pair[0].i; i++) {
        optstr[base_pair[i].i] = '(';
        optstr[base_pair[i].j] = ')';
    }

    // show original amino acids
    for (int i = 0; i < aalen; i++) {
        cout << aaseq[i] << "  ";
    }
    cout << endl;
    // check amino acids of desinged DNA
    int j = 0;
    for (unsigned int i = 1; i < optseq.size(); i = i + 3) {
        char aa = codon_table.c2a(n2i[optseq[i]], n2i[optseq[i + 1]], n2i[optseq[i + 2]]);
        cout << aa << "  ";
        if (aaseq[j] != aa) {
            cerr << j + 1 << "-th amino acid differs:" << aaseq[j] << ":" << aa << endl;
        }
        j++;
    }
    cout << endl;
    // Display of optimal bases and secondary structure
    optseq.erase(0, 1);
    optstr.erase(0, 1);
    cout << optseq << endl;
    cout << optstr << endl;
    cout << "MFE:" << float(MFE) / 100 << " kcal/mol" << endl;
}

// Get all the possibilities at each position. This is essentially the process
// by which we transform the AA sequence to a nucleotide [possibility] sequence.
vector<vector<int>> getPossibleNucleotide(string const & aaseq, codon &codon_table, map<char, int> const & n2i,
                                          string const & excludedCodons) {

    vector<vector<int>> v;

    int nuclen = aaseq.size() * 3;
    v.resize(nuclen + 1);

    for (unsigned int i = 0; i < aaseq.size(); i++) {
        int nucpos = i * 3 + 1;
        char aa = aaseq[i];
        // cout << "test: " << aa << endl;
        // AMW TODO: should get codons from the codon table, filtering out
        // exclusions, *once* rather than generating new vector<string>s for each
        // aa position. This function is only called once, but it's still silly.
        vector<string> codons = codon_table.getCodons(aa, excludedCodons);
        for (int k = 0; k < 3; k++) { // for each codon position
            nucpos += k;

            if (aa == 'L' && k == 1) {
                v[nucpos].push_back(n2i.at('V'));
                v[nucpos].push_back(n2i.at('W'));
            } else if (aa == 'R' && k == 1) {
                v[nucpos].push_back(n2i.at('X'));
                v[nucpos].push_back(n2i.at('Y'));
            } else {
                bool flg_A, flg_C, flg_G, flg_U;
                flg_A = flg_C = flg_G = flg_U = false;

                for (unsigned int j = 0; j < codons.size(); j++) { // for each codon corresp. to each aa
                    char nuc = codons[j][k];
                    if (nuc == 'A' && flg_A == 0) {
                        v[nucpos].push_back(n2i.at(nuc));
                        flg_A = 1;
                    } else if (nuc == 'C' && flg_C == 0) {
                        v[nucpos].push_back(n2i.at(nuc));
                        flg_C = 1;
                    } else if (nuc == 'G' && flg_G == 0) {
                        v[nucpos].push_back(n2i.at(nuc));
                        flg_G = 1;
                    } else if (nuc == 'U' && flg_U == 0) {
                        v[nucpos].push_back(n2i.at(nuc));
                        flg_U = 1;
                    }
                }
            }
            nucpos -= k;
        }
    }
    return v;
}

void allocate_arrays(
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

void allocate_F2(int len, vector<int> const & indx, int w, vector<vector<int>> &pos2nuc, vector<vector<vector<int>>> & f2) {
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

inline auto getMemoryUsage(const string &fname) -> int {
    // cout << fname << endl;
    ifstream ifs(fname.c_str());
    if (ifs) {
        string line;
        while (getline(ifs, line)) { // $B9T$NFI$_9~$_(B
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