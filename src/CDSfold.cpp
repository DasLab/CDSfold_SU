/*
 * CDSfold.cpp

 *
 *  Created on: Sep 2, 2014
 *      Author: terai
 */
#define MIN2(A, B) ((A) < (B) ? (A) : (B))
#define MAX2(A, B) ((A) > (B) ? (A) : (B))
#define TURN 3
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <unistd.h>

extern "C" {
#include <cctype>
#include "fold.h"
#include "fold_vars.h"
#include <climits>
#include <cmath>
#include "params.h"
#include "part_func.h"
#include <cstdio>
#include <cstdlib>
#include "utils.h"
}

#include "AASeqConverter.hpp"
#include "CDSfold.hpp"
#include "CDSfold_rev.hpp"
#include "codon.hpp"
#include "fasta.hpp"
//#include <algorithm>
//#include <sys/time.h>
//#include <sys/resource.h>

int BP_pair[5][5] =
    /* _  A  C  G  U  */
    {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 5}, {0, 0, 0, 1, 0}, {0, 0, 2, 0, 3}, {0, 6, 0, 4, 0}};
int rtype[7] = {0, 2, 1, 4, 3, 6, 5};

//#define MAXLOOP 20
#define noGUclosure 0
int test;

using namespace std;
auto main(int argc, char *argv[]) -> int {
    // printf("%d\n%ld", INT_MAX, LONG_MAX);
    int W = 0;            // -w
    string exc = "";      // -e
    int m_disp = 0;       // -M
    int m_estim = 0;      // -U
    int rand_tb_flg = 0;  //-R
    int rev_flg = 0;      // -r
    int part_opt_flg = 0; // -f and -t
    int opt_fm = 0;       // -f
    int opt_to = 0;       // -t
    // get options
    {
        int opt;
        while ((opt = getopt(argc, argv, "w:e:f:t:rMUR")) != -1) {
            switch (opt) {
            case 'w':
                W = atoi(optarg);
                break;
            case 'e':
                exc = string(optarg);
                break;
            case 'M':
                m_disp = 1;
                break;
            case 'U':
                m_estim = 1;
                break;
            case 'R':
                rand_tb_flg = 1;
                break;
            case 'r':
                rev_flg = 1;
                break;
            case 'f':
                opt_fm = atoi(optarg);
                part_opt_flg = 1;
                break;
            case 't':
                opt_to = atoi(optarg);
                part_opt_flg = 1;
                break;
            }
        }
    }
    // exit(0);

    // -R オプションに関するチェック
    if (rand_tb_flg) {
        if (W != 0 || exc != "" || m_disp || m_estim || rev_flg || part_opt_flg) {
            cerr << "The -R option must not be used together with other options." << endl;
            return 0;
        }
    }

    if (part_opt_flg) {
        // 部分最適化が指定された。
        if (opt_fm == 0 || opt_to == 0) {
            cerr << "The -f and -t option must be used together." << endl;
            exit(1);
        }
        if (opt_fm < 1) {
            cerr << "The -f value must be 1 or more." << endl;
            exit(1);
        }
        if (opt_to < 1) {
            cerr << "The -t value must be 1 or more." << endl;
            exit(1);
        }
        if (opt_to < opt_fm) {
            cerr << "The -f value must be smaller than -t value." << endl;
            exit(1);
        }
    }

    // int energy = E(5, 1, 1, 2, "ATGCATGC");
    // int energy = E_Hairpin(5, 1, 1, 2, "ATGCATGC");

    map<char, int> n2i = make_n2i();
    std::array<char, 20> i2n;
    make_i2n(i2n);

    array<int, 20> i2r;
    array<int, 100> ii2r;
    make_i2r(i2r);
    make_ii2r(ii2r);

    AASeqConverter conv;

    codon codon_table;
    // codon_table.Table();

    const char dummy_str[10] = "XXXXXXXXX";
    // int NCflg = 1;
    //	int TB_CHK_flg = 0;
    //	int preHPN_flg = 1;
    int TEST = 1;
    int DEPflg = 1;
    int NCflg = 0;
    //	char *NucDef =
    //"*AUGGGUCUUCCAGUGUCAUUACGAGCUGACACCAUUCGAGAUUUAUUACUUGGUGUCAGCUCGAUAAUGACCUGGAAGACCCUUGCUCUUGUGUUAGCUGUGAUCAAUCUCAAGAAUCUGCCACUAGUGUGGCACCCGGGGGAUCCUCAUUUCCCCCGGGGGAAGGCGCUGGUGACGCAUACGGGCAAACCCACUCAUCCGGUGUUUGUCCCGUAUGCGAUCACCAGUCGCACUCCGAUUCUUGAGACUGAUUACAACUUUCACAAGAGCAAUUCCACGUAUUUUAGCGAUUUGGAUAUU";
    const char NucDef[] =
        "*AUGGAGGGGAUUGUCACGGGAGAUCGGCUUGCUUGCGUGGCGCUUCAUGGAAGCUCUUUGCUCCAUGAAGCGUCCGUAAGCAAGUAUACCGAUAUCCCGGGCAUUCUCC"
        "UCCAAUACAUCGAUGAAUUUCCCCUCACUGAUAUUGCCGCGCACGCGCCACGCGAGGCGUGGCAAAGCCUGUGCGAACAGGCGAUCUGUAUCGUCCAUCAUAUUAGCGAC"
        "CGGGGCAUCCUCAAUGAGGAUGUUAAAACCCGGUCGCUGACGAUACAGAUCAACAGUGAGGGGAUGUUCAAGAUGUUUAUG";
    // string tmp_def =
    // "*AUGGCCCCCAUACAGCAGAAGGCACUAAUCAACUGCGAUAUGGGGGAAGCUUACGGGAACUGGGCCUGCGGCCCAGAUCUCGAGCUCCUCCCCAUGAUCGACAUCGCCAACGUGGCGUGUGGAUUUCAUGGGGGGGAUCCAUUAAUAAUGAUGGAAACGGUGCGCAACUGUAAAGCGCACAAUGUGCGCAUAGGGGCGCACCCUGGCCUCCCGGACCUGCAGGGGUUCGGGAGGCGGGAGAUGAAACUCUCCCCUGAAGAGCUCACCGCCAUGACUAUUUAUCAGGUGGGAGCUCUUCAG";
    // char *NucDef = "*AUGAGUCUGGCGUGCAUGGCCAAGUAG";
    // char *NucDef = "*AUGUCUCUCGCGUGCAUGGCCAAGUGA";
    // char *NucDef = "*AUGUCUUUAGCCUGUAUGGCUAAAUAA";
    // cout << "optind is " << optind << endl;
    // const char *NucDef = tmp_def.c_str();
    fasta all_aaseq(argv[optind]); // get all sequences

    cout << "W = " << W << endl;
    cout << "e = " << exc << endl;
    do {
        char *aaseq = all_aaseq.getSeq();
        int aalen = all_aaseq.getSeqLen();

        if (aalen <= 2) {
            cerr << "The amino acid sequence is too short.\n";
            exit(1);
        }

        int n_inter = 0; //今の実装では、n_inter=1 or 2となる。
        int ofm[100];
        int oto[100];
        if (part_opt_flg) {
            // 部分最適化が指定された。
            if (opt_fm == 0 || opt_to == 0) {
                cerr << "The -f and -t option must be used together." << endl;
                exit(1);
            }
            if (opt_fm < 1) {
                cerr << "The -f value must be 1 or more." << endl;
                exit(1);
            }
            if (opt_to < 1) {
                cerr << "The -t value must be 1 or more." << endl;
                exit(1);
            }
            if (opt_to < opt_fm) {
                cerr << "The -f value must be smaller than -t value." << endl;
                exit(1);
            }
            if (opt_to > aalen) {
                opt_to = aalen;
            }

            // Creating partial inverse optimization information
            if (rev_flg) { // Structural removal of the specified area
                ofm[0] = (opt_fm - 1) * 3 + 1;
                oto[0] = opt_to * 3;
                n_inter = 1;
            } else { // Structural stabilization of the specified area
                int l = 0;
                if (opt_fm != 1) {
                    ofm[l] = 1;
                    oto[l++] = (opt_fm - 1) * 3;
                }
                if (opt_to != aalen) {
                    ofm[l] = opt_to * 3 + 1;
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

        int nuclen = aalen * 3;
        int w_tmp;

        if (W == 0) {
            w_tmp = nuclen;
        } else if (W < 10) {
            cerr << "W must be more than 10"
                 << "(you used " << W << ")" << endl;
            exit(1);
        } else if (W > nuclen) {
            w_tmp = nuclen;
        } else {
            w_tmp = W;
        }

        //		w_tmp = 50;// test!
        //		vector<vector<vector<string> > >  substr = conv.getBases(string(aaseq),8, exc);
        vector<vector<vector<string>>> substr = conv.getOriginalBases(string(aaseq), exc);
        vector<vector<int>> Dep1;
        vector<vector<int>> Dep2;
        float ptotal_Mb_base = 0;

        if (m_estim) {
            ptotal_Mb_base = 2 + nuclen * 0.006956;
        } else {
            Dep1 = conv.countNeighborTwoBase(string(aaseq), exc);
            Dep2 = conv.countEveryOtherTwoBase(string(aaseq), exc);
        }

        //		cout << ptotal_Mb_base << endl;

        //		pid_t pid2 = getpid();
        //		stringstream ss2;
        //		ss2 << "/proc/" << pid2 << "/status";
        //		int m2 = getMemoryUsage(ss2.str());
        //		cout << "Memory(VmRSS): "  << float(m2)/1024 << " Mb" << endl;
        // cout << "Estimate: "  << ptotal_Mb_base << " Mb" << endl;
        // exit(0);

        // map<string, int> predefHPN_E;
        map<string, int> predefHPN_E = conv.getBaseEnergy();
        //		vector<vector<vector<vector<pair<int, string> > > > > predefHPN =
        // conv.calcQueryOriginalBaseEnergy(string(aaseq), "");
        vector<vector<vector<vector<pair<int, string>>>>> predefHPN;
        // vector<vector<vector<vector<pair<int, string> > > > > predefHPN =
        // conv.calcQueryOriginalBaseEnergy(string(aaseq), "");
        /*
                        for(unsigned int i = 1; i < predefHPN.size(); i++){
                                cout << ">>position " << i << endl;
                                for(unsigned int l = 0; l < predefHPN[i].size(); l++){
                                        for(unsigned int li = 0; li < predefHPN[i][l].size(); li++){
                                                for(unsigned int rj = 0; rj < predefHPN[i][l][li].size(); rj++){
                                                        if(predefHPN[i][l][li][rj].first != INF)
                                                                continue;
                                                        cout << predefHPN[i][l][li][rj].second << ":" <<
           predefHPN[i][l][li][rj].first << endl;
                                                }
                                        }
                                }
                        }
                        cout << INF << endl;
                        exit(0);
        */

        stack sector[500];

        // vector<int> NucConst = createNucConstraint(NucDef, nuclen, n2i);

        vector<int> NucConst;
        if (NCflg) {
            NucConst = createNucConstraint(NucDef, nuclen, n2i);
        }
        //		for(int i = 1; i <= nuclen; i++)
        //		printf("%d %d\n", i, NucConst[i]);

        // createNucConstraint

        cout << aaseq << endl;
        //		cout << aalen << endl;

        vector<vector<int>> pos2nuc = getPossibleNucleotide(aaseq, aalen, codon_table, n2i, exc);
        //		vector<vector<int> > pos2nuc = getPossibleNucleotide(aaseq, aalen, codon_table, n2i, 'R');
        //		showPos2Nuc(pos2nuc, i2n);
        //		exit(0);
        vector<int> indx(nuclen+1, 0);

        set_ij_indx(indx, nuclen, w_tmp);
        // set_ij_indx(indx, nuclen);

        string optseq;
        optseq.resize(nuclen + 1, 'N');
        optseq[0] = ' ';

        string optseq_org;
        optseq_org.resize(nuclen + 1, 'N');
        optseq_org[0] = ' ';

        //	  int ***C, ***Mbl, ***Mbr, ***Mbb, ***M, ***F, ***Fbr, ***tFbr;
        vector<vector<vector<int>>> C, M, F, F2;
        vector<array<array<int, 4>, 4>> DMl, DMl1, DMl2;
        vector<int> chkC, chkM;
        vector<bond> base_pair;

        //		int n_inter = 1;

        if (m_estim) {
            //		allocate_arrays(nuclen, indx, pos2nuc, pos2nuc, &C, &M, &F);
            float ptotal_Mb_alloc = predict_memory(nuclen, w_tmp, pos2nuc);
            float ptotal_Mb = ptotal_Mb_alloc + ptotal_Mb_base;
            cout << "Estimated memory usage: " << ptotal_Mb << " Mb" << endl;
            return 0;
        }

        paramT *P = scale_parameters();
        update_fold_params();

        //		rev_flg = 0;
        //		if(rev_flg && num_interval == 0){
        if (rev_flg && !part_opt_flg) {
            // reverse mode
            string optseq_rev = rev_fold_step1(aaseq, aalen, codon_table, exc);
            //			rev_fold_step2(&optseq_rev, aaseq, aalen, codon_table, exc, ofm, oto, 1);
            rev_fold_step2(&optseq_rev, aaseq, aalen, codon_table, exc);
            fixed_fold(optseq_rev, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
            free(P);
            break; // returnすると、実行時間が表示されなくなるためbreakすること。
        }

        //		allocate_arrays(nuclen, indx, pos2nuc, pos2nuc, &C, &M, &F);
        allocate_arrays(nuclen, indx, w_tmp, pos2nuc, C, M, F, DMl, DMl1, DMl2, chkC, chkM, base_pair);
        if (rand_tb_flg) {
            allocate_F2(nuclen, indx, w_tmp, pos2nuc, F2);
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
                int ij = getIndx(i, j, w_tmp, indx);

                chkC[ij] = INF;
                chkM[ij] = INF;

                for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                    int L_nuc = pos2nuc[i][L];
                    if (NCflg == 1 && i2r[L_nuc] != NucConst[i]) {
                        continue;
                    }
                    for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {
                        int R_nuc = pos2nuc[j][R];
                        if (NCflg == 1 && i2r[R_nuc] != NucConst[j]) {
                            continue;
                        }
                        // L-R pair must be filtered
                        //						if(j-1==1){
                        //							cout << i << ":" << j << " " << L_nuc << "-"
                        //<< R_nuc << " " << Dep1[ii2r[L_nuc*10+R_nuc]][i] <<endl;
                        //						}
                        if (DEPflg && j - i == 1 && i <= nuclen - 1 && Dep1[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
                            continue;
                        } // nuclen - 1はいらないのでは？
                        if (DEPflg && j - i == 2 && i <= nuclen - 2 && Dep2[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
                            continue;
                        }

                        C[ij][L][R] = INF;
                        M[ij][L][R] = INF;
                        if (rand_tb_flg)
                            F2[ij][L][R] = 0;
                    }
                }
            }
        }

        //		cout << "TEST" << M[13][0][0] << endl;
        // main routine
        for (int l = 5; l <= nuclen; l++) {
            if (l > w_tmp)
                break;
            cout << "process:" << l << endl;

            //	  for(int l = 5; l <= 5; l++){
            for (int i = 1; i <= nuclen - l + 1; i++) {
                int j = i + l - 1;

                int opt_flg_ij = 1;
                if (part_opt_flg) {
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
                    if (NCflg == 1 && i2r[L_nuc] != NucConst[i]) {
                        continue;
                    }
                    //					cout << "ok" << endl;

                    for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {

                        int R_nuc = pos2nuc[j][R];

                        if (NCflg == 1 && i2r[R_nuc] != NucConst[j]) {
                            continue;
                        }

                        // int ij = indx[j] + i;
                        int ij = getIndx(i, j, w_tmp, indx);

                        C[ij][L][R] = INF;
                        M[ij][L][R] = INF;
                        //						cout << i << " " << j << ":" << M[ij][L][R] <<
                        // endl;

                        int type = BP_pair[i2r[L_nuc]][i2r[R_nuc]];

                        if (!type || !opt_flg_ij) {
                            C[ij][L][R] = INF;
                        } else {
                            // hairpin
                            if ((l == 5 || l == 6 || l == 8) && TEST) {
                                for (unsigned int s = 0; s < substr[i][l].size(); s++) {
                                    string hpn = substr[i][l][s];
                                    int hL_nuc = n2i[hpn[0]];
                                    int hL2_nuc = n2i[hpn[1]];
                                    int hR2_nuc = n2i[hpn[l - 2]];
                                    int hR_nuc = n2i[hpn[l - 1]];
                                    if (hL_nuc != i2r[L_nuc])
                                        continue;
                                    if (hR_nuc != i2r[R_nuc])
                                        continue;

                                    if (NCflg == 1) {
                                        string s1 = string(NucDef).substr(i, l);
                                        if (hpn != s1)
                                            continue;
                                    }

                                    //									cout << hpn <<
                                    // endl;
                                    //
                                    if (DEPflg && L_nuc > 4 && Dep1[ii2r[L_nuc * 10 + hL2_nuc]][i] == 0) {
                                        continue;
                                    } // Dependencyをチェックした上でsubstringを求めているので,
                                      // hpnの内部についてはチェックする必要はない。
                                    if (DEPflg && R_nuc > 4 && Dep1[ii2r[hR2_nuc * 10 + R_nuc]][j - 1] == 0) {
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
                                    if (NCflg == 1 && i2r[L2_nuc] != NucConst[i + 1]) {
                                        continue;
                                    }
                                    // if(chkDep2){continue:}
                                    for (unsigned int R2 = 0; R2 < pos2nuc[j - 1].size(); R2++) {
                                        int R2_nuc = pos2nuc[j - 1][R2];
                                        if (NCflg == 1 && i2r[R2_nuc] != NucConst[j - 1]) {
                                            continue;
                                        }

                                        if (DEPflg && Dep1[ii2r[L_nuc * 10 + L2_nuc]][i] == 0) {
                                            continue;
                                        }
                                        if (DEPflg && Dep1[ii2r[R2_nuc * 10 + R_nuc]][j - 1] == 0) {
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

                                        // check predefined hairpin energy
                                        // if((l == 5 || l == 6 || l == 8) && preHPN_flg == 1){
                                        //  if(predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].second != ""){
                                        //		if(NCflg == 1){
                                        //			string s1 = string(NucDef).substr(i, l);
                                        //			if(predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].second ==
                                        // s1){
                                        //				//C[ij][L][R] = MIN2(C[ij][L][R],
                                        // predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].first);
                                        // C[ij][L][R] = predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].first; // Note that
                                        // predefined hairpin is forced when it is found
                                        //			}
                                        //			}
                                        //		else{
                                        //			//一つ内側の塩基とのDependencyをチェックする。
                                        //			string s1 =
                                        // predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].second; 			int preL2_nuc
                                        // = n2i[s1[1]]; 			int preR2_nuc = n2i[s1[s1.size()-2]];
                                        //										//
                                        // cout << s1 << endl; 			if(DEPflg &&
                                        // Dep1[ii2r[L_nuc*10+preL2_nuc]][i] ==
                                        // 0){continue;} 			if(DEPflg && Dep1[ii2r[preR2_nuc*10+R_nuc]][j-1]
                                        // == 0){continue;} 			C[ij][L][R] =
                                        // predefHPN[i][l][i2r[L_nuc]][i2r[R_nuc]].first; // Note that predefined hairpin
                                        // is forced when it is found
                                        //		}
                                        //										//
                                        // exit(0);
                                        //	}
                                        //}
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

                                    int pq = getIndx(p, q, w_tmp, indx);

                                    for (unsigned int Lp = 0; Lp < pos2nuc[p].size(); Lp++) {
                                        int Lp_nuc = pos2nuc[p][Lp];
                                        if (NCflg == 1 && i2r[Lp_nuc] != NucConst[p]) {
                                            continue;
                                        }

                                        if (DEPflg && p == i + 1 && Dep1[ii2r[L_nuc * 10 + Lp_nuc]][i] == 0) {
                                            continue;
                                        }
                                        if (DEPflg && p == i + 2 && Dep2[ii2r[L_nuc * 10 + Lp_nuc]][i] == 0) {
                                            continue;
                                        }

                                        for (unsigned int Rq = 0; Rq < pos2nuc[q].size(); Rq++) { // nucleotide for p, q
                                            int Rq_nuc = pos2nuc[q][Rq];
                                            if (NCflg == 1 && i2r[Rq_nuc] != NucConst[q]) {
                                                continue;
                                            }

                                            if (DEPflg && q == j - 1 && Dep1[ii2r[Rq_nuc * 10 + R_nuc]][q] == 0) {
                                                continue;
                                            }
                                            if (DEPflg && q == j - 2 && Dep2[ii2r[Rq_nuc * 10 + R_nuc]][q] == 0) {
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
                                                if (NCflg == 1 && i2r[L2_nuc] != NucConst[i + 1]) {
                                                    continue;
                                                }

                                                if (DEPflg && Dep1[ii2r[L_nuc * 10 + L2_nuc]][i] == 0) {
                                                    continue;
                                                }

                                                for (int const R2_nuc : pos2nuc[j - 1]) {
                                                    if (NCflg == 1 && i2r[R2_nuc] != NucConst[j - 1]) {
                                                        continue;
                                                    }

                                                    if (DEPflg && Dep1[ii2r[R2_nuc * 10 + R_nuc]][j - 1] == 0) {
                                                        continue;
                                                    }

                                                    for (int const Lp2_nuc : pos2nuc[p - 1]) {
                                                        if (NCflg == 1 && i2r[Lp2_nuc] != NucConst[p - 1]) {
                                                            continue;
                                                        }

                                                        if (DEPflg && Dep1[ii2r[Lp2_nuc * 10 + Lp_nuc]][p - 1] == 0) {
                                                            continue;
                                                        }
                                                        if (p == i + 2 && L2_nuc != Lp2_nuc) {
                                                            continue;
                                                        } // check when a single nucleotide between i and p, this
                                                          // sentence confirm the dependency between Li_nuc and Lp2_nuc
                                                        if (DEPflg && i + 3 == p &&
                                                            Dep1[ii2r[L2_nuc * 10 + Lp2_nuc]][i + 1] == 0) {
                                                            continue;
                                                        } // check dependency between i+1, p-1 (i,X,X,p)

                                                        for (int const Rq2_nuc : pos2nuc[q + 1]) {
                                                            if (q == j - 2 && R2_nuc != Rq2_nuc) {
                                                                continue;
                                                            } // check when a single nucleotide between q and j,this
                                                              // sentence confirm the dependency between Rj_nuc and
                                                              // Rq2_nuc

                                                            if (NCflg == 1 && i2r[Rq2_nuc] != NucConst[q + 1]) {
                                                                continue;
                                                            }

                                                            if (DEPflg && Dep1[ii2r[Rq_nuc * 10 + Rq2_nuc]][q] == 0) {
                                                                continue;
                                                            }
                                                            if (DEPflg && q + 3 == j &&
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
                                if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                                    continue;
                                }

                                if (DEPflg && Dep1[ii2r[L_nuc * 10 + Li1_nuc]][i] == 0) {
                                    continue;
                                }

                                for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                                    int Rj1_nuc = pos2nuc[j - 1][Rj1];
                                    if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                        continue;
                                    }

                                    if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + R_nuc]][j - 1] == 0) {
                                        continue;
                                    }
                                    // if(DEPflg && j-i == 2 && i <= nuclen - 2 && Dep2[ii2r[L_nuc*10+R_nuc]][i] ==
                                    // 0){continue;}
                                    if (DEPflg && (j - 1) - (i + 1) == 2 &&
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
                            if (NCflg == 1 && i2r[Li1_nuc] != NucConst[i + 1]) {
                                continue;
                            }
                            if (DEPflg && Dep1[ii2r[L_nuc * 10 + Li1_nuc]][i] == 0) {
                                continue;
                            }

                            // int energy_M = M[indx[j]+i+1][Li1][R]+P->MLbase;
                            int energy_M = M[getIndx(i + 1, j, w_tmp, indx)][Li1][R] + P->MLbase;
                            M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
                        }

                        // create M[ij] from M[i][j-1]
                        for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                            int Rj1_nuc = pos2nuc[j - 1][Rj1];
                            if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                                continue;
                            }
                            if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + R_nuc]][j - 1] == 0) {
                                continue;
                            }

                            // int energy_M = M[indx[j-1]+i][L][Rj1]+P->MLbase;
                            int energy_M = M[getIndx(i, j - 1, w_tmp, indx)][L][Rj1] + P->MLbase;
                            M[ij][L][R] = MIN2(energy_M, M[ij][L][R]);
                        }

                        /* modular decomposition -------------------------------*/
                        for (int k = i + 2 + TURN; k <= j - TURN - 1; k++) { // Is this correct?
                            // cout << k << endl;
                            for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                                int Rk1_nuc = pos2nuc[k - 1][Rk1];
                                if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                                    continue;
                                }
                                // if(DEPflg && k == i + 2 && Dep1[ii2r[L_nuc*10+Rk1_nuc]][k-1] == 0){ continue;} //
                                // dependency between i and k - 1(=i+1) if(DEPflg && k == i + 3 &&
                                // Dep2[ii2r[L_nuc*10+Rk1_nuc]][k-1] == 0){ continue;} // dependency between i and k -
                                // 1(=i+2)

                                for (unsigned int Lk = 0; Lk < pos2nuc[k].size(); Lk++) {
                                    int Lk_nuc = pos2nuc[k][Lk];
                                    if (NCflg == 1 && i2r[Lk_nuc] != NucConst[k]) {
                                        continue;
                                    }
                                    if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                                        continue;
                                    } // dependency between k - 1 and k
                                    // if(DEPflg && (k-1) - i + 1 == 2 && Dep2[ii2r[Rk1_nuc*10+L_nuc]][k-1] == 0){
                                    // continue;} // dependency between i and k - 1

                                    // cout << i << " " << k-1 << ":" << M[indx[k-1]+i][L][Rk1] << "," << k << " " << j
                                    // << ":" << M[indx[j]+k][Lk][R] << endl; int energy_M =
                                    // M[indx[k-1]+i][L][Rk1]+M[indx[j]+k][Lk][R];
                                    int energy_M = M[getIndx(i, k - 1, w_tmp, indx)][L][Rk1] +
                                                   M[getIndx(k, j, w_tmp, indx)][Lk][R];
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
                for (unsigned int L = 0; L < 4; L++) {
                    DMl[j][L].fill(INF);
                    // AMW TODO: may be a better way to do this. Might want to set all of DMl to infs I think?
                    // fill(DMl[j][L], DMl[j][L] + 4, INF);
                }
            }
        }

        // Fill F matrix
        // F[1], as well as the rest of F, is default-initialized to 0.

        for (unsigned int L1 = 0; L1 < pos2nuc[1].size(); L1++) {
            int L1_nuc = pos2nuc[1][L1];
            if (NCflg == 1 && i2r[L1_nuc] != NucConst[1]) {
                continue;
            }

            for (int j = 2; j <= nuclen; j++) {

                //				int opt_flg_1 = 1;
                //				for(int I = 0; I < n_inter; I++){
                //					if(ofm[I] <= 1 && oto[I] >= 1){
                //						opt_flg_1 = 0;
                //						break;
                //					}
                //				}
                //				int opt_flg_j = 1;
                //				for(int I = 0; I < n_inter; I++){
                //					if(ofm[I] <= j && oto[I] >= j){
                //						opt_flg_j = 0;
                //						break;
                //					}
                //				}

                for (unsigned int Rj = 0; Rj < pos2nuc[j].size(); Rj++) {
                    int Rj_nuc = pos2nuc[j][Rj];
                    if (NCflg == 1 && i2r[Rj_nuc] != NucConst[j]) {
                        continue;
                    }

                    if (DEPflg && j == 2 && Dep1[ii2r[L1_nuc * 10 + Rj_nuc]][1] == 0) {
                        continue;
                    }
                    if (DEPflg && j == 3 && Dep2[ii2r[L1_nuc * 10 + Rj_nuc]][1] == 0) {
                        continue;
                    }

                    F[j][L1][Rj] = INF;

                    int type_L1Rj = BP_pair[i2r[L1_nuc]][i2r[Rj_nuc]];
                    if (type_L1Rj) {
                        //						if(opt_flg_1 && opt_flg_j){
                        int au_penalty = 0;
                        if (type_L1Rj > 2)
                            au_penalty = P->TerminalAU;
                        if (j <= w_tmp)
                            F[j][L1][Rj] =
                                MIN2(F[j][L1][Rj], C[getIndx(1, j, w_tmp, indx)][L1][Rj] + au_penalty); // recc 1
                        // F[j][L1][Rj] = MIN2(F[j][L1][Rj], C[indx[j] + 1][L1][Rj] + au_penalty); // recc 1
                        //						}
                    }

                    // create F[j] from F[j-1]
                    for (unsigned int Rj1 = 0; Rj1 < pos2nuc[j - 1].size(); Rj1++) {
                        int Rj1_nuc = pos2nuc[j - 1][Rj1];
                        if (NCflg == 1 && i2r[Rj1_nuc] != NucConst[j - 1]) {
                            continue;
                        }
                        if (DEPflg && Dep1[ii2r[Rj1_nuc * 10 + Rj_nuc]][j - 1] == 0) {
                            continue;
                        }

                        F[j][L1][Rj] = MIN2(F[j][L1][Rj], F[j - 1][L1][Rj1]); // recc 2
                    }

                    // create F[j] from F[k-1] and C[k][j]
                    // for (int k = 2; k <= j - TURN - 1; k++) { // Is this correct?
                    for (int k = MAX2(2, j - w_tmp + 1); k <= j - TURN - 1; k++) { // Is this correct?

                        //						int opt_flg_k = 1;
                        //						for(int I = 0; I < n_inter; I++){
                        //							if(ofm[I] <= k && oto[I] >= k){
                        //								opt_flg_k = 0;
                        //								break;
                        //							}
                        //						}

                        for (unsigned int Rk1 = 0; Rk1 < pos2nuc[k - 1].size(); Rk1++) {
                            int Rk1_nuc = pos2nuc[k - 1][Rk1];
                            if (NCflg == 1 && i2r[Rk1_nuc] != NucConst[k - 1]) {
                                continue;
                            }
                            if (DEPflg && k == 3 && Dep1[ii2r[L1_nuc * 10 + Rk1_nuc]][1] == 0) {
                                continue;
                            } // dependency between 1(i) and 2(k-1)
                            if (DEPflg && k == 4 && Dep2[ii2r[L1_nuc * 10 + Rk1_nuc]][1] == 0) {
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

                                int au_penalty = 0;
                                if (type_LkRj > 2)
                                    au_penalty = P->TerminalAU;
                                // int kj = indx[j] + k;
                                int kj = getIndx(k, j, w_tmp, indx);

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
            if (NCflg == 1 && i2r[L_nuc] != NucConst[1]) {
                continue;
            }
            for (unsigned int R = 0; R < pos2nuc[nuclen].size(); R++) {
                int R_nuc = pos2nuc[nuclen][R];
                if (NCflg == 1 && i2r[R_nuc] != NucConst[nuclen]) {
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

        if (rand_tb_flg) {
            // Fill F2 matrix
            for (int l = 5; l <= nuclen; l++) {
                if (l > w_tmp)
                    break;
                cout << "process F2:" << l << endl;

                for (int i = 1; i <= nuclen - l + 1; i++) {
                    int j = i + l - 1;

                    for (unsigned int L = 0; L < pos2nuc[i].size(); L++) {
                        int L_nuc = pos2nuc[i][L];

                        for (unsigned int R = 0; R < pos2nuc[j].size(); R++) {
                            int R_nuc = pos2nuc[j][R];
                            int ij = getIndx(i, j, w_tmp, indx);

                            F2[ij][L][R] = 0;

                            int type = BP_pair[i2r[L_nuc]][i2r[R_nuc]];

                            // from i, j-1 -> i, j
                            for (unsigned int R1 = 0; R1 < pos2nuc[j - 1].size(); R1++) {
                                int R1_nuc = pos2nuc[j - 1][R1];
                                if (DEPflg && Dep1[ii2r[R1_nuc * 10 + R_nuc]][j - 1] == 0) {
                                    continue;
                                }
                                int ij1 = getIndx(i, j - 1, w_tmp, indx);
                                F2[ij][L][R] = MIN2(F2[ij][L][R], F2[ij1][L][R1]);
                            }
                            // from i-1, j -> i, j
                            for (unsigned int L1 = 0; L1 < pos2nuc[i + 1].size(); L1++) {
                                int L1_nuc = pos2nuc[i + 1][L1];
                                if (DEPflg && Dep1[ii2r[L_nuc * 10 + L1_nuc]][i] == 0) {
                                    continue;
                                }
                                int i1j = getIndx(i + 1, j, w_tmp, indx);
                                F2[ij][L][R] = MIN2(F2[ij][L][R], F2[i1j][L1][R]);
                            }

                            // from C
                            int au_penalty = 0;
                            if (type > 2)
                                au_penalty = P->TerminalAU;
                            if (j - i + 1 <= w_tmp) {
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
                                        if (DEPflg && Dep1[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
                                            continue;
                                        } // dependency between k - 1 and k

                                        int energy = F2[getIndx(i, k - 1, w_tmp, indx)][L][Rk1] +
                                                     F2[getIndx(k, j, w_tmp, indx)][Lk][R];
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

        if (rand_tb_flg) {
            backtrack2(&optseq, &*sector, base_pair, C, M, F2, indx, minL, minR, P, NucConst, pos2nuc, NCflg, i2r,
                       nuclen, w_tmp, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, DEPflg, predefHPN, predefHPN_E, substr,
                       n2i, NucDef);
        } else {
            backtrack(&optseq, &*sector, base_pair, C, M, F, indx, minL, minR, P, NucConst, pos2nuc, NCflg, i2r,
                      nuclen, w_tmp, BP_pair, i2n, rtype, ii2r, Dep1, Dep2, DEPflg, predefHPN, predefHPN_E, substr, n2i,
                      NucDef);
        }

        //塩基Nの修正
        for (int i = 1; i <= nuclen; i++) {
            if (optseq[i] == 'N') {
                for (unsigned int R = 0; R < pos2nuc[i].size();
                     R++) { // check denendency with the previous and next nucleotide
                    int R_nuc = pos2nuc[i][R];
                    if (NCflg == 1 && i2r[R_nuc] != NucConst[i]) {
                        continue;
                    }

                    if (i != 1 && optseq[i - 1] != 'N') { // check consistensy with the previous nucleotide
                        int R_prev_nuc = n2i[optseq[i - 1]];
                        if (DEPflg && Dep1[ii2r[R_prev_nuc * 10 + R_nuc]][i - 1] == 0) {
                            continue;
                        }
                    }
                    if (i != nuclen && optseq[i + 1] != 'N') { // check consistensy with the next nucleotide
                        int R_next_nuc = n2i[optseq[i + 1]];
                        if (DEPflg && Dep1[ii2r[R_nuc * 10 + R_next_nuc]][i] == 0) {
                            continue;
                        }
                    }
                    if (i < nuclen - 1 && optseq[i + 2] != 'N') { // check consistensy with the next nucleotide
                        int R_next_nuc = n2i[optseq[i + 2]];
                        if (DEPflg && Dep2[ii2r[R_nuc * 10 + R_next_nuc]][i] == 0) {
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
            char aa = codon_table.c2a(n2i[optseq[i]], n2i[optseq[i + 1]], n2i[optseq[i + 2]]);
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

        if (part_opt_flg == 1) {
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

                string part_optseq = rev_fold_step1(&part_aaseq[0], part_aalen, codon_table, exc);
                rev_fold_step2(&part_optseq, &part_aaseq[0], part_aalen, codon_table, exc);
                // combine optseq_rev and optseq
                j = 1;
                for (int i = ofm[I]; i <= oto[I]; i++) {
                    optseq[i] = part_optseq[j++];
                }
            }
            fixed_fold(optseq, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
            // fixed_fold(optseq, indx, w_tmp, predefHPN_E, BP_pair, P, aaseq, codon_table);
        }

        if (m_disp) {
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

    } while (all_aaseq.next());

    clock_t end = clock();
    float sec = (double)end / CLOCKS_PER_SEC;
    float min = sec / 60;
    cout << "Runing time: " << min << " minutes" << endl;

    //	struct rusage r;
    //	if (getrusage(RUSAGE_SELF, &r) != 0) {
    //			/*Failure*/
    //	}
    //	printf("Memory usage: %ld Mb\n", r.ru_maxrss/1024);

    //	printf("Memory usage: %ld Mb\n", r.ru_maxrss/1024);

    return 0;
}

void backtrack(string *optseq, stack *sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f,
               vector<int> const & indx, const int &initL, const int &initR, paramT *const P, const vector<int> &NucConst,
               const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
               int const (&BP_pair)[5][5], array<char, 20> const &i2n, int *const &rtype, array<int, 100> const &ii2r,
               vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, int &DEPflg,
               vector<vector<vector<vector<pair<int, string>>>>> &, map<string, int> &predefE,
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

void backtrack2(string *optseq, stack *sector, vector<bond> & base_pair, vector<vector<vector<int>>> const &c, vector<vector<vector<int>>> const &m, vector<vector<vector<int>>> const &f2,
                vector<int> const & indx, const int &initL, const int &initR, paramT *const P, const vector<int> &NucConst,
                const vector<vector<int>> &pos2nuc, const int &NCflg, array<int, 20> const &i2r, int const &length, int const &w,
                int const (&BP_pair)[5][5], array<char, 20> const &i2n, int *const &rtype, array<int, 100> const &ii2r,
                vector<vector<int>> &Dep1, vector<vector<int>> &Dep2, int &DEPflg,
                vector<vector<vector<vector<pair<int, string>>>>> &, map<string, int> &predefE,
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
    stack sector[500];
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

void fixed_fold(string optseq, vector<int> const & indx, const int &w, map<string, int> &predefE, const int (&BP_pair)[5][5],
                paramT *P, char *aaseq, const codon& codon_table) {
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


vector<vector<int>> getPossibleNucleotide(char *aaseq, int aalen, codon &codon_table, map<char, int> &n2i,
                                          string excludedCodons) {

    vector<vector<int>> v;

    int nuclen = aalen * 3;
    v.resize(nuclen + 1);

    for (int i = 0; i < aalen; i++) {
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
                v[nucpos].push_back(n2i['V']);
                v[nucpos].push_back(n2i['W']);
            } else if (aa == 'R' && k == 1) {
                v[nucpos].push_back(n2i['X']);
                v[nucpos].push_back(n2i['Y']);
            } else {
                bool flg_A, flg_C, flg_G, flg_U;
                flg_A = flg_C = flg_G = flg_U = false;

                for (unsigned int j = 0; j < codons.size(); j++) { // for each codon corresp. to each aa
                    char nuc = codons[j][k];
                    if (nuc == 'A' && flg_A == 0) {
                        v[nucpos].push_back(n2i[nuc]);
                        flg_A = 1;
                    } else if (nuc == 'C' && flg_C == 0) {
                        v[nucpos].push_back(n2i[nuc]);
                        flg_C = 1;
                    } else if (nuc == 'G' && flg_G == 0) {
                        v[nucpos].push_back(n2i[nuc]);
                        flg_G = 1;
                    } else if (nuc == 'U' && flg_U == 0) {
                        v[nucpos].push_back(n2i[nuc]);
                        flg_U = 1;
                    }
                }
            }
            nucpos -= k;
        }
    }
    return v;
}
