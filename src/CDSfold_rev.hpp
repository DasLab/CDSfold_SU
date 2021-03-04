#pragma once

#include "codon.hpp"
#include "Util.hpp"

extern "C" {
#include "energy_const.h"
}

struct Ntable {
    int A;
    int C;
    int G;
    int U;

  public:
    Ntable() {
        A = 0;
        C = 0;
        G = 0;
        U = 0;
    }
};

struct Ctable {
    int AU;
    int GC;
    int GU;

  public:
    Ctable() {
        AU = 0;
        GC = 0;
        GU = 0;
    }
};

inline void addNtable(Ntable &N, char n) {
    switch (n) {
    case 'A':
        N.A++;
        break;
    case 'C':
        N.C++;
        break;
    case 'G':
        N.G++;
        break;
    case 'U':
        N.U++;
        break;
    default:
        cerr << "Unexpected nucleotide" << n << endl;
    }
}
inline void subtNtable(Ntable &N, char n) {
    switch (n) {
    case 'A':
        N.A--;
        break;
    case 'C':
        N.C--;
        break;
    case 'G':
        N.G--;
        break;
    case 'U':
        N.U--;
        break;
    default:
        cerr << "Unexpected nucleotide" << n << endl;
    }
}

inline void addCtable(Ctable &C, char n1, char n2) {
    if ((n1 == 'A' && n2 == 'U') || (n1 == 'U' && n2 == 'A')) {
        C.AU++;
    } else if ((n1 == 'G' && n2 == 'C') || (n1 == 'C' && n2 == 'G')) {
        C.GC++;
    } else if ((n1 == 'G' && n2 == 'U') || (n1 == 'U' && n2 == 'G')) {
        C.GU++;
    }
}

inline void showNtable(Ntable N) {
    cout << "A=" << N.A << endl;
    cout << "C=" << N.C << endl;
    cout << "G=" << N.G << endl;
    cout << "U=" << N.U << endl;
}
inline void showCtable(Ctable C) {
    cout << "AU=" << C.AU << endl;
    cout << "GC=" << C.GC << endl;
    cout << "GU=" << C.GU << endl;
}

inline auto calcPseudoEnergy(const Ntable &N, const Ctable &C) -> float {
    float energy = 0;

    energy = N.A * N.U * -1;
    energy += N.G * N.C * -3.12;
    energy += N.G * N.U * -1;

    //	cout << "chk1:" << energy << endl;

    energy -= C.AU * -1;
    //	cout << "chk2-1:" << energy << endl;
    energy -= C.GC * -3.12;
    //	cout << "chk2-2:" << energy << ":" << C.GC << endl;
    energy -= C.GU * -1;
    //	cout << "chk2-3:" << energy << endl;

    return energy;
}


inline auto countNtable(string &seq, int F) -> Ntable {
    Ntable N;
    N.A = 0;
    N.C = 0;
    N.G = 0;
    N.U = 0;
    for (unsigned int i = F; i < seq.size(); i++) {
        switch (seq[i]) {
        case 'A':
            N.A++;
            break;
        case 'C':
            N.C++;
            break;
        case 'G':
            N.G++;
            break;
        case 'U':
            N.U++;
            break;
        default:
            cerr << "Unexpected nucleotide" << seq[i] << endl;
        }
    }
    return N;
}

inline auto countCtable(string &seq, int F) -> Ctable {
    Ctable C;
    C.AU = 0;
    C.GC = 0;
    C.GU = 0;

    for (unsigned int i = F; i < seq.size() - 1; i++) {
        for (unsigned int j = i + 1; j <= MIN2(i + 3, seq.size() - 1); j++) {
            // cout << i << ";" << j << endl;
            addCtable(C, seq[i], seq[j]);
        }
    }

    return C;
}

void rev_fold_step2(string *optseq_r, const char *aaseq, const int aalen, codon &codon_table,
                    const string &exc_codons);