#include <iostream>
#include <fstream>
#include <cstring> // to use strcpy for gcc-4.3 or later

#include "fasta.hpp"

fasta::fasta(const char *fname) { // DP$B%^%H%j%/%9$NF0E*%a%b%j3NJ](B

    initP();
    numSeq = 0;

    ifstream ifs(fname);
    if (ifs) {
        string line;
        string tmp_desc;
        string tmp_seq;
        // int id = 0;
        getline(ifs, line);
        if (line[0] != '>') {
            cerr << "Invalid fasta format." << endl;
            exit(1);
        } else {
            line.erase(0, 1); // ">"$B$r:o=|(B
            tmp_desc = line;
        }
        while (getline(ifs, line)) { // $B9T$NFI$_9~$_(B
            // cout << line << endl;
            if (line[0] == '>') {
                eachseq e;
                // cout << "ok1" << endl;
                e.seq = tmp_seq;
                e.desc = tmp_desc;
                data.push_back(e); //
                numSeq++; // $BG[Ns?t(B

                line.erase(0, 1); // ">"$B$r:o=|(B
                tmp_desc = line;
            } else {
                tmp_seq += line;
            }
        }
        // $B0lHV:G8e$NG[Ns(B
        eachseq e;
        e.seq = tmp_seq;
        e.desc = tmp_desc;
        data.push_back(e);
    } else {
        if (fname == nullptr) {
            cerr << "Error: no input file" << endl;
        } else {
            cerr << "Error: cannot open file(" << fname << ")" << endl;
        }
        exit(1);
    }
}

fasta::~fasta() {}

void fasta::printP() { cout << p << endl; }
