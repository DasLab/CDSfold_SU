#include <iostream>
#include <fstream>
#include <cstring> // to use strcpy for gcc-4.3 or later

#include "fasta.hpp"

fasta::fasta(const char *fname) { // DPマトリクスの動的メモリ確保

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
            line.erase(0, 1); // ">"を削除
            tmp_desc = line;
        }
        while (getline(ifs, line)) { // 行の読み込み
            // cout << line << endl;
            if (line[0] == '>') {
                eachseq e;
                // cout << "ok1" << endl;
                e.seq = tmp_seq;
                e.desc = tmp_desc;
                data.push_back(e); //
                numSeq++; // 配列数

                line.erase(0, 1); // ">"を削除
                tmp_desc = line;
            } else {
                tmp_seq += line;
            }
        }
        // 一番最後の配列
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
