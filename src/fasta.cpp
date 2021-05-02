#include <iostream>
#include <fstream>
#include <cstring> // to use strcpy for gcc-4.3 or later

#include "fasta.hpp"

/* Constructor for the fasta class representing an amino acid
 * sequence or an RNA strand
 */
fasta::fasta(const char *fname) {

    initP();               /* set index into sequences to 0 */
    numSeq_ = 0;            /* initialize number of sequences */

    ifstream ifs(fname);
   
    /* file can be opened - try to parse it */
    if (ifs) {
        string line;
        string tmp_desc;
        string tmp_seq;
       
        /* read the first line, which should contain a comment starting w/ '<' */
        getline(ifs, line);
        if (line[0] != '>') {
            cerr << "Invalid fasta format." << endl;
            exit(1);
        } else {
            line.erase(0, 1);           /* drop the comment char '<' */
            tmp_desc = line;            /* save the description */
        }
        
        /* repeatedly loop reading sequences from the file */
        while (getline(ifs, line)) {
             
            if (line[0] == '>') {        /* hit another description line - save the sequence */
                eachseq e;
                e.seq_ = tmp_seq;
                e.desc_ = tmp_desc;       /* save the description we've built up */
                data_.push_back(e);
                numSeq_++; 

                line.erase(0, 1);
                tmp_desc = line;
            } else {                   /* line w/ the sequence - build up the sequence */
                tmp_seq += line;
            }
        }
        /* finished parsing - save the rest of the sequence */ 
        eachseq e;
        e.seq_ = tmp_seq;
        e.desc_ = tmp_desc;
        data_.push_back(e);
    
    /* file could not be opened */ 
    } else {
        if (fname == nullptr) {
            cerr << "Error: no input file" << endl;
        } else {
            cerr << "Error: cannot open file(" << fname << ")" << endl;
        }
        exit(1);
    }
}


void fasta::printP() { cout << p_ << endl; }
