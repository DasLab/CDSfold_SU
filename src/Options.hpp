#pragma once

#include <string>

struct Options {
    unsigned int max_bp_distance = 0;            // -w
    std::string codons_excluded = "";      // -e
    bool show_memory_use = false;       // -M
    bool estimate_memory_use = false;      // -U
    bool random_backtrack = false;  //-R
    bool maximize_mfe = false;      // -r
    bool partial_opt = false; // -f and -t
    unsigned int opt_from = 0;       // -f
    unsigned int opt_to = 0;       // -t
    bool DEPflg = true;
    bool nucleotide_constraints = false;
};
