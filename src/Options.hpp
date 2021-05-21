/* Options.hpp
 *
 * Struct for holding different command line options; populated by the main
 * function in CDSfold by parsing the command line arguments. Please refer to the
 * README for details on each command line argument
 */



#pragma once

#include <string>
#include "constants.hpp"

struct Options {
    unsigned int max_bp_distance = 0;         // -w  max distance between base pairs
    std::string codons_excluded = "";         // -e  comma separated codons to avoid
    bool show_memory_use = false;             // -M
    bool estimate_memory_use = false;         // -U
    bool random_backtrack = false;            // -R
    bool maximize_mfe = false;                // -r  perform heuristic pseudo-MFE max
    bool partial_opt = false;                 // -f and -t
    unsigned int opt_from = 0;                // -f  start position for stable/unstable
    unsigned int opt_to = 0;                  // -t  end position for stable/unstable
    bool DEPflg = true;
    bool nucleotide_constraints = false;
    bool fixed_seed = false;
    float temp = DEFAULT_TEMP;                // -C temperature in Celcius
    float jitter = 0;                         // -j range (0 to 1)
};
