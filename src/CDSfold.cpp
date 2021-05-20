/*
 * CDSfold.cpp

 *
 *  Created on: Sep 2, 2014
 *      Author: terai
 */
#define TURN 3
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <unistd.h>

extern "C" {
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
}

#include "CDSfold_rev.hpp"
#include "fasta.hpp"
#include "Problem.hpp"
#include "Options.hpp"

#define noGUclosure 0

using namespace std;
auto main(int argc, char *argv[]) -> int {
    // printf("%d\n%ld", INT_MAX, LONG_MAX);
    // get options
    Options options;
    
    // TODO this doesn't need to be in brackets
    {
        int opt;
        while ((opt = getopt(argc, argv, "w:e:f:j:t:C:rMURs")) != -1) {
            switch (opt) {
            case 'w':
                options.max_bp_distance = atoi(optarg);
                break;
            case 'e':
                options.codons_excluded = string(optarg);
                break;
            case 'M':
                options.show_memory_use = true;
                break;
            case 'U':
                options.estimate_memory_use = true;
                break;
            case 'R':
                options.random_backtrack = true;
                break;
            case 'r':
                options.maximize_mfe = true;
                break;
            case 'f':
                options.opt_from = atoi(optarg);
                options.partial_opt = true;
                break;
            case 't':
                options.opt_to = atoi(optarg);
                options.partial_opt = true;
                break;
            case 's':
                options.fixed_seed = true;
                break;
            case 'C':
                options.temp = atof(optarg);
                break;
            case 'j':
                options.jitter = atof(optarg);
                break;
            }
        }
    }
    // exit(0);

    // Check that -R option isn't used with other options
    if (options.random_backtrack) {
        if (options.max_bp_distance != 0 || options.codons_excluded != "" || options.show_memory_use || options.estimate_memory_use || options.maximize_mfe || options.partial_opt) {
            cerr << "The -R option must not be used together with other options." << endl;
            return 0;
        }
    }

    // Check that partial optimization options are used correctly
    if (options.partial_opt) {
        if (options.opt_from == 0 || options.opt_to == 0) {
            cerr << "The -f and -t option must be used together." << endl;
            exit(1);
        }
        if (options.opt_from < 1) {
            cerr << "The -f value must be 1 or more." << endl;
            exit(1);
        }
        if (options.opt_to < 1) {
            cerr << "The -t value must be 1 or more." << endl;
            exit(1);
        }
        if (options.opt_to < options.opt_from) {
            cerr << "The -f value must be smaller than -t value." << endl;
            exit(1);
        }
    }

    // Check that jitter is in correct range
    if (options.jitter >= 1.0 or options.jitter < 0.0) {
        cerr << "The -j option must be between 0 and 1" << endl;
        exit(1);
    }

    fasta all_aaseq(argv[optind]);       // get all sequences
    cout << "-w [max base pair distance] = " << options.max_bp_distance << endl;   // display max distance
    cout << "-e [excluded codons]        = " << options.codons_excluded << endl;   // display excluded codons

    do {
        string aaseq = string(all_aaseq.getSeq());

        // OK, here's where we begin.
        Problem problem(options, aaseq);
        problem.calculate();


    } while (all_aaseq.next());

    clock_t end = clock();
    float sec = (double)end / CLOCKS_PER_SEC;
    float min = sec / 60;
    cout << "Running time: " << min << " minutes" << endl;

    //	struct rusage r;
    //	if (getrusage(RUSAGE_SELF, &r) != 0) {
    //			/*Failure*/
    //	}
    //	printf("Memory usage: %ld Mb\n", r.ru_maxrss/1024);

    //	printf("Memory usage: %ld Mb\n", r.ru_maxrss/1024);

    return 0;
}
