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

#include "CDSfold.hpp"
#include "CDSfold_rev.hpp"
#include "fasta.hpp"
#include "Problem.hpp"
#include "Options.hpp"

//#define MAXLOOP 20
#define noGUclosure 0

using namespace std;
auto main(int argc, char *argv[]) -> int {
    // printf("%d\n%ld", INT_MAX, LONG_MAX);
    // get options
    Options options;
    {
        int opt;
        while ((opt = getopt(argc, argv, "w:e:f:t:rMUR")) != -1) {
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
            }
        }
    }
    // exit(0);

    // -R Check for options
    if (options.random_backtrack) {
        if (options.max_bp_distance != 0 || options.codons_excluded != "" || options.show_memory_use || options.estimate_memory_use || options.maximize_mfe || options.partial_opt) {
            cerr << "The -R option must not be used together with other options." << endl;
            return 0;
        }
    }

    if (options.partial_opt) {
        // Partial optimization was specified.
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

    fasta all_aaseq(argv[optind]); // get all sequences

    cout << "W = " << options.max_bp_distance << endl;
    cout << "e = " << options.codons_excluded << endl;
    do {
        string aaseq = string(all_aaseq.getSeq());

        // OK, here's where we begin.
        Problem problem(options, aaseq);
        problem.calculate();


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
