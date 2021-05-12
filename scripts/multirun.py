# multirun.py
#
# Script for running CDS fold with a given set of arguments multiple times
# to generate a variety of possible sequences
#
# The arguments for running CDSfold are specified by two formats. The short
# form is a list of lists, with each sublist holding the arguments as a
# single string, path to the .faa file, and the number of times to run:
# [["-w 20 -e "CCA", "example/random_seq.faa", 3]]
#
# The full form has a single list for each run: 
# [["-w 20 -e CCA", "example/random_seq.faa"],
#  ["-w 20 -e CCA", "example/random_seq.faa"],
#  ["-w 20 -e CCA", "example/random_seq.faa"]]


import sys
import os
import argparse
import util

CDS_HOME = os.getenv("CDS_HOME")

# list of CDSfold commands to run if the -l argument is specified
RUN_SET = [["-w 20 -e CCA", "example/mev.faa", 3],
           ["-w 10", "example/mev.faa", 1],
           ["", "example/mev.faa", 10]]

def expand_run_set(run_set):
    '''Expands a passed run_set in short form to the format expected by 
    run_cds_fold_many and returns the result'''

    # expand the run set to the format expected by run_cds_fold_many
    full_run_set = []
    for run in run_set:
        run_spec = run[0:-1]
        for i in range(run[-1]):
            full_run_set.append(run_spec)

    return full_run_set

def get_parser():
    '''create a parser for the CDSfold multirun script'''
    parser = argparse.ArgumentParser(description="Script for running CDSfold multiple times.")

    parser.add_argument("cds_args", nargs='?', type=str, 
                        help='quote enclosed CDSfold options', default="")
    parser.add_argument("faa_path", nargs='?', type=str, 
                        help='path to faa file relative to CDS_HOME', default="")
    parser.add_argument("--repeat", help="number of times to run CDSfold")
    parser.add_argument("-l", action="store_true", 
                        help="use list specified in the script instead of cmd line")
    
    return parser

if __name__ == "__main__":
    
    # location of cdsfold binary to use 
    test_info = {
        "name": "default",
        "bin_path": os.path.join(CDS_HOME, "src/CDSfold"),
        "out_path": None
    } 
    
    parser = get_parser()
    args = parser.parse_args() 

    # use the run set at the top of the file
    if args.l:
        run_set = RUN_SET
    # use the command line arguments
    else:
        if not args.repeat:
            repeat = 1
        else:
            repeat = args.repeat
        run_set = [[args.cds_args, args.faa_path, repeat]]


    full_run_set = expand_run_set(run_set)                        # expand to full form
    raw_results = util.run_cdsfold_many(full_run_set, test_info)  # run CDSfold
    results = util.clean_results(raw_results, full_run_set)       # clean the results
    print(results)
