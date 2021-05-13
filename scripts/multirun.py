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
import json

CDS_HOME = os.getenv("CDS_HOME")

# list of CDSfold commands to run if the -l argument is specified
RUN_SET = [["-w 20 -e CCA", "example/mev.faa", 3],
           ["-w 10", "example/mev.faa", 1],
           ["-w 20", "example/mev.faa", 1],
           ["-w 30 -e CCA, UAG", "example/mev.faa", 1],
           ["-w 40", "example/mev.faa", 1],
           ["", "example/random_seq.faa", 1]]

def expand_run_set(run_set):
    '''Expands a passed run_set in short form to the format expected by 
    run_cds_fold_many and returns the result'''

    # expand the run set to the format expected by run_cds_fold_many
    full_run_set = []
    for run in run_set:
        run_spec = run[0:-1]
        for i in range(int(run[-1])):
            full_run_set.append(run_spec)

    return full_run_set

def dump_eterna_schema(results, outname):
    '''dump the results of the run to a json file according to the eterna schema'''
    # reformat to group sequences that encode the same amino acid chain
    sequences = {}      # map from sequence to puzzle id - puzzle IDs assigned sequentially
    eterna_schema = []  # eterna puzzle sequence to dump
    
    # look for all unique AA sequences and add empty dicts for each to eterna_schema 
    for key in results.keys():
        for runs in results[key]:
            # iterate through amino acids in this .faa file
            for i, seq in enumerate(runs["aa_seqs"]):
                # new amino acid sequence found
                if seq not in sequences.keys():
                    sequences[seq] = len(sequences)
                    # add empty puzzleID entry corresponding to new AA sequence
                    new_design = {"puzzleID": sequences[seq], "designs": []}
                    eterna_schema.append(new_design)

    # loop through again, now adding the sequences to eterna_schema
    for key in results.keys():
        for run in results[key]:
            # iterate through amino acids in this .faa file
            for i, seq in enumerate(run["aa_seqs"]):
                puzzle_dict = eterna_schema[sequences[seq]]   # lookup index of aa sequence in eterna
                new_seq = {
                    "sequence": run["sequences"][i], 
                    "title": util.cmd_to_str(run["cmd"]), 
                    "body": ""
                }
                puzzle_dict["designs"].append(new_seq)

    data = json.dumps(eterna_schema, indent=4)
    with open(outname, "w") as f:
        f.write(data)

    return

def dump_results(results, outname):
    '''dump the results of the run to a json file'''
    
    outname = outname.replace(".json", "_full.json")
    data = json.dumps(results, indent=4) 
    with open(outname, "w") as f:
        f.write(data)

    return

def get_parser():
    '''create a parser for the CDSfold multirun script'''
    parser = argparse.ArgumentParser(description="Script for running CDSfold multiple times.")

    parser.add_argument("cds_args", nargs='?', type=str, 
                        help='quote enclosed CDSfold options', default="")
    parser.add_argument("faa_path", nargs='?', type=str, 
                        help='path to faa file relative to CDS_HOME', default="")
    parser.add_argument("--repeat", type=int, help="number of times to run CDSfold")
    parser.add_argument("--out", help="name of output file - requires .json extension")
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

    if not args.out:
        outname = "cds_multirun.json"
    else:
        outname = args.out

    full_run_set = expand_run_set(run_set)                        # expand to full form
    raw_results = util.run_cdsfold_many(full_run_set, test_info)  # run CDSfold
    results = util.clean_results(raw_results, full_run_set, full=False) # clean the results
    dump_results(results, outname)           # save full results
    dump_eterna_schema(results, outname)
