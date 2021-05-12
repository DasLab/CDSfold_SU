# various utility functions related to CDSfold

import os
from os import path
import sys
import subprocess

CDS_HOME = os.getenv("CDS_HOME")
FOLD_OFFSET = -1          # offsets of the the fold, seq, and aa_seq relative
SEQ_OFFSET = -2           # to the line w/ MFE
AA_SEQ_OFFSET = -3

RUNTIME_LINE_NUM = -2     # line number of output w/ run time

def run_cdsfold_many(test_suite, test_info, quiet=False, save=False):
    """
    run all tests comparing the behavior of the latest CDSfold binary with either
    a saved gold standard output or the results of running a reference binary

    inputs: test_suite (list) - list of tests to run
            test_info (dict) - dictionary with name of tests, path to binary, and
                               output path
            quiet (bool) - suppress prints
            save (bool) - save outputs to file
    """

    # results dictionary mapping name of binary run to list of output strings,
    # one string for each test
    results = {}
    results[test_info["name"]] = []

    if not quiet: 
        print("\n==== Executing: {} ==== ".format(test_info["name"])) 
    
    for i, test in enumerate(test_suite):
        if not quiet: 
            print("Running test: {}".format(test))
       
        # set up the command for running the test
        if len(test[0]) == 0:
            # no command line args
            cmd = [test_info["bin_path"], 
                   os.path.join(CDS_HOME, test[1])]
        else:
            # command line args used
            cmd = [test_info["bin_path"],
                   test[0],
                   os.path.join(CDS_HOME, test[1])]
        
        # run command and capture output
        output = subprocess.run(cmd, capture_output=True)
        results[test_info["name"]].append(output.stdout.decode("utf-8"))
       
        # save output to file
        if save:
            # file to save output in 
            save_path = test_info["out_path"].format(i)
            
            # delete file it if already exists
            if os.path.exists(save_path): 
                os.remove(save_path)  # clear the file being output to
    
            if not quiet:
                print("Saving test output to {}".format(save_path))
            
            with open(save_path, "w") as f:
                # result is the last added to list of results
                f.writelines(["Executed command: {}\n".format(str(test)), 
                              results[test_info["name"]][-1]])
    
    # dump the results to file if specified by command line arguments

    return results

def parse_mfe(mfe_line):
    '''parses a line with the mean free energy and returns the
    MFE in kcal/mol as a float'''
    mfe = mfe_line.split()[0]
    return float(mfe[4:])

def parse_runtime(runtime_line):
    '''parses a line with the runtime and returns the runtime in
    seconds as a float ''' 
    runtime_list = runtime_line.split()
    # runtime is the element before "minutes" 
    runtime = runtime_list[runtime_list.index("minutes") - 1]
    
    return float(runtime) * 60 

def clean_results(raw_results, test_suite, full=True):
    '''clean the results from run_all_tests. Converts the raw_results to
    a dictionary of the form:
    {
        "test_bin": [
        {
            'cmd': "< command line arguments CDSfold invoked with >", 
            'run_time': [float, float, float, ...],  # run time for each seq
            'sequence': [str, str, str, ...],  # each sequence from this file
            'fold':     [str, str, str, ...],  # each fold pattern
            'lines': []
        },
        {
            test_1 info...
        }
        ]
    }
    inputs: raw_results - results from run_cdsfold_many
            test_suite - full form list of commands run
            full (bool) - save the entire cdsfold output to dictionary
    '''

    results = {}

    # iterate through the test binaries
    for key in raw_results.keys():
        results[key] = []
        # iterate through tests with the same binary
        for i in range(len(test_suite)):
            test_result = {"sequences": [], "folds": [], "mfes": [], "aa_seqs": []} 
            test_result["cmd"] = str(test_suite[i])
            lines = raw_results[key][i].split("\n")
            
            if full:
                # save all the output lines only if specified
                test_result["lines"] = lines
            
            # potentially multiple sequences in same CDSfold run - parse for all the 
            # sequence specific information
            for i, line in enumerate(lines):
                if "MFE:" in line: 
                    test_result["mfes"].append(parse_mfe(lines[i]))
                    test_result["folds"].append(lines[i + FOLD_OFFSET])
                    test_result["sequences"].append(lines[i + SEQ_OFFSET])
                    test_result["aa_seqs"].append(lines[i + AA_SEQ_OFFSET])
            
            # parse the runtime for the entire run 
            test_result["Run time (s)"] = parse_runtime(lines[RUNTIME_LINE_NUM])
            results[key].append(test_result)

    return results
