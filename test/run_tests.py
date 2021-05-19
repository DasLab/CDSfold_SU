import os
from os import path
import subprocess
import argparse
import sys

CDS_HOME = os.getenv("CDS_HOME")
sys.path.append(os.path.join(CDS_HOME, "scripts"))

import util


def diff_output(results, args):
    '''diffs the results from the reference and CDS binaries.
    inputs: results - dictionary storing the results from running each binary
            test_suite - list of tests run'''
 
    if not args.q:
        print("\n==== Comparing reference and cds results ====")

    if len(results) != 2:
        raise ValueError("Expected exactly 2 elements in results; found {}".format(len(results)))

    mismatches = 0                   # number of mismatched tests

    # iterate through the tests
    for i, test in enumerate(test_suite):
        cds_result = results["cds"][i]   # convert string to list of lines
        ref_result = results["ref"][i] 

        # iterate through lines in result
        for j in range(min(len(cds_result["lines"]), len(ref_result["lines"]))):
            # lines don't match 
            if cds_result["lines"][j] != ref_result["lines"][j]:
                # non-matching run times do not count as mismatches 
                if "Running time" in cds_result["lines"][j] and \
                   "Running time" in ref_result["lines"][j]:
                    continue
                else:
                    print("Test {} failed. Command: {}".format(i, cds_result["cmd"]))
                    mismatches += 1
                    break
    
    print("\n========= Test Summary ===========")
    print("Results match reference for {} out of {} tests\n".format(
          len(test_suite) - mismatches, len(test_suite)))
  
    if mismatches > 0:
        raise ValueError("{} tests failed to match reference".format(mismatches))

    return

def get_parser():
    '''creates a parse for parsing command line arguments'''

    # create parser
    parser = argparse.ArgumentParser(description="Test harness for CDSfold")
    
    # use full test suite
    parser.add_argument("-f", action="store_true", help="use full test suite")

    # save the test results to a file
    parser.add_argument("-s", action="store_true", help="save output results to file")

    # run in verbose mode
    parser.add_argument("-q", action="store_true", help="verbose mode")
    
    return parser

def print_summary(args, num_tests):
    '''prints a summary of the tests and arguments
    inputs: args - args object with the parsed cmd line arguments
            num_tests - number of tests in the test suite'''

    summary_str = "\n=== Starting CDSFold test suite. Running {} tests -".format(num_tests)
    
    if args.q:
        summary_str += " Quiet mode -"
    else:
        summary_str += " Verbose mode -"
    if args.s:
        summary_str += " Saving outputs"

    summary_str += " ==="
    print(summary_str)

    return

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    # list of simple tests
    simple_tests = [
        ['', 'example/mev.faa'],
        ['-w 20', 'example/mev.faa'],
        ['-e GGU', 'example/mev.faa'],
        ['-w 20 -e GGC', 'example/mev.faa'],
        ['-r', 'example/mev.faa'],
        ['-r -w 20', 'example/mev.faa'], 
        ['-r -e CCA', 'example/mev.faa'],
        ['-r -w 20 -e CCA', 'example/mev.faa'],
        ['-R -s', 'example/mev.faa'],
        ['', 'example/random_seq.faa'],
        ['-w 10', 'example/random_seq.faa'],
    ]

    # list of extra tests
    extra_tests = []

    # construct list of tests to perform
    test_suite = []
    if args.f:
        test_suite = simple_tests + extra_tests
    else:
        test_suite = simple_tests

    print_summary(args, len(test_suite))
    
    # package each binary being run with the associated paths
    test_binaries = [
        {"name": "ref",
         "bin_path": os.path.join(CDS_HOME, "test/bin/CDSfoldRef"),
         "out_path": os.path.join(CDS_HOME, 'test/output/ref_output_{}.txt'),
        },
        {"name": "cds",
         "bin_path": os.path.join(CDS_HOME, "test/bin/CDSfoldLatest"),
         "out_path": os.path.join(CDS_HOME, 'test/output/cds_output_{}.txt'),
        },
    ]

    # dictionary holding all raw results in form 
    # {"test_name_1": [run1output, run2output, ...], "test_name_2" : []}

    raw_results = {}
    for test in test_binaries:
        test_result = util.run_cdsfold_many(test_suite, test, quiet=args.q, save=args.s)
        raw_results = {**raw_results, **test_result}

    # put the results in a more useful format
    results = util.clean_results(raw_results, test_suite)
    diff_output(results, args)
