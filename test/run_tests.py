import os
from os import path
import subprocess
import argparse
    
CDS_HOME = os.getenv("CDS_HOME")

def run_all_tests(test_suite, args):
    """
    run all tests comparing the behavior of the latest CDSfold binary with either
    a saved gold standard output or the results of running a reference binary

    inputs: test_suite (list) - list of tests to run
            args - object with command line arguments
    """
    # package each binary being run with the associated paths
    test_binaries = [
        {"name": "ref",
         "bin_path": os.path.join(CDS_HOME, "test/bin/CDSfoldRef"),
         "out_path": os.path.join(CDS_HOME, 'test/output/ref_output_{}.txt'),
        },
        {"name": "cds",
         "bin_path": os.path.join(CDS_HOME, "test/bin/CDSfold"),
         "out_path": os.path.join(CDS_HOME, 'test/output/cds_output_{}.txt'),
        },
    ]

    # results dictionary mapping name of binary run to list of output strings,
    # one string for each test
    results = {}
    for test_bin in test_binaries:
        results[test_bin["name"]] = []

    for test_bin in test_binaries:
        if not args.q: 
            print("\n==== Executing: {} ==== ".format(test_bin["name"])) 
        
        for i, test in enumerate(test_suite):
            if not args.q: 
                print("Running test: {}".format(test))
           
            # set up the command for running the test
            if len(test[0]) == 0:
                # no command line args
                cmd = [test_bin["bin_path"], 
                       os.path.join(CDS_HOME, test[1])]
            else:
                # command line args used
                cmd = [test_bin["bin_path"],
                       test[0],
                       os.path.join(CDS_HOME, test[1])]
            
            # run command and capture output
            output = subprocess.run(cmd, capture_output=True)
            results[test_bin["name"]].append(output.stdout.decode("utf-8"))
           
            # save output to file
            if args.s:
                # file to save output in 
                save_path = test_bin["out_path"].format(i)
                
                # delete file it if already exists
                if os.path.exists(save_path): 
                    os.remove(save_path)  # clear the file being output to
       
                if not args.q:
                    print("Saving test output to {}".format(save_path))
                
                with open(save_path, "w") as f:
                    # result is the last added to list of results
                    f.writelines(["Executed command: {}\n".format(str(test)), 
                                  results[test_bin["name"]][-1]])
       
        # dump the results to file if specified by command line arguments

    return results

def clean_results(raw_results, test_suite):
    '''clean the results from run_all_tests. Converts the raw_results to
    a dictionary of the form:
    {
        test_bin: [
        {
            'cmd': "", 
            'run_time': str, 
            'lines': []
        },
        {
            test_1 info...
        }
        ]
    }'''

    results = {}

    # iterate through the tests
    for key in raw_results.keys():
        results[key] = []
        for i in range(len(test_suite)):
            test_result = {} 
            test_result["cmd"] = str(test_suite[i])
            test_result["lines"] = raw_results[key][i].split("\n")
            
            # look for line with run time 
            for line in test_result["lines"]:
                if "Running time" in line:
                    runtime_list = line.split()
                    # runtime is the element before "minutes" 
                    runtime = runtime_list[runtime_list.index("minutes") - 1]
                    test_result["Run time (s)"] = float(runtime) * 60 

            results[key].append(test_result)

    return results

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
    raw_results = run_all_tests(test_suite, args)
    results = clean_results(raw_results, test_suite)
    diff_output(results, args)
