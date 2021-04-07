import os
from os import system, path
    
CDS_HOME = os.getenv("CDS_HOME")

def run_all_tests(fn):
    """
    Conduct simple tests of the CDSfold executable, accumulating one logfile.
    Tests currently hit "regular" mode, bp width limit, codon exclusion,
    both, and 'random backtrack' mode.

    Just added tests of the 'reverse' mode, which engages in FE maximization.
    """
    cds_path = os.path.join(CDS_HOME, "src/CDSfold")
    
    tests = [                   # list of tests format [cmd line args, sequence file]
        ['', 'example/mev.faa'],
        ['-w 20', 'example/mev.faa'],
        ['-e GGU', 'example/mev.faa'],
        ['-w 20 -e GGC', 'example/mev.faa'],
        ['-r', 'example/mev.faa'],
        ['-r -w 20', 'example/mev.faa'], 
        ['-r -e CCA', 'example/mev.faa'],
        ['-r -w 20 -e CCA', 'example/mev.faa'],
        ['-R -s', 'example/mev.faa'],
    #    ['', 'example/random_seq.faa'],
    #    ['-w 10', 'example/random_seq.faa'],
    ]

    system("rm -rf {}".format(fn))  # clear the file being output to
    for test in tests:
        # run the test and append to the test file
        system(cds_path + ' {} {} >> {}'.format(test[0], os.path.join(CDS_HOME, test[1]), fn))

def diff_output():
    '''diffs the gold output and the newly generated output'''
    
    with open(os.path.join(CDS_HOME, 'test/output/gold_standard_output.txt')) as f, \
    open(os.path.join(CDS_HOME, 'test/output/new_output.txt')) as g:
        for ii, (line1, line2) in enumerate(zip(f.readlines(), g.readlines())):
            # Filter a run-time line -- currently off for development to check
            # if we improve or degrade performance.
            # if 'Runing' in line1 and 'Runing' in line2: continue
            if line1 != line2:
                print('{} --- {}'.format(ii, line1.strip()))
                print('{} +++ {}'.format(ii, line2.strip()))

if __name__ == '__main__':
    if path.exists(os.path.join(CDS_HOME, 'test/output/gold_standard_output.txt')):
        run_all_tests(os.path.join(CDS_HOME, 'test/output/new_output.txt'))
        diff_output()
    else:
        print('no regression reference exists, creating it!')
        run_all_tests(os.path.join(CDS_HOME, 'test/output/gold_standard_output.txt'))

        
