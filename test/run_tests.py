from os import system, path

def run_all_tests(fn):
    """
    Conduct simple tests of the CDSfold executable, accumulating one logfile.
    Tests currently hit "regular" mode, bp width limit, codon exclusion,
    both, and 'random backtrack' mode.

    Just added tests of the 'reverse' mode, which engages in FE maximization.
    """

    system('./src/CDSfold example/mev.faa > {}'.format(fn))
    system('./src/CDSfold -w 20 example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -e GGU example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -w 20 -e GGC example/mev.faa >> {}'.format(fn))


    system('./src/CDSfold -r example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -r -w 20 example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -r -e CCA example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -r -w 20 -e CCA example/mev.faa >> {}'.format(fn))

    system('./src/CDSfold -R -s example/mev.faa >> {}'.format(fn))

def diff_output():
    with open('test/output/gold_standard_output.txt') as f, open('test/output/new_output.txt') as g:
        for ii, (line1, line2) in enumerate(zip(f.readlines(), g.readlines())):
            # Filter a run-time line -- currently off for development to check
            # if we improve or degrade performance.
            # if 'Runing' in line1 and 'Runing' in line2: continue
            if line1 != line2:
                print('{} --- {}'.format(ii, line1.strip()))
                print('{} +++ {}'.format(ii, line2.strip()))

if __name__ == '__main__':
    if path.exists('test/output/gold_standard_output.txt'):
        run_all_tests('test/output/new_output.txt')
        diff_output()
    else:
        print('no regression reference exists, creating it!')
        run_all_tests('test/output/gold_standard_output.txt')
        