from os import system, path

def run_all_tests(fn):
    #with open(fn, 'w') as f:
    system('./src/CDSfold example/mev.faa > {}'.format(fn))
    system('./src/CDSfold -w 20 example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -e GGU example/mev.faa >> {}'.format(fn))
    system('./src/CDSfold -w 20 -e GGC example/mev.faa >> {}'.format(fn))

def diff_output():
    with open('test/output/gold_standard_output.txt') as f, open('test/output/new_output.txt') as g:
        for ii, (line1, line2) in enumerate(zip(f.readlines(), g.readlines())):
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
        