# Run everything in the "toy peptide benchmark"
import pandas as pd
from os import system

df = pd.read_csv('toy_peptide_benchmark.csv')
protein_seqs = df["Protein Sequence"].unique()
print("Will evaluate {} sequences.".format(len(protein_seqs)))

for seq in protein_seqs:
    print(seq)
    with open('tmp.faa', 'w') as ofasta:
        ofasta.write('>test\n{}\n'.format(seq))
    system('./src/CDSfold tmp.faa')

system('rm tmp.faa')

