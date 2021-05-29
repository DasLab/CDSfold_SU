How to work with Pfam data
--------------------------

The big Pfam alignment file is a flat text file with a bunch
of headings. I transformed it to a dataframe, selecting one
sequence per family. A script for processing that alignment
file is contained in this folder, `consume_pfam_seed.py`.

In theory, then, you could obtain the corresponding mRNA 
sequence via a df.apply invocation, like:

```
def cds_fold(row):
    # Add a * to represent the stop codon.
    s = row['Protein_Seq'] + '*' 

    filename = '{}.{}.fasta'.format(
        row['Pfam_Accession']
        row['Uniprot_Accession_Sequence_Positions'].replace('/', '__')
    )
    with open(filename, 'w') as f:
        f.write('>foo\n')
        f.write('{}\n'.format(s))
    return subprocess.check_output(['CDSfoldLatest','filename']).decode().splitlines()[-4])

df = pd.read_csv('one_seq_with_mrna.cv')

# Default temperature is 37; pass -C 60 for 'Hot Vienna'
# conditions
df['Opt_mRNA_CDSfold_Vienna2_37'] = df.apply(lambda row: cds_fold(row), axis=1)
```

but this relies on a pretty stable compute environment where 
each process can get a significant memory allocation, especially
if you choose to use `parallel_apply` or `dask` or something.

Similarly, we wanted to record the runtime as well -- a different
line in the output -- so it would be convenient if we simply
retained logfile output. So instead we used conventional job
distribution on a Slurm cluster and `cds_fold` looked for a 
logfile and, if it couldn't find one, `sbatch`ed a job if one
hadn't been queued. These details are pretty compute-environment
specific! 

Anyway, there are some roadblocks you are going to run into with
this dataset. I mentioned that I selected one sequence per Pfam.
In a few cases, there was only one sequence that I could have
chosen -- in a couple, there were none at all! The reason is that
Pfam uses ambiguous amino acids in some of their seed alignment
sequences, and I did not want to bother to work out the best way
to handle them. These characters are not particularly complicated;
essentially I/L, D/N, E/Q are sometimes interchangeable, for example.
But in terms of how to handle generating a fair data set, it wasn't
clear it was my place to give some Pfams two (or twenty!) very
similar examples (and it would complicate train/test splits, if
that's your thing). Of course, I could have just picked one
representative ssequence per Pfam -- declared J means either I or
L , but not both, for example -- and pretended this issue didn't
exist, but then where would we be. So I just skipped sequences with
that property.


We have included a python script `consume_pfam_seed_total.py` that
produces a csv with *every* sequence from the seed alignment, rather
than just one per Pfam with no ambiguous characters. I have not 
versioned the CSV because it's ~270MB, but this should work on
the Pfam-A.seed available from

```
http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz
```

and, if you're feeling flush, it'll work with the equivalent 
`Pfam-A.full.gz` as well, and with small alterations to the other
reduced sets released on the website. Transform to your own taste!
