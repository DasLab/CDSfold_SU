# read in seed file and obtain a data structure mapping
# ONE code per family to a sequence. write that structure
# out

import pandas as pd
import matplotlib.pyplot as plt

objs = []


with open('Pfam-A.seed', encoding='latin-1') as f:
    obj = {}
    for line in f:
        if line[:7] == '#=GF AC':
            obj = {'Pfam_Accession': line[8:].strip()}
        elif line[:7] == '#=GF ID':
            obj['One_Word_Pfam_Descriptor'] = line[8:].strip()
        elif line[:7] == '#=GF DE':
            obj['Pfam_Description_String'] = line[8:].strip()
        elif line[0] == '#': continue
        elif line.strip() == '//': continue
        else:
            line = line.split()
            try:
                potential_seq = line[1].replace('.', '').upper().strip()
            except:
                print(line)
                quit()
                
            # potentially skip any sequences with forbidden ambiguous
            # characters here: B, X, Z, O, J, U
            obj['Uniprot_Accession_Sequence_Positions'] = line[0].strip()
            obj['Protein_Seq'] = potential_seq
    
            objs.append({k: v for k,v in obj.items()})

df = pd.DataFrame.from_dict(objs)
df['Protein_Seq_Length'] = df.apply(lambda row: len(row['Protein_Seq']), axis=1)

df.to_csv('all_seed_seq.csv')
