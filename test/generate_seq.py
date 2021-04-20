#
# randomly generate a an amino acid sequence and save it in a .faa file
#
import sys
import os
import random
from datetime import datetime

random.seed(datetime.now())     # use the current time as random seed

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

def gen_seq(seq_length, num_sequences):
    '''randomly generate an amino acid sequence of a given length'''
        
    aa_string = ""
    for i in range(num_sequences):
        aa_string += "> Random sequence no. {}\n".format(i) 
        for i in range(seq_length):
            aa_string += amino_acids[random.randint(0, len(amino_acids) - 1)]

        aa_string += "*\n\n"

    # dump the string to a file
    with open("random_seq.faa", "w") as f:
        f.write(aa_string)

    print("Generated {} random sequences with length {}".format(num_sequences, seq_length))

if __name__ == "__main__":
    # first argument - length of each sequence
    # second argument - number of randomly generated sequences

    gen_seq(int(sys.argv[1]), int(sys.argv[2]))
