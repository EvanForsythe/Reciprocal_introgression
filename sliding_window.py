# import modules
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import re
import glob
import sys

print('Arguments: ', str(sys.argv))
#Create an argument parser object
parser = argparse.ArgumentParser(description= "Sliding window analysis")
#add arguments
parser.add_argument('seq_file', help = 'Sequence file name')
parser.add_argument('-w', '--win_len', type = int, required = True, help = 'Window length')
parser.add_argument('-o', '--outgroup', type = str, required = True, help = 'Name of the outgroup')
parser.add_argument('-1', '--P1', type = str, required = True, help = 'Name of species P1')
parser.add_argument('-2', '--P2', type = str, required = True, help = 'Name of species P2')
parser.add_argument('-3', '--P3', type = str, required = True, help = 'Name of species P3')
#Define the parser
args = parser.parse_args() 

#Store arguments
seq_file = args.seq_file
win_len = args.win_len
outgroup = args.outgroup
P1 = args.P1
P2 = args.P2
P3 = args.P3
#set base path
sequence_base_folder_path = seq_file
sequence_files = glob.glob(sequence_base_folder_path)

df = pd.DataFrame(columns = ['Window number', 'Window start site', 'Window stop site', 'Number of ABBA sites', 'Number of BABA sites', 'D-statistic for window'])

r1 = re.compile(P1)
r2 = re.compile(P2)
r3 = re.compile(P3)
rO = re.compile(outgroup)
for f in sequence_files:
        print(f)
        #set the full path to the file to read in
        sequence_file_path = f

        #Read the sequence file, storing it as a variable
        sequence_data = SeqIO.to_dict(SeqIO.parse(sequence_file_path, 'fasta'))

        keys = list(sequence_data.keys())

        #Print the number of sequences (rows) in the alignment
        num_sequences = len(keys)
        print("Number of sequences:", num_sequences)

        #Print the number of sites (columns) in the alignment
        alignment = AlignIO.read(sequence_file_path, 'fasta')
        num_sites = alignment.get_alignment_length()
        print("Number of sites:", num_sites)

        #Calculate how many windows fit within the full sequence
        n_windows = num_sites // args.win_len
        print("Number of windows:", n_windows)

        P1 = next(filter(r1.match,keys))
        P2 = next(filter(r2.match,keys))
        P3 = next(filter(r3.match,keys))
        O = next(filter(rO.match,keys))

        print("P1:", P1)
        print("P2:", P2)
        print("P3:", P3)
        print("O:", O)
        #Create a loop to iterate through the sites in the alignment
        for j in range(n_windows):
                start = j * win_len
                stop = start + win_len 
                ABBA = 0
                BABA = 0
                for i in range(start, stop):
                        if sequence_data[P1][i] == "-" or sequence_data[P2][i] == "-"  or sequence_data[P3][i] == "-" or sequence_data[O][i] == "-":
                                pass #ignores gaps
                        elif  sequence_data[P1][i]  == sequence_data[P2][i]  == sequence_data[P3][i] == sequence_data[O][i]:
                                pass #Non-polymorphic
                        elif sequence_data[O][i] == sequence_data[P1][i]  and sequence_data[P2][i] == sequence_data[P3][i]:
                                ABBA += 1
                        elif sequence_data[O][i] == sequence_data[P2][i] and sequence_data[P1][i] == sequence_data[P3][i]:
                                BABA += 1
                        
                #Calculate the D-statistic
                #print(f"D = ({ABBA} - {BABA}) / ({ABBA} + {BABA})")
                if ABBA + BABA > 0:
                        D = (ABBA - BABA) / (ABBA + BABA)
                else:
                        D = "NA"

                #CSV file
                df.loc[j] = {'Window number' :j + 1, 'Window start site' :start + 1, 'Window stop site':stop, 'Number of ABBA sites': ABBA, 'Number of BABA sites':BABA, 'D-statistic for window':D}
        output_file = f + '.csv'
        df.to_csv(output_file, index=False)
        df.dropna()
