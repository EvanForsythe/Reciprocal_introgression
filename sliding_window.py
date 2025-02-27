# import modules
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import re
import sys

# Debug print of command-line arguments for troubleshooting
print('Arguments: ', str(sys.argv))

#################################
# Parse Command-Line Arguments
#################################

#Create an argument parser object
parser = argparse.ArgumentParser(description= "Sliding window analysis")

# add arguments
parser.add_argument('-j', '--job_name', type=str, required=True, help='Unique job name for this run')
parser.add_argument('-w', '--win_len', type = int, required = False, default = 10000, help = 'Window length for sliding window analysis')
parser.add_argument('-o', '--outgroup', type = str, required = False, default = 'Outgroup', help = 'Name of the outgroup')
parser.add_argument('-1', '--P1', type = str, required = False, default = 'Pop1', help = 'Name of species P1')
parser.add_argument('-2', '--P2', type = str, required = False, default = 'Pop2', help = 'Name of species P2')
parser.add_argument('-3', '--P3', type = str, required = False, default = 'Pop3', help = 'Name of species P3')

# Parse and assign the command-line arguments
args = parser.parse_args() 
job_name = args.job_name
win_len = args.win_len
outgroup = args.outgroup
P1 = args.P1
P2 = args.P2
P3 = args.P3

#########################################
# Set Up File Paths and Data Structures
#########################################

# Define output folder and input sequence file (FASTA)
output_folder = f"output_{job_name}"
seq_file = os.path.join(output_folder, f"{job_name}.fa")

print(f"Analyzing the fasta file located at: {seq_file}\n")

# Define column names for the sliding window results DataFrame
col1 = 'Window_Number'
col2 = 'Window_Start_Site'
col3 = 'Window_Stop_Site'
col4 = 'Number_of_ABBA_Sites'
col5 = 'Number_of_BABA_Sites'
col6 = 'Number_of_AABB_Sites'
col7 = 'D_Statistic'

# Initialize an empty DataFrame to store results for each window
df = pd.DataFrame(columns = [col1,col2,col3,col4,col5,col6,col7])

# Compile regular expressions for matching species names
r1 = re.compile(P1)
r2 = re.compile(P2)
r3 = re.compile(P3)
rO = re.compile(outgroup)

print("Created dataframe with the following column names:")
print(df)

####################################
# Read Sequence Data and Alignment
####################################

# Parse the FASTA file and store sequences in a dictionary
sequence_data = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))

# Get the number of sequences (rows) directly from the dictionary
num_sequences = len(sequence_data)
print(f"Number of sequences:{num_sequences}")

# Get the number of sites (columns) from the length of the first sequence
# Assuming all sequences are aligned and of equal length
num_sites = len(next(iter(sequence_data.values())).seq)
print(f"Number of sites:{num_sites}")

# Deterimne if the number fo windows that fit into the alignment based on win_len
n_windows = num_sites // args.win_len
print(f"Number of windows:{n_windows}")

# Get the list of keys (sequence names)
keys = list(sequence_data.keys())

print(f"Sequence IDs {keys}")

# Identify the sequences corresponding to each species using regex matching
P1 = next(filter(r1.match,keys))
P2 = next(filter(r2.match,keys))
P3 = next(filter(r3.match,keys))
O = next(filter(rO.match,keys))
print("P1:", P1)
print("P2:", P2)
print("P3:", P3)
print("O:", O)

################################
# Sliding Window Analysis Loop
################################

print("Starting sliding window analysis...\n")

# Loop over each window and compute counts of informative sites for D-statistic
for j in range(n_windows):

    if j % 100 == 0:
        print(f"Finshed {j} windows...")

    # Calculate start and stop positions for the current window    
    start = j * win_len
    stop = start + win_len

    # Initialize counts for different site patterns 
    ABBA = 0
    BABA = 0
    AABB = 0
    AAAA = 0 # Count for non-polymorphic sites

    # Loop over each site in the window
    for i in range(start, stop):
        # Skip sites with gaps in any sequence
        if sequence_data[P1][i] == "-" or sequence_data[P2][i] == "-"  or sequence_data[P3][i] == "-" or sequence_data[O][i] == "-":
            pass
        # Count sites where all sequences are identical (non-polymorphic) 
        elif  sequence_data[P1][i]  == sequence_data[P2][i]  == sequence_data[P3][i] == sequence_data[O][i]:
            AAAA += 1
        # Identify ABBA pattern: Outgroup = P1 and P2 = P3
        elif sequence_data[O][i] == sequence_data[P1][i]  and sequence_data[P2][i] == sequence_data[P3][i]:
            ABBA += 1
        # Identify BABA pattern: Outgroup = P2 and P1 = P3
        elif sequence_data[O][i] == sequence_data[P2][i] and sequence_data[P1][i] == sequence_data[P3][i]:
            BABA += 1
        # Identify AABB pattern: Outgroup = P3 and P1 = P3
        elif sequence_data[O][i] == sequence_data[P3][i] and sequence_data[P1][i] == sequence_data[P2][i]:
            AABB += 1
         
    # Calculate the D-statistic for the current window if informative sites exist
    if ABBA + BABA > 0:
        D = (ABBA - BABA) / (ABBA + BABA)
    else:
        D = 0
        
    # Record the window number, positions, site counts, and computed D-statistic in the DataFrame as a list
    df.loc[j] = [j + 1, start + 1, stop, ABBA, BABA, AABB, D]

# Print the column names for verification
print("Column names:", df.columns)

########################################
# Check for Sufficient Variable Sites
########################################
if df['Number_of_AABB_Sites'].median() < 10:
    print("ERROR: Detected low # of variable sites in sequences. Revise parameters used in Int_sim.py and/or sliding_window.py. Stopping...")
    sys.exit(1)


###############################################
# Save Results to CSV
###############################################

print("Saving windows df to the csv file.")

# Define output filename for sliding window results
output_file = os.path.join(output_folder, f"{job_name}_windows_with_d_stat.csv")
# Write the DataFrame to CSV without the index column
df.to_csv(output_file, index=False)
