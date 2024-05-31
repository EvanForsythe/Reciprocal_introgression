# Load needed modules
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import re
import glob
import sys
import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# Read in user arguments
parser = argparse.ArgumentParser(description= "Paramscan")

parser.add_argument('-i', '--indir', type = str, required = True, help = 'Full path to the output directory (directory should contain files *runstats.csv. String should end in "/")')

args = parser.parse_args()

indir = args.indir

# Use glob.glob to get a list of all the friles in the folder that end in "*runstats.csv"
in_files=glob.glob(indir+"*runstats.csv")

#print(in_files)

# Start an empty dataframe that we can add to with each roudn of the loop
col0 = 'replicate'
col1 = 'no_int_d_avg'
col2 = 'p3p2_d_avg'
col3 = 'p2p3_d_avg'
col4 = 'recip_d_avg'
col5 = 'threshold'
df = pd.DataFrame(columns = [col0,col1,col2,col3,col4,col5])






j = 0

# Loop through all the files and read each in as a dataframe
for file in in_files:
    file_df = pd.read_csv(file, header = None)
# inside of the loop pull out the imporatant information from the file

    df.loc[j] = {col0 :os.path.splitext(file)[0], col1 : file_df.iloc[0,1], col2 :file_df.iloc[0,2], col3: file_df.iloc[0,3], col4 :file_df.iloc[0,4],col5 :file_df.iloc[0,5]}
    j += 1

#sort the dataframe
df = df.sort_values(by = 'replicate')

output_file = 'paramscandata.csv'
df.to_csv(output_file, index=False)
# replicate #column 0
# threashold #column 5
# p3p2_d # column 2
# p2p3_d # column 3
# recip_d #column 4
# no_int_d #column 1

# Add all of these numbers to the empty dataframe (one row per 'run')



# Plot the data in the dataframe
def plot_function(d_stat_column, threshold_column, dataframe, title):
    plt.figure(figsize=(10.0,5.0))
    sns.swarmplot(data=dataframe, x=threshold_column, y=d_stat_column, hue='replicate', dodge=True)
    plt.ylim(-1,1)
    plt.title(title)
    plt.xlabel('Threshold')
    plt.ylabel('D-statistic')
    #plt.legend()
    plt.savefig('paramscan.pdf')
    
plot_function('p3p2_d_avg', 'threshold', df, 'p3p2 Introgression')
