#!/usr/bin/env python3

#Based on examples from: https://tskit.dev/tutorials/introgression.html

#The full manual for thet tskit API
#https://tskit.dev/tskit/docs/stable/python-api.html?highlight=fasta#tskit.TreeSequence.as_fasta

import sys
import random
import collections
import subprocess
import matplotlib.pyplot as plt
import msprime
import numpy as np
import dataclasses
import tskit
import os
import re
import argparse
import math
import pandas as pd
import shutil as sh

#Set wd
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#os.chdir('/Users/esforsythe/Documents/OSU/Work/Research/Reciprocal_introgression/Reciprocal_introgression/')

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Script for simulating introgression')

#Add arguments
parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of this script. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-s', '--Seq_len', type=int, metavar='', required=False, default=10000000, help='Specify an interger to set length of total simulateed alignment (default = 10000000')
parser.add_argument('-p', '--Prop_int', type=float, metavar='', required=False, default=0.2, help='Specify the proportion of genome to be introgressed with each introgression event (default = 0.2)')
parser.add_argument('-m', '--Mut_rate', type=float, metavar='', required=False, default=0.0000001, help='Specify the mutation rate (default = 0.0000001)')
parser.add_argument('-r', '--Recomb_rate', type=float, metavar='', required=False,default=0.0000000001, help='Specify the recomb rate (default = 0.0000000001)')
parser.add_argument('-n', '--Ne', type=int, metavar='', required=False, default=10000, help='Specify the effective pop size (Ne) (default = 10000)')

#Time arguments
parser.add_argument('-1', '--t_int', type=int, metavar='', required=False, default=40000 , help='Time of introgression (years ago) (default = 40000)')
parser.add_argument('-2', '--t_sp12', type=int, metavar='', required=False,default=80000 , help='Time of first most recent speciation (years ago) (default = 80000)')
parser.add_argument('-3', '--t_sp123', type=int, metavar='', required=False,default=120000 , help='Time of second most recent speciation (default = 120000)')
parser.add_argument('-4', '--t_sp1234', type=int, metavar='', required=False, default=200000 , help='Time of third  most recent speciation (default = 200000)')



#Define the parser
args = parser.parse_args()

JOBname=args.JOBname
Seq_len=args.Seq_len
Prop_int=args.Prop_int
Mut_rate=args.Mut_rate
Recomb_rate=args.Recomb_rate
Ne=args.Ne
t_int=args.t_int
t_sp12=args.t_sp12
t_sp123=args.t_sp123
t_sp1234=args.t_sp1234



#Pop1=Africa
#Pop2=Eurasia
#POp3=Neanderthal
#Pop4=Chimanzee

#Get list of taxa
taxa_names=["Pop1", "Pop2", "Pop3", "Outgroup"]

#set up a demographic history 

#Setup the simulations
#time_units = 1000 / 25  # Conversion factor for kya to generations
demography = msprime.Demography()


#Loop through taxa and add each as a population
for t_name in taxa_names:
    demography.add_population(name=t_name, initial_size=Ne)
    print(f'Adding population: {t_name}')


#t_int=40000
#t_sp12=80000
#t_sp123=120000
#t_sp1234=200000

# introgression 50 kya
demography.add_mass_migration(
    time=t_int, source="Pop2", dest="Pop3", proportion=Prop_int)

#opposite direction introgression
# introgression 50 kya
demography.add_mass_migration(
    time=t_int, source="Pop3", dest="Pop2", proportion=Prop_int)

# Speciation event
demography.add_mass_migration(
    time=t_sp12, source="Pop2", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=t_sp123, source="Pop3", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=t_sp1234, source="Pop1", dest="Outgroup", proportion=1)



ts = msprime.sim_ancestry(
    recombination_rate=Recomb_rate,
    sequence_length=Seq_len, 
    samples=[
        msprime.SampleSet(1, ploidy=1, population="Pop1"),
        msprime.SampleSet(1, ploidy=1, population="Pop2"),
        msprime.SampleSet(1, ploidy=1, population="Pop3"),
        msprime.SampleSet(1, ploidy=1, population="Outgroup"),
    ],
    demography = demography,
    record_migrations=True,  # Needed for tracking segments.
)


# Generate mutations on the tree sequence
ts_mutes = msprime.sim_mutations(ts, rate=Mut_rate, random_seed=None)


#write a fasta file
ts_mutes.write_fasta(JOBname+".fa", reference_sequence=tskit.random_nucleotides(ts.sequence_length))

### Pseudo code for editing file:
#Create file handle for the fasta file that we wrote (open for reading)

#Create file handle for a new file (open for 'appending')

#Loop through and read each line of original file
fasta_read_filename = JOBname+".fa"
fasta_write_filename = JOBname+".fa.tmp"
fasta_read_handle = open(fasta_read_filename, "r")
fasta_write_handle = open(fasta_write_filename, "a")
#Create an empty dictionary
seq_dict = {}

#Loop through the line in the file
for line in fasta_read_handle:
    if line.startswith(">"):
        #Use the replae method to replace (use your dictionary)
        id_temp = line.strip() #Removes "\n"
        id_clean = id_temp.replace(">", "") #Removes ">" by replacing with nothing.
        id_new = ''
        if id_clean == "n0":
                id_new = "Pop1"
        elif id_clean == "n1":
                id_new = "Pop2"
        elif id_clean == "n2":
                id_new = "Pop3"
        elif id_clean == "n3":
                id_new = "Outgroup"
        else:
                id_new = id_clean 
        fasta_write_handle.write(">" + id_new + "\n")
    else:
       fasta_write_handle.write(line)

sh.move(fasta_write_filename, fasta_read_filename)


## Track the tracts that underwent migration

def get_migrating_tracts(ts, dest_pop):
    dest_id = [p.id for p in ts.populations() if p.metadata['name']==dest_pop][0]
    #print(dest_id)
    #print(dest_pop)
    migrating_tracts = []
    # Get all tracts that migrated into the destination population
    for migration in ts.migrations():
        #print(migration.dest)
        if migration.dest == dest_id:
            migrating_tracts.append((int(migration.left), int(migration.right)))
            #print(migrating_tracts)
    return np.array(migrating_tracts) 

#Get the tracts from Pop3 -> Pop2
migrating_pop3_to_pop2 = get_migrating_tracts(ts, "Pop3")
print(migrating_pop3_to_pop2)

#print("hello")
#Get the tracts from Pop2 -> Pop3
migrating_pop2_to_pop3 = get_migrating_tracts(ts, "Pop2")
print(migrating_pop2_to_pop3)

#Get the overlap (reciprocal introgression)
def find_overlap_intervals(arr1, arr2):
    overlap_intervals = []
    
    for interval1 in arr1:
        for interval2 in arr2:
            # Check if the intervals overlap
            if interval1[1] >= interval2[0] and interval1[0] <= interval2[1]:
                # Calculate the overlapping region
                overlap_start = max(interval1[0], interval2[0])
                overlap_stop = min(interval1[1], interval2[1])
                overlap_intervals.append([overlap_start, overlap_stop])
    
    #create dataframe
    overlap_intervals_df = pd.DataFrame(overlap_intervals)
    
    #remove any duplicate rows from the dataframe
    overlap_intervals_df_nodup = overlap_intervals_df.drop_duplicates()

    #return the array
    return np.array(overlap_intervals_df_nodup)

recip_introgression = find_overlap_intervals(migrating_pop3_to_pop2, migrating_pop2_to_pop3)

#Output to CSV
col1 = "Introgression_Type"
col2 = "Start_Site"
col3 = "Stop_Site"
df = pd.DataFrame(columns = [col1,col2,col3])

for migration in migrating_pop3_to_pop2:
	 df.loc[len(df)] = {col1: "pop3 to pop2", col2: migration[0], col3: migration[1]}

for migration in migrating_pop2_to_pop3:
	df.loc[len(df)] = {col1: "pop2 to pop3", col2: migration[0], col3: migration[1]}

for migration in recip_introgression:
	df.loc[len(df)] = {col1: "Recip", col2: migration[0], col3: migration[1]}

output_file = JOBname + ".csv"
df.to_csv(output_file, index=False)
