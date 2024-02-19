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
parser.add_argument('-s', '--Seq_len', type=int, metavar='', required=True, help='Specify an interger to set length of total simulateed alignment (e.g. 20000000') 
parser.add_argument('-p', '--Prop_int', type=float, metavar='', required=True, help='Specify the proportion of genome to be introgressed with each introgression event (e.g. 0.2)')
parser.add_argument('-m', '--Mut_rate', type=float, metavar='', required=True, help='Specify the mutation rate (e.g. 0.0000001)')
parser.add_argument('-r', '--Recomb_rate', type=float, metavar='', required=True, help='Specify the recomb rate (e.g. 0.00000001)')
parser.add_argument('-n', '--Ne', type=int, metavar='', required=True, help='Specify the effective pop size (Ne) (e.g. 10000)')


#Define the parser
args = parser.parse_args()

JOBname=args.JOBname
Seq_len=args.Seq_len
Prop_int=args.Prop_int
Mut_rate=args.Mut_rate
Recomb_rate=args.Recomb_rate
Ne=args.Ne

'''
#Arguments for testing
JOBname="reconfig"
Seq_len=2000000
Prop_int=0.2
Mut_rate=0.0000001
Recomb_rate=0.00000001
Ne=10000
'''
#Pop1=Africa
#Pop2=Eurasia
#POp3=Neanderthal
#Pop4=Chimanzee

#Get list of taxa
taxa_names=["Pop1", "Pop2", "Pop3", "Outgroup"]

#set up a highly simplified demographic history of human+Neanderthal demography and simulate a single chromosome of 20Mb in length

#Assign sequence length
#REMOVE####sequence_length=Seq_len


#Setup the simulations
time_units = 1000 / 25  # Conversion factor for kya to generations
demography = msprime.Demography()


#Loop through taxa and add each as a population
for t_name in taxa_names:
    demography.add_population(name=t_name, initial_size=Ne)
    print(f'Adding population: {t_name}')


# introgression 50 kya
demography.add_mass_migration(
    time=1000 * time_units, source="Pop2", dest="Pop3", proportion=Prop_int)

#opposite direction introgression
# introgression 50 kya
demography.add_mass_migration(
    time=1000 * time_units, source="Pop3", dest="Pop2", proportion=Prop_int)

# Speciation event
demography.add_mass_migration(
    time=2000 * time_units, source="Pop2", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=3000 * time_units, source="Pop3", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=5000 * time_units, source="Pop1", dest="Outgroup", proportion=1)



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



'''
#Commands for replacing seq ids
replace_cmd0="sed -i '' 's/n0/Africa/' "+JOBname+".fa"
replace_cmd1="sed -i '' 's/n1/Eurasia/' "+JOBname+".fa"
replace_cmd2="sed -i '' 's/n2/Neanderthal/' "+JOBname+".fa"
replace_cmd3="sed -i '' 's/n3/Chimpanzee/' "+JOBname+".fa"

#Edit the fasta file
#Run all the commands (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search(JOBname, replace_cmd0): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd0, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
else:
	print(f"unable to run the following find/replace command: {replace_cmd0}")

if re.search(JOBname, replace_cmd1): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd1, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
else:
	print(f"unable to run the following find/replace command: {replace_cmd1}")

if re.search(JOBname, replace_cmd2): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd2, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
else:
	print(f"unable to run the following find/replace command: {replace_cmd2}")

if re.search(JOBname, replace_cmd3): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd3, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
else:
	print(f"unable to run the following find/replace command: {replace_cmd3}")
'''

## Track the tracts that underwent migration

def get_migrating_tracts(ts, dest_pop):
    dest_id = [p.id for p in ts.populations() if p.metadata['name']==dest_pop][0]
    migrating_tracts = []
    # Get all tracts that migrated into the neanderthal population
    for migration in ts.migrations():
        #print(migration.dest)
        if migration.dest == dest_id:
            migrating_tracts.append((migration.left, migration.right))
    return np.array(migrating_tracts) 

#Get the tracts from Pop3 -> Pop2
migrating_pop3_to_pop2 = get_migrating_tracts(ts, "Pop3")

#Get the tracts from Pop2 -> Pop3
migrating_pop2_to_pop3 = get_migrating_tracts(ts, "Pop2")

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
    
    return np.array(overlap_intervals)

recip_introgression = find_overlap_intervals(migrating_pop3_to_pop2, migrating_pop2_to_pop3)


#Create a plot of introgression tracts
#Locations of vert lines
first_quart = math.floor(Seq_len*0.25)
halfway = math.floor(Seq_len*0.5)
third_quart = math.floor(Seq_len*0.75)



fig = plt.figure(figsize=(10.0,5.0))

### plot the introgressed tracts pop3 -> pop2
plt.hlines(
    [1] * len(migrating_pop3_to_pop2), migrating_pop3_to_pop2[:,0], migrating_pop3_to_pop2[:,1], color="C0", lw=10, label="pop3 -> pop2 introgression")

### plot the introgressed tracts E -> N
plt.hlines(
    [2] * len(migrating_pop2_to_pop3), migrating_pop2_to_pop3[:,0], migrating_pop2_to_pop3[:,1], color="C1", lw=10, label="pop2 -> pop3 introgression")


### plot the introgressed tracts E -> N
plt.hlines(
    [3] * len(recip_introgression), recip_introgression[:,0], recip_introgression[:,1], color="C2", lw=10, label="Reciprocal introgression")

plt.axvline(x=first_quart, color='b', linestyle='-')
plt.axvline(x=halfway, color='b', linestyle='-')
plt.axvline(x=third_quart, color='b', linestyle='-')


#Format plot
plt.title(f"Introgressed tracks")
plt.xlabel("Genomic position")
plt.ylim(0, 4)
plt.yticks([])
plt.legend()
#plt.show()

fig.savefig(JOBname+".pdf")

#Output to CSV
col1 = "Introgression Type"
col2 = "Start Site"
col3 = "Stop Site"
df = pd.DataFrame(columns = [col1,col2,col3])

for migration in migrating_pop3_to_pop2:
	 df.loc[len(df)] = {col1: "pop3 to pop2", col2: migration[0], col3: migration[1]}

for migration in migrating_pop2_to_pop3:
	df.loc[len(df)] = {col1: "pop2 to pop3", col2: migration[0], col3: migration[1]}

for migration in recip_introgression:
	df.loc[len(df)] = {col1: "Recip", col2: migration[0], col3: migration[1]}

output_file = JOBname + ".csv"
df.to_csv(output_file, index=False)
