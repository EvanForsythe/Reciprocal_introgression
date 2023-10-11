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

#Set wd
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#os.chdir('/Users/esforsythe/Documents/OSU/Work/Research/Reciprocal_introgression/')

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Script for simulating introgression')

#Add arguments
parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of this script. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-s', '--Seq_len', type=int, metavar='', required=True, help='Specify an interger to set length of total simulateed alignment (e.g. 20000000') 
parser.add_argument('-p', '--Prop_int', type=float, metavar='', required=True, help='Specify the proportion of genome to be introgressed with each introgression event (e.g. 0.2') 

#Define the parser
args = parser.parse_args()

JOBname=args.JOBname
Seq_len=args.Seq_len
Prop_int=args.Prop_int


#JOBname="MYTEST"
#Seq_len=20000000
#Prop_int=0.2


#Get list of taxa
taxa_names=["Africa", "Eurasia", "Neanderthal", "Chimpanzee"]

#set up a highly simplified demographic history of human+Neanderthal demography and simulate a single chromosome of 20Mb in length

#Assign sequence length
sequence_length=Seq_len


#Setup the simulations
time_units = 1000 / 25  # Conversion factor for kya to generations
demography = msprime.Demography()
# The same size for all populations; highly unrealistic!
Ne = 10**4
demography.add_population(name="Africa", initial_size=Ne)
demography.add_population(name="Eurasia", initial_size=Ne)
demography.add_population(name="Neanderthal", initial_size=Ne)
demography.add_population(name="Chimpanzee", initial_size=Ne)

# introgression 50 kya
demography.add_mass_migration(
    time=50 * time_units, source='Eurasia', dest='Neanderthal', proportion=Prop_int)

#opposite direction introgression
# introgression 50 kya
demography.add_mass_migration(
    time=50 * time_units, source='Neanderthal', dest='Eurasia', proportion=Prop_int)


# Eurasian 'merges' backwards in time into Africa population, 70 kya
demography.add_mass_migration(
    time=70 * time_units, source='Eurasia', dest='Africa', proportion=1)

# Neanderthal 'merges' backwards in time into African population, 300 kya
demography.add_mass_migration(
    time=300 * time_units, source='Neanderthal', dest='Africa', proportion=1)

# Africa 'merges' backwards in time into Chimp population, TBD kya
demography.add_mass_migration(
    time=500 * time_units, source='Africa', dest='Chimpanzee', proportion=1)

ts = msprime.sim_ancestry(
    recombination_rate=1e-9,
    sequence_length=sequence_length, 
    samples=[
        msprime.SampleSet(1, ploidy=1, population='Africa'),
        msprime.SampleSet(1, ploidy=1, population='Eurasia'),
        # Neanderthal sample taken 30 kya
        msprime.SampleSet(1, ploidy=1, time=30 * time_units, population='Neanderthal'),
        msprime.SampleSet(1, ploidy=1, population='Chimpanzee'),
    ],
    demography = demography,
    record_migrations=True,  # Needed for tracking segments.
)


# Generate mutations on the tree sequence
ts_mutes = msprime.sim_mutations(ts, rate=1e-8, random_seed=None)


#write a fasta file
ts_mutes.write_fasta(JOBname+".fa", reference_sequence=tskit.random_nucleotides(ts.sequence_length))


#Commands for replacing seq ids
replace_cmd0="sed -i '' 's/n0/Africa/' "+JOBname+".fa"
replace_cmd1="sed -i '' 's/n1/Eurasia/' "+JOBname+".fa"
replace_cmd2="sed -i '' 's/n2/Neanderthal/' "+JOBname+".fa"
replace_cmd3="sed -i '' 's/n3/Chimpanzee/' "+JOBname+".fa"

#Edit the fasta file
#Run all the commands (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search(JOBname, replace_cmd0): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd0, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
if re.search(JOBname, replace_cmd1): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd1, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
if re.search(JOBname, replace_cmd2): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd2, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
if re.search(JOBname, replace_cmd3): #Check if cmd contains expected string and run if so
    subprocess.call(replace_cmd3, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


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

#Get the tracts from N -> E
migrating_nead_to_euro = get_migrating_tracts(ts, "Neanderthal")

#Get the tracts from N -> E
migrating_euro_to_nean = get_migrating_tracts(ts, "Eurasia")

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

recip_introgression = find_overlap_intervals(migrating_nead_to_euro, migrating_euro_to_nean)



#Create a plot of introgression tracts


fig = plt.figure(figsize=(8.0,8.0))

### plot the introgressed tracts N -> E
plt.hlines(
    [1] * len(migrating_nead_to_euro), migrating_nead_to_euro[:,0], migrating_nead_to_euro[:,1], color="C0", lw=10, label="N -> E introgression")

### plot the introgressed tracts E -> N
plt.hlines(
    [2] * len(migrating_euro_to_nean), migrating_euro_to_nean[:,0], migrating_euro_to_nean[:,1], color="C1", lw=10, label="E -> N introgression")


### plot the introgressed tracts E -> N
plt.hlines(
    [3] * len(recip_introgression), recip_introgression[:,0], recip_introgression[:,1], color="C2", lw=10, label="Reciprocal introgression")


#Format plot
plt.title(f"Introgressed tracks")
plt.xlabel("Genomic position")
plt.ylim(0, 5)
plt.yticks([])
plt.legend()
plt.show()

fig.savefig(JOBname+".pdf")



