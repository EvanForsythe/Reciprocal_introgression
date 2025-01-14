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
import itertools


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
parser.add_argument('-m', '--Mut_rate', type=float, metavar='', required=False, default=0.000001, help='Specify the mutation rate (default = 0.000001)')
parser.add_argument('-r', '--Recomb_rate', type=float, metavar='', required=False,default=0.000000001, help='Specify the recomb rate (default =.000000001)')
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

# Create output folder
output_folder = f"output_{JOBname}"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#Pop1=Africa
#Pop2=Eurasia
#POp3=Neanderthal
#Pop4=Chimanzee

#Get list of taxa
taxa_names=["Pop1", "Pop2", "Pop3", "Outgroup"]

# Track introgression events manually
introgression_events = []

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
#Track info about introgression event
introgression_events.append({"time": t_int, "source": "Pop2", "dest": "Pop3", "proportion": Prop_int})

#opposite direction introgression
# introgression 50 kya
demography.add_mass_migration(
    time=t_int, source="Pop3", dest="Pop2", proportion=Prop_int)
#Track info about introgression event
introgression_events.append({"time": t_int, "source": "Pop1", "dest": "Ghost", "proportion": Prop_int})


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


# Simulation statistics
# Count the number of tracts (haplotype blocks) and calculate their average length
total_length = 0
num_tracts = 0

for migration in ts.migrations():
    tract_length = migration.right - migration.left  # Calculate tract length
    total_length += tract_length
    num_tracts += 1

# Calculate the average length
average_length = total_length / num_tracts if num_tracts > 0 else 0

# Calculate total sequence divergence (average pairwise divergence)
def calculate_total_divergence(tree_sequence):
    divergence = 0
    pair_count = 0

    # Iterate over pairs of samples
    for i, j in itertools.combinations(range(tree_sequence.num_samples), 2):
        divergence += tree_sequence.divergence([[i], [j]])  # Pairwise divergence
        pair_count += 1

    average_divergence = divergence / pair_count if pair_count > 0 else 0
    return average_divergence


total_divergence = calculate_total_divergence(ts_mutes)
print(f"Total sequence divergence (average pairwise divergence): {total_divergence:.6f} substitutions/site")

#Get the 'diagnostic SNPs' that indicate introgression
def count_introgressed_mutations(ts, introgression_events):
    """
    Counts the number of mutations that occurred along branches in introgressed regions,
    calculates the average length of introgressed tracts, and computes the total
    length of non-overlapping introgressed tracts.

    Args:
        ts (TreeSequence): The simulated tree sequence.
        introgression_events (list): List of tracked introgression events.

    Returns:
        tuple: A tuple containing:
            - int: Total number of introgressed tracts.
            - float: Average length of introgressed tracts.
            - int: Total number of "diagnostic" mutations in introgressed tracts.
            - float: Total length of non-overlapping introgressed tracts.
    """
    # Generate introgressed tracts based on introgression events
    introgressed_tracts = []
    for migration in ts.migrations():
        for event in introgression_events:
            source_id = next((pop.id for pop in ts.populations() if pop.metadata.get("name") == event["source"]), None)
            dest_id = next((pop.id for pop in ts.populations() if pop.metadata.get("name") == event["dest"]), None)
            if source_id is not None and dest_id is not None:
                if (event["time"] - 1 <= migration.time <= event["time"] + 1 and
                        migration.source == source_id and migration.dest == dest_id):
                    introgressed_tracts.append((migration.left, migration.right))

    # Sort and merge overlapping tracts
    introgressed_tracts.sort()
    merged_tracts = []
    if introgressed_tracts:
        current_start, current_end = introgressed_tracts[0]
        for start, end in introgressed_tracts[1:]:
            if start <= current_end:  # Overlap
                current_end = max(current_end, end)
            else:  # No overlap
                merged_tracts.append((current_start, current_end))
                current_start, current_end = start, end
        merged_tracts.append((current_start, current_end))  # Add the last tract

    # Calculate the total length of introgressed tracts
    total_tract_length = sum(end - start for start, end in merged_tracts)

    # Calculate the average tract length
    num_tracts = len(merged_tracts)
    average_tract_length = total_tract_length / num_tracts if num_tracts > 0 else 0

    # Loop through mutations and check if they are in introgressed regions
    diagnostic_mutations = 0
    for mutation in ts.mutations():
        site_position = ts.site(mutation.site).position
        for left, right in merged_tracts:
            if left <= site_position < right:
                diagnostic_mutations += 1
                break

    return num_tracts, average_tract_length, diagnostic_mutations, total_tract_length


#Call the function
num_int_tracts, average_int_tract_length, num_diagnostic_mutations, total_introgressed_length = count_introgressed_mutations(ts_mutes, introgression_events)

print(f"Number of introgressed tracts: {num_int_tracts}")
print(f"Average length of int tract: {average_int_tract_length}")
print(f"Number of 'diagnostic' mutations in introgressed tracts: {num_diagnostic_mutations}")
print(f"Total length of non-overlapping introgressed tracts: {total_introgressed_length}")

# write to a quant file
quant_log_file = "Sim_stats_log.tsv"

#Create the quantitative data log file
if not os.path.isfile(quant_log_file):
	with open(quant_log_file, "a") as f:
		f.write("JOBname\tSeq_len\tProp_int\tMut_rate\tRecomb_rate\tnum_tracts\taverage_length\ttotal_divergence\tnum_int_tracts\taverage_int_tract_length\ttotal_introgressed_length\tnum_diagnostic_mutations\n")

#Write to the quant file
with open (quant_log_file, "a") as f:
	f.write(f"{JOBname}\t{Seq_len}\t{Prop_int}\t{Mut_rate}\t{Recomb_rate}\t{num_tracts}\t{average_length}\t{total_divergence}\t{num_int_tracts}\t{average_int_tract_length}\t{total_introgressed_length}\t{num_diagnostic_mutations}\n")


## END simulation statistics



#write a fasta file
fasta_filename = os.path.join(output_folder, f"{JOBname}.fa")
ts_mutes.write_fasta(fasta_filename, reference_sequence=tskit.random_nucleotides(ts.sequence_length))

### Pseudo code for editing file:
#Create file handle for the fasta file that we wrote (open for reading)

#Create file handle for a new file (open for 'appending')

#Loop through and read each line of original file
#fasta_read_filename = JOBname+".fa"
#fasta_write_filename = JOBname+".fa.tmp"
#fasta_read_handle = open(fasta_read_filename, "r")
#fasta_write_handle = open(fasta_write_filename, "a")
#Create an empty dictionary
seq_dict = {}

#Loop through the line in the file
fasta_read_filename = fasta_filename
fasta_write_filename = fasta_filename + ".tmp"
with open(fasta_read_filename, 'r') as fasta_read_handle, open(fasta_write_filename, "a") as fasta_write_handle:
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

#print("hello")
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
    
    #create dataframe
    overlap_intervals_df = pd.DataFrame(overlap_intervals)
    
    #remove any duplicate rows from the dataframe
    overlap_intervals_df_nodup = overlap_intervals_df.drop_duplicates()

    #return the array
    return np.array(overlap_intervals_df_nodup)

recip_introgression = find_overlap_intervals(migrating_pop3_to_pop2, migrating_pop2_to_pop3)

#Check for the presence of all three introgression types

if len(migrating_pop3_to_pop2) == 0 or len(migrating_pop2_to_pop3) == 0 or len(recip_introgression) == 0:
    print("Error: Not all types of introgression (pop3 to pop2, pop2 to pop3, reciprocal) exist. Revise parameters.")
    sys.exit()
#Output to CSV
col1 = "Introgression_Type"
col2 = "Start_Site"
col3 = "Stop_Site"

csv_filename = os.path.join(output_folder, JOBname + ".csv")
df = pd.DataFrame(columns = [col1,col2,col3])

for migration in migrating_pop3_to_pop2:
	 df.loc[len(df)] = {col1: "pop3 to pop2", col2: migration[0], col3: migration[1]}

for migration in migrating_pop2_to_pop3:
	df.loc[len(df)] = {col1: "pop2 to pop3", col2: migration[0], col3: migration[1]}

for migration in recip_introgression:
	df.loc[len(df)] = {col1: "Recip", col2: migration[0], col3: migration[1]}

output_file = JOBname + ".csv"
df.to_csv(csv_filename, index=False)
