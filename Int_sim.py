#!/usr/bin/env python3

#Based on examples from: https://tskit.dev/tutorials/introgression.html

#The full manual for thet tskit API
#https://tskit.dev/tskit/docs/stable/python-api.html?highlight=fasta#tskit.TreeSequence.as_fasta

import sys
import msprime
import numpy as np
import tskit
import os
import argparse
import pandas as pd
import shutil as sh
import itertools


# Set wd
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

# Set up an argumanet parser to accept command-line inputs
parser = argparse.ArgumentParser(description='Script for simulating introgression')

# Define required and optional arguments
parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of this script. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-s', '--Seq_len', type=int, metavar='', required=False, default=10000000, help='Specify an interger to set length of total simulateed alignment (default = 10000000')
parser.add_argument('-p', '--Prop_int', type=float, metavar='', required=False, default=0.2, help='Specify the proportion of genome to be introgressed with each introgression event (default = 0.2)')
parser.add_argument('-m', '--Mut_rate', type=float, metavar='', required=False, default=0.0000005, help='Specify the mutation rate (default = 0.000001)')
parser.add_argument('-r', '--Recomb_rate', type=float, metavar='', required=False,default=0.0000000005, help='Specify the recomb rate (default =.000000001)')
parser.add_argument('-n', '--Ne', type=int, metavar='', required=False, default=10000, help='Specify the effective pop size (Ne) (default = 10000)')

# Define timing of events (in years ago) for introgression and speciation events
parser.add_argument('-1', '--t_int', type=int, metavar='', required=False, default=40000 , help='Time of introgression (years ago) (default = 40000)')
parser.add_argument('-2', '--t_sp12', type=int, metavar='', required=False,default=80000 , help='Time of first most recent speciation (years ago) (default = 80000)')
parser.add_argument('-3', '--t_sp123', type=int, metavar='', required=False,default=120000 , help='Time of second most recent speciation (default = 120000)')
parser.add_argument('-4', '--t_sp1234', type=int, metavar='', required=False, default=200000 , help='Time of third  most recent speciation (default = 200000)')



#Define the parser
args = parser.parse_args()

# Assign parsed arguments to variables
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

# Create an output folder for this simulation run
output_folder = f"output_{JOBname}"
if not os.path.exists(output_folder):
     os.makedirs(output_folder)

#################################
# SET UP DEMOGRAPHIC MODEL 
#################################

# Define taxon names for simulation
taxa_names=["Pop1", "Pop2", "Pop3", "Outgroup"]

print(f"Running simulations from the following populations {str(taxa_names)}...\n")

# Initialize the demographic model with msprime.Demography
demography = msprime.Demography()

# Loop through taxa and add each as a population to the demography
for t_name in taxa_names:
    demography.add_population(name=t_name, initial_size=Ne)
    print(f'Adding population: {t_name}')

# Define introgression events (bidirectional migration between Pop2 and Pop3)
introgression_events = []

# Add mass migration event: Pop2 -> Pop3
demography.add_mass_migration(time=t_int, source="Pop2", dest="Pop3", proportion=Prop_int)
introgression_events.append({"time": t_int, "source": "Pop2", "dest": "Pop3", "proportion": Prop_int})

# Add mass migration event: Pop3 -> Pop2
demography.add_mass_migration(time=t_int, source="Pop3", dest="Pop2", proportion=Prop_int)
introgression_events.append({"time": t_int, "source": "Pop3", "dest": "Pop2", "proportion": Prop_int})

# Define speciation events by adding mass migrations to merge populations
demography.add_mass_migration(time=t_sp12, source="Pop2", dest="Pop1", proportion=1)
demography.add_mass_migration(time=t_sp123, source="Pop3", dest="Pop1", proportion=1)
demography.add_mass_migration(time=t_sp1234, source="Pop1", dest="Outgroup", proportion=1)

####################################
# SIMULATE GENEOLGY AND MUTATIONS 
####################################

print(f"Simulating geneology and mutations...\n")

# Simulate tree sequence with recombination using msprime
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


# Generate mutations along the simulated tree sequence
ts_mutes = msprime.sim_mutations(ts, rate=Mut_rate, random_seed=None)

#################################
# COMPUTE SIMULATION STATISTICS
#################################

# Count the number of migration tracts (haplotype blocks) and calculate their average length
total_length = 0
num_tracts = 0

for migration in ts.migrations():
    tract_length = migration.right - migration.left  # Calculate tract length
    total_length += tract_length
    num_tracts += 1

# Calculate the average length
average_length = total_length / num_tracts if num_tracts > 0 else 0

print(f"{num_tracts} total tracts simulated with an average length of {average_length}.\n")

# Function to calculate total sequence divergence (average pairwise divergence)
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
print(f"Total sequence divergence (average pairwise divergence): {total_divergence:.6f} substitutions/site\n")

######################################################
# Count Introgressed Mutations and Tract Statistics
######################################################

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
            if start <= current_end:  # Overlapping tract; extend current tract
                current_end = max(current_end, end)
            else:  # No overlap
                merged_tracts.append((current_start, current_end))
                current_start, current_end = start, end
        merged_tracts.append((current_start, current_end))  # Append the last tract

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


print(f"Tallying statistics about simulations...\n")

# Calculate statistics related to introgressed regions
num_int_tracts, average_int_tract_length, num_diagnostic_mutations, total_introgressed_length = count_introgressed_mutations(ts_mutes, introgression_events)

print(f"Number of introgressed tracts: {num_int_tracts}")
print(f"Average length of int tract: {average_int_tract_length}")
print(f"Number of 'diagnostic' mutations in introgressed tracts: {num_diagnostic_mutations}")
print(f"Total length of non-overlapping introgressed tracts: {total_introgressed_length}\n")

# Write simulation statistics to a quantitative log file
quant_log_file = "Sim_stats_log.tsv"
if not os.path.isfile(quant_log_file):
     print("Creating Sim_stats_log.tsv...\n")
     with open(quant_log_file, "a") as f:
          f.write("JOBname\tSeq_len\tProp_int\tMut_rate\tRecomb_rate\tnum_tracts\taverage_length\ttotal_divergence\tnum_int_tracts\taverage_int_tract_length\ttotal_introgressed_length\tnum_diagnostic_mutations\n")

print("Adding to Sim_stats_log.tsv...\n")
with open (quant_log_file, "a") as f:
     f.write(f"{JOBname}\t{Seq_len}\t{Prop_int}\t{Mut_rate}\t{Recomb_rate}\t{num_tracts}\t{average_length}\t{total_divergence}\t{num_int_tracts}\t{average_int_tract_length}\t{total_introgressed_length}\t{num_diagnostic_mutations}\n")


##################################
# Write Simulated Data to FASTA
##################################

print("Writing simulation data to fasta file...\n")

# Define the filename for the FASTA output
fasta_filename = os.path.join(output_folder, f"{JOBname}.fa")
# Write the simulated sequence with mutations to a FASTA file
ts_mutes.write_fasta(fasta_filename, reference_sequence=tskit.random_nucleotides(ts.sequence_length))

#######################################
# Edit FASTA Headers for Readability
#######################################

print("Editing seq IDs in fasta file...\n")

# Define file names for reading and writing
fasta_read_filename = fasta_filename
fasta_write_filename = fasta_filename + ".tmp"

# Open the original FASTA file for reading and a temporary file for writing the updated headers
with open(fasta_read_filename, 'r') as fasta_read_handle, open(fasta_write_filename, "a") as fasta_write_handle:
    for line in fasta_read_handle:
        if line.startswith(">"):
            # Clean and update the sequence ID
            id_temp = line.strip() # Removes "\n"
            id_clean = id_temp.replace(">", "") # Removes ">" by replacing with nothing.
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


########################################################
# Track Migration Tracts and Export Introgression Info
########################################################
def get_migrating_tracts(ts, dest_pop):
    """
    Extracts migration tracts for a given destination population
    """
    dest_id = [p.id for p in ts.populations() if p.metadata['name']==dest_pop][0]
    migrating_tracts = []
    for migration in ts.migrations():
        if migration.dest == dest_id:
            migrating_tracts.append((int(migration.left), int(migration.right)))
    return np.array(migrating_tracts) 

print("Finding overlapping tracts (i.e. reciprocal introgression tracts)\n")

def find_overlap_intervals(arr1, arr2):
    """
    Finds overlapping intervals between two sets of tracts
    Returns the overlapping regions as a numpy array
    """
    overlap_intervals = []
    for interval1 in arr1:
        for interval2 in arr2:
            # Check if the intervals overlap
            if interval1[1] >= interval2[0] and interval1[0] <= interval2[1]:
                # Calculate the overlapping region
                overlap_start = max(interval1[0], interval2[0])
                overlap_stop = min(interval1[1], interval2[1])
                overlap_intervals.append([overlap_start, overlap_stop])
    
    # Remove duplicate intervals by converting to a DataFrame and dropping duplicates
    overlap_intervals_df = pd.DataFrame(overlap_intervals)
    overlap_intervals_df_nodup = overlap_intervals_df.drop_duplicates()
    return np.array(overlap_intervals_df_nodup)

# Get the migration tracts from Pop3 -> Pop2
migrating_pop3_to_pop2 = get_migrating_tracts(ts, "Pop2")

# Get the migration tracts from Pop2 -> Pop3
migrating_pop2_to_pop3 = get_migrating_tracts(ts, "Pop3")

# Get the reciprocal introgression (overlapping) tracts
recip_introgression = find_overlap_intervals(migrating_pop3_to_pop2, migrating_pop2_to_pop3)

# Verify that all types of introgression are present; exit if any type is missing
if len(migrating_pop3_to_pop2) == 0 or len(migrating_pop2_to_pop3) == 0 or len(recip_introgression) == 0:
    print("Error: Not all types of introgression (pop3 to pop2, pop2 to pop3, reciprocal) exist. Revise parameters. Quitting...")
    sys.exit()

print("Creating a dataframe with simulation info...\n")

# Write introgression tract information to CSV
output_file= os.path.join(output_folder, f"{JOBname}_introgression_info.csv")
df = pd.DataFrame(columns = ["Introgression_Type", "Start_Site", "Stop_Site"])

# Add tracts for each introgression type to the DataFrame
for migration in migrating_pop3_to_pop2:
	 df.loc[len(df)] = {"Introgression_Type": "pop3 to pop2", "Start_Site": migration[0], "Stop_Site": migration[1]}

for migration in migrating_pop2_to_pop3:
	df.loc[len(df)] = {"Introgression_Type": "pop2 to pop3", "Start_Site": migration[0], "Stop_Site": migration[1]}

for migration in recip_introgression:
	df.loc[len(df)] = {"Introgression_Type": "Recip", "Start_Site": migration[0], "Stop_Site": migration[1]}

print("writing the dataframe as a csv file...\n")

df.to_csv(output_file, index=False)

print("Done with Int_sim run!")
