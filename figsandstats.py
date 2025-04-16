#import modules
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Create an argument parser to accept command-line inputs
parser = argparse.ArgumentParser(description="Figs and Stats")

# Add job name argument
parser.add_argument('-j', '--job_name', type=str, required=True, help='Unique job name for this run')

# Get the arguments from the parser and store as variables
args = parser.parse_args()
job_name = args.job_name

# Define output directory and input file paths
output_folder = f"output_{job_name}"
sim_file = os.path.join(output_folder, f"{job_name}_introgression_info.csv")
win_file = os.path.join(output_folder, f"{job_name}_windows_with_d_stat.csv")

# Define the quantitative log file
quant_log_file = "Quant_results_log.tsv"

# Create the quantitative data log file with column headers if it doesn't exist
if not os.path.isfile(quant_log_file):
    print("Creating quant file...\n")
    with open(quant_log_file, "a") as f:
        f.write("Job_name\tMean_no_int\tMean_2to3\tMean_3to3\tMean_recip\tMedian_no_int\tMedian_2to3\tMedian_3to2\tMedian_recip\tMeans_test\tMedians_test\n")

# Read in CSV files as DataFrames
print(f"Reading in simulation data from {sim_file}\n")
tracts_df = pd.read_csv(sim_file)
print(f"Reading in window data from {win_file}\n")
win_df = pd.read_csv(win_file)

# Explicitly cast to integers to avoid any possible mismatch issues (float and int) later.
tracts_df['Start_Site'] = tracts_df['Start_Site'].astype(int)
tracts_df['Stop_Site'] = tracts_df['Stop_Site'].astype(int)
win_df['Window_Start_Site'] = win_df['Window_Start_Site'].astype(int)
win_df['Window_Stop_Site'] = win_df['Window_Stop_Site'].astype(int)

# Sort the introgression dataframe by start site
tracts_df_sorted = tracts_df.sort_values(by='Start_Site')
current_site = 0

# Collect new rows in a list to concatenate all at once later
new_rows = []

print(f"Adding no-introgression tracts between introgressed regions...\n")

# Add "No Introgression" tracts where gaps exist between introgressed regions
for index, row in tracts_df_sorted.iterrows():
    # If there is a gap between the current site and the next introgressed tract
    if row['Start_Site'] > current_site:
        # Create a new "No Introgression" row that fills the gap
        new_row = pd.DataFrame({'Introgression_Type': ['No_Int'], 
                                'Start_Site': [current_site], 
                                'Stop_Site': [row['Start_Site'] - 1]})
        new_rows.append(new_row)
    
    # Update the current site counter to be **past** the current tract's stop position
    if row['Stop_Site'] > current_site:
        current_site = row['Stop_Site'] + 1

# Ensure the last segment of the genome is also marked if no introgression occurs at the end
last_row = win_df.iloc[-1]
if last_row['Window_Stop_Site'] > current_site:
    new_row = pd.DataFrame({'Introgression_Type': ['No_Int'],
                            'Start_Site': [int(current_site)],
                            'Stop_Site': [int(last_row['Window_Stop_Site'])]})
    new_rows.append(new_row)

# Concatenate all new rows at once
tracts_df_sorted = pd.concat([tracts_df_sorted] + new_rows, ignore_index=True)

# Re-sort the dataframe to ensure tracts are in order
tracts_df_re_sorted = tracts_df_sorted.sort_values(by='Start_Site')

# Initialize all windows explicitly as No_Int
win_df['Introgression_Type'] = 'NA'

# Initialize a list to store window subsets and D-stat averages
Intwindows = []


print(f"assigning introgression type to windows...\n")

print(tracts_df_re_sorted.head())

# Loop through the windows and pull assign introgression type

for index in win_df.index:
    temp_row=win_df.iloc[index]
    temp_start=temp_row["Window_Start_Site"]
    temp_stop=temp_row["Window_Stop_Site"]
    temp_midpoint=(temp_start+temp_stop)//2

    # Find the row where Start_Site <= number and Stop_Site > number
    matching_rows = tracts_df_re_sorted[
        (tracts_df_re_sorted['Start_Site'] <= temp_midpoint) &
        (tracts_df_re_sorted['Stop_Site'] > temp_midpoint)
    ]

    # Determine the output based on how many matches we found
    if len(matching_rows) == 1:
        temp_int_type = matching_rows.iloc[0]['Introgression_Type']
    elif len(matching_rows) > 1:
        # Prioritize 'Recip' if it exists
        if 'Recip' in matching_rows['Introgression_Type'].values:
            temp_int_type = 'Recip'
        else:
            temp_int_type="Unknown"
            print(f"Warning: couldn't find correct tract for window (index:{index})")
    else:
        temp_int_type = None  # Or raise an error if no matching interval found

    win_df.loc[index, 'Introgression_Type'] = temp_int_type

# Save the updated windows file
win_df.to_csv(os.path.join(output_folder, f"{job_name}_windows_with_int_info.csv"), index=False)

##########################
### Plot Introgressed Tracts ###
##########################

print(f"Creating tracts figure...\n")

#Create color palette for use in both figures
color_palette = {"No_Int": "C3", "pop2 to pop3": "C1", "pop3 to pop2": "C0", "Recip": "C2"}

# Get the start and stop sites of each type
migrating_pop3_to_pop2 = tracts_df_re_sorted[tracts_df_re_sorted['Introgression_Type'] == 'pop3 to pop2'][['Start_Site', 'Stop_Site']].to_numpy()
migrating_pop2_to_pop3 = tracts_df_re_sorted[tracts_df_re_sorted['Introgression_Type'] == 'pop2 to pop3'][['Start_Site', 'Stop_Site']].to_numpy()
recip_introgression = tracts_df_re_sorted[tracts_df_re_sorted['Introgression_Type'] == 'Recip'][['Start_Site', 'Stop_Site']].to_numpy()
no_introgression = tracts_df_re_sorted[tracts_df_re_sorted['Introgression_Type'] == 'No_Int'][['Start_Site', 'Stop_Site']].to_numpy()

# Create a figure with fixed size
fig = plt.figure(figsize=(10.0,5.0))

# Set y-axis limit to control spacing of plotted elements
plt.ylim(-1,6)

### plot the no introgression tracts (red)
plt.hlines(
	[4] * len(no_introgression), no_introgression[:,0], no_introgression[:,1], color = "C3", lw=10, label ="No Introgression" )

### plot the introgressed tracts pop3 -> pop2 (blue)
plt.hlines(
	[3] * len(migrating_pop3_to_pop2), migrating_pop3_to_pop2[:,0], migrating_pop3_to_pop2[:,1], color="C0", lw=10, label="pop3 -> pop2 introgression") 

### plot the introgressed tracts pop2 -> pop3 (orange)
plt.hlines(
	[2] * len(migrating_pop2_to_pop3), migrating_pop2_to_pop3[:,0], migrating_pop2_to_pop3[:,1], color="C1", lw=10, label="pop2 -> pop3 introgression")


### plot the recip introgressed tracts (green)
plt.hlines(
	[1] * len(recip_introgression), recip_introgression[:,0], recip_introgression[:,1], color="C2", lw=10, label="Reciprocal introgression")
	
#plt.plot(valid_windows_df['Window_Start_Site'], valid_windows_df['D_Statistic'], '-', color = 'blue', linewidth = 1)     
# Plot the continuous line in a neutral color (e.g., gray)
plt.plot(win_df['Window_Start_Site'], win_df['D_Statistic'], '-', color = 'gray', linewidth = 1)     

# Plot points colored by Introgression_Type
for intro_type, color in color_palette.items():
    # Filter the dataframe for the current introgression type
    subset_df = win_df[win_df['Introgression_Type'] == intro_type]
    
    # Plot the points for this subset
    plt.scatter(subset_df['Window_Start_Site'], 
                subset_df['D_Statistic'], 
                color=color, 
                label=intro_type, 
                s=3)  # s controls the size of the points

plt.axhline(y=0, color='r', linestyle='-')

#Format plot
plt.title(f"Introgressed tracks")
plt.xlabel("Genomic position")
plt.text(-0.025, 0.57, 'Tracts', transform=plt.gca().transAxes, rotation=90, va='center')
plt.text(-0.025, 0.15, 'D-Statistic', transform=plt.gca().transAxes, rotation=90, va='center')

plt.yticks([])

plt.legend()

fig_file = os.path.join(output_folder, f"{job_name}_figsandstats.pdf")
fig.savefig(fig_file)

plt.close()
plt.rcdefaults()


print(f"Creating violin plot figure...\n")

# Plotting

category_order = ["No_Int", "pop2 to pop3", "pop3 to pop2", "Recip"]


sns.violinplot(x="Introgression_Type", y="D_Statistic", data=win_df, hue="Introgression_Type",
               order=category_order, palette=color_palette, cut=0, bw_adjust=0.5, alpha=0.7)
plt.ylim(-1, 1)
plt.axhline(y=0, color='r', linestyle='--')
violin_file = os.path.join(output_folder, f"{job_name}_violinplot.pdf")
plt.savefig(violin_file)
plt.close()

# Compute Statistics + Pass/Fail
avg_vals = win_df.groupby('Introgression_Type')['D_Statistic'].mean()
med_vals = win_df.groupby('Introgression_Type')['D_Statistic'].median()

print("Mean values:\n")
print(avg_vals)

print("\nMedian values:\n")
print(med_vals)


# Check if reciprocal conditions are met
means_test = "Passed_means" if avg_vals.get('pop3 to pop2', 0) > 0 and avg_vals.get('pop2 to pop3', 0) > 0 and avg_vals.get('Recip', 0) < 0 else "Failed_means"
medians_test = "Passed_medians" if med_vals.get('pop3 to pop2', 0) > 0 and med_vals.get('pop2 to pop3', 0) > 0 and med_vals.get('Recip', 0) < 0 else "Failed_medians"

print("Writing a line to Quant log...\n")

# Write results to quant log file
with open(quant_log_file, "a") as f:
    f.write(f"{job_name}\t{avg_vals.get('No_Int', np.nan)}\t{avg_vals.get('pop3 to pop2', np.nan)}\t{avg_vals.get('pop2 to pop3', np.nan)}\t{avg_vals.get('Recip', np.nan)}\t"
            f"{med_vals.get('No_Int', np.nan)}\t{med_vals.get('pop3 to pop2', np.nan)}\t{med_vals.get('pop2 to pop3', np.nan)}\t{med_vals.get('Recip', np.nan)}\t"
            f"{means_test}\t{medians_test}\n")

print("Done with analysis!")