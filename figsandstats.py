#import modules
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

# Create an argument parser to accept command-line inputs
parser = argparse.ArgumentParser(description= "Figs and Stats")

# Add job name argument
parser.add_argument('-j', '--job_name', type = str, required = True, help = 'Unique job name for this run')

# Get the arguments from the parser and store as variables
args = parser.parse_args()
job_name = args.job_name

# Define output directory and input file paths
output_folder = f"output_{job_name}"
sim_file = os.path.join(output_folder, f"{job_name}_introgression_info.csv")
win_file = os.path.join(output_folder, f"{job_name}_windows_with_d_stat.csv")


# Define the quantitative log file
quant_log_file = "Quant_results_log.tsv"

try:
	# Create the quantitative data log file with column headers if it doesn't exist
	if not os.path.isfile(quant_log_file):
		with open(quant_log_file, "a") as f:
			f.write("Job_name\tMean_no_int\tMean_2to3\tMean_3to2\tMean_recip\tMedian_no_int\tMedian_2to3\tMedian_3to2\tMedian_recip\tMeans_test\tMedians_test\n")


	# Read in CSV files as dataframe
	sim_df = pd.read_csv(sim_file)
	win_df = pd.read_csv(win_file)

	# Sort the introgression dataframe by start site
	sim_df = sim_df.sort_values(by = 'Start_Site')

	# Create a counter to track genomic positions
	current_site = 0

	# Add "No Introgression" tracts where gaps exist between introgressed regions
	for index, row in sim_df.iterrows():
		# If there is a gap between the current site and the next introgressed tract
		if row['Start_Site'] > current_site:
			# Create a new "No Introgression" row that fills the gap
			new_row = pd.DataFrame({'Introgression_Type': ['No_Int'], # Marks the tract as non-introgressed
			 						'Start_Site': [int(current_site)], # The gap starts from the current position
			  						'Stop_Site': [int(row['Start_Site']) -1]}) # Ends just before the next tract
					
			# Append the new "No Introgression" row to the dataframe
			sim_df = pd.concat([sim_df, new_row], ignore_index = True)

		# Update the current site counter to be **past** the current tract's stop position
		if row['Stop_Site'] > current_site:
			current_site = row['Stop_Site'] + 1
		

	# Ensure the last segment of the genome is also marked if no introgression occurs at the end

	last_row = win_df.iloc[-1] # Get the last row of the window dataframe

	# If there is a remaining portion at the end that isn't covered by introgressed tracts
	if last_row['Window_Stop_Site'] > current_site:
			# Create a final "No Introgression" row to fill this last portion
			new_row = pd.DataFrame({'Introgression_Type': ['No_Int'],
									 'Start_Site': [int(current_site)],  # Starts from the last tracked position
									  'Stop_Site': [int(last_row['Window_Stop_Site'])]})  # Ends at the last window's stop site

			# Append this final "No Introgression" tract to the dataframe
			sim_df = pd.concat([sim_df, new_row], ignore_index = True)

	# Re-sort the dataframe to ensure tracts are in order
	sim_df = sim_df.sort_values(by = 'Start_Site')
	
	# Create a new dataframe to store windows with introgression information
	new_windows_df=pd.DataFrame(columns=win_df.columns.tolist()+['Introgression_Type'])

	# Create lists to store computed values
	Davglist = []   # Stores average D-statistics for each introgression block
	Intwindows = []  # Stores window numbers associated with introgressed tracts

	# Assign introgression type to corresponding windows
	for ind, row in sim_df.iterrows():
		temp_int_type = row['Introgression_Type']
		temp_block_start = row['Start_Site']
		temp_block_stop = row['Stop_Site']

		# Select windows that fall within introgression blocks
		filteredtemp = win_df[(win_df['Window_Start_Site'] >= temp_block_start) & 
							(win_df['Window_Stop_Site'] <= temp_block_stop)]

		# Skip empty results
		if filteredtemp.empty:
			print(f"No matching windows for tract: {row}")
			Davglist.append(np.nan)
			continue

		# Create a copy to avoid SettingWithCopyWarning
		filteredtemp = filteredtemp.copy()
		# Assign introgression type to the filtered windows
		filteredtemp['Introgression_Type'] = temp_int_type

		# Append rows to the new dataframe
		new_windows_df = pd.concat([new_windows_df, filteredtemp], ignore_index=True)
		
		# Computer average D-statistic for the windows iwthin the introgression tract
		Dstatslist = filteredtemp[['D_Statistic']].to_numpy()

		
		# Check if there are items in the list and take average if so
		if len(Dstatslist) < 1:
			Daverage = np.nan
		else:
			Daverage = np.mean(Dstatslist)
		
		# Add D average values to list
		Davglist.append(Daverage)
		
		# Track window numbers associated with introgression
		for i in list(filteredtemp['Window_Number']):
			if i not in Intwindows:
				Intwindows.append(i)


	# Save the updated windows file
	new_windows_df.to_csv(os.path.join(output_folder, f"{job_name}_windows_with_int_info.csv"), index=False)

	# Add the average D-stat values for each tract to the dataframe
	sim_df['Average_Dstat_for_windows_in_tract'] = Davglist
	
	# Save the updated introgression file
	output_file = os.path.join(output_folder, f"{job_name}_figsandstats.csv") 
	sim_df.to_csv(output_file, index=False)



	##########################
	### Computer Statistics ###
	##########################

	# Calculate the mean and median D-statistics for each introgression type
	avg_noint = np.nanmean(sim_df.loc[sim_df['Introgression_Type'] == 'No_Int']['Average_Dstat_for_windows_in_tract'])
	avg_32 = np.nanmean(sim_df.loc[sim_df['Introgression_Type'] == 'pop3 to pop2']['Average_Dstat_for_windows_in_tract'])
	avg_23 = np.nanmean(sim_df.loc[sim_df['Introgression_Type'] == 'pop2 to pop3']['Average_Dstat_for_windows_in_tract'])
	avg_recip = np.nanmean(sim_df.loc[sim_df['Introgression_Type'] == 'Recip']['Average_Dstat_for_windows_in_tract'])

	med_noint = np.nanmedian(sim_df.loc[sim_df['Introgression_Type'] == 'No_Int']['Average_Dstat_for_windows_in_tract'])
	med_32 = np.nanmedian(sim_df.loc[sim_df['Introgression_Type'] == 'pop3 to pop2']['Average_Dstat_for_windows_in_tract'])
	med_23 = np.nanmedian(sim_df.loc[sim_df['Introgression_Type'] == 'pop2 to pop3']['Average_Dstat_for_windows_in_tract'])
	med_recip = np.nanmedian(sim_df.loc[sim_df['Introgression_Type'] == 'Recip']['Average_Dstat_for_windows_in_tract'])


	

	# Check if reciprocal conditions are met
	means_test = "Passed_means" if avg_32 > 0 and avg_23 > 0 and avg_recip < 0 else "Failed_means"
	medians_test = "Passed_medians" if med_32 > 0 and med_23 > 0 and med_recip < 0 else "Failed_medians"


	# Write results to median and mean log files
	with open (quant_log_file, "a") as f:
		f.write(f"{job_name}\t{avg_noint}\t{avg_32}\t{avg_23}\t{avg_recip}\t{med_noint}\t{med_32}\t{med_23}\t{med_recip}\t{means_test}\t{medians_test}\n")



	##########################
	### Plot Introgressed Tracts ###
	##########################

	# Get the start and stop sites of each type
	migrating_pop3_to_pop2 = sim_df[sim_df['Introgression_Type'] == 'pop3 to pop2'][['Start_Site', 'Stop_Site']].to_numpy()
	migrating_pop2_to_pop3 = sim_df[sim_df['Introgression_Type'] == 'pop2 to pop3'][['Start_Site', 'Stop_Site']].to_numpy()
	recip_introgression = sim_df[sim_df['Introgression_Type'] == 'Recip'][['Start_Site', 'Stop_Site']].to_numpy()
	no_introgression = sim_df[sim_df['Introgression_Type'] == 'No_Int'][['Start_Site', 'Stop_Site']].to_numpy()

	# Create a new column to store the middle point of a window
	win_df['Average_Site'] = (win_df['Window_Start_Site'] + win_df['Window_Stop_Site']) / 2



	# Create a figure with fixed size
	fig = plt.figure(figsize=(10.0,5.0))
	
	# Set y-axis limit to control spacing of plotted elements
	plt.ylim(-1,6)

	### plot the no introgression tracts (red)
	plt.hlines(
		[4] * len(no_introgression), no_introgression[:,0], no_introgression[:,1], color = "C3", lw=10, label =" No Introgression" )

	### plot the introgressed tracts pop3 -> pop2 (blue)
	plt.hlines(
		[3] * len(migrating_pop3_to_pop2), migrating_pop3_to_pop2[:,0], migrating_pop3_to_pop2[:,1], color="C0", lw=10, label="pop3 -> pop2 introgression") 

	### plot the introgressed tracts pop2 -> pop3 (orange)
	plt.hlines(
		[2] * len(migrating_pop2_to_pop3), migrating_pop2_to_pop3[:,0], migrating_pop2_to_pop3[:,1], color="C1", lw=10, label="pop2 -> pop3 introgression")


	### plot the recip introgressed tracts (green)
	plt.hlines(
		[1] * len(recip_introgression), recip_introgression[:,0], recip_introgression[:,1], color="C2", lw=10, label="Reciprocal introgression")
		

	plt.plot(win_df['Window_Start_Site'], win_df['D_Statistic'], '-', color = 'blue', linewidth = 1)     
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

	
	###################
	### Violin Plot ###
	###################

	valid_windows_df = new_windows_df.dropna(subset=['Introgression_Type'])

	# Generate the plot
	sns.violinplot(x = "Introgression_Type", y = "D_Statistic", data = valid_windows_df, split = False)
	plt.ylim(-1,1)
	plt.axhline(y=0, color='r', linestyle = '--')

	#Create a file handle
	violin_file = os.path.join(output_folder, f"{job_name}_violinplot.pdf")
	#Save the figure
	plt.savefig(violin_file)
	plt.close()


except Exception as e:
	#Log NA values in case of failure
	with open(quant_log_file, "a") as f:
		f.write(f"{job_name}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tFailed_means\tFailed_medians\n")

	print(f"Error occurred: {e}")
	sys.exit(1)