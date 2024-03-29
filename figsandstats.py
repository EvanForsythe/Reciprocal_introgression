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
parser = argparse.ArgumentParser(description= "Figs and Stats")

parser.add_argument('-s', '--sim_file', type = str, required = True, help = 'Full path to the Int_sim output csv')
parser.add_argument('-w', '--win_file', type = str, required = True, help = 'Full path to the Sliding_window  output csv')

args = parser.parse_args()
sim_file = args.sim_file
win_file = args.win_file


#Read in CSV files as dataframe
sim_df = pd.read_csv(sim_file)
print(sim_df)
win_df = pd.read_csv(win_file)
print(win_df)

Davglist = []

#Loop through rows in dataframe
for ind, row in sim_df.iterrows():

	temp_block_start = sim_df['Start Site'][ind]
	temp_block_stop = sim_df['Stop Site'][ind]

	#n=0

	#while temp_block_start >= win_df['Window_Start_Site'][n]:
		#temp_win_start = win_df['Window_Start_Site'][n]
		#start_index = n
		#n += 1

	#n = 0
	#while temp_block_stop > win_df['Window_Stop_Site'][n]:
		#n += 1
		#end_index = n
		#temp_win_stop = win_df['Window_Stop_Site'][n]

	#print(temp_block_start)
	#print(temp_block_stop)

	filteredtemp = win_df[(win_df['Window_Start_Site'] >= temp_block_start) & (win_df['Window_Stop_Site'] <= temp_block_stop)]
	
	Dstatslist = filteredtemp[['D-Statistic']].to_numpy()
	
	#Dstatslist = []
	#for x in range(start_index,(end_index + 1)):
		#Dstatslist.append(win_df['D-Statistic'][x])

	
	#print(temp_block_start)
	#print(temp_block_stop)
	#print(temp_win_start)
	#print(temp_win_stop)
	#print(start_index)
	#print(end_index)
	#print(Daverage)
	#print(Dstatslist)
	if not len(Dstatslist):
		Daverage = 0
	else:
		Daverage = np.mean(Dstatslist)

	Davglist.append(Daverage)
sim_df['Average_Dstat_for_windows_in_tract'] = Davglist
	
sim_df = sim_df.sort_values(by = 'Start Site')

output_file = os.path.splitext(sim_file)[0] + '_figsandstats.csv'

#output_file = 'figs_and_stats_output_file.csv' 
sim_df.to_csv(output_file, index=False)


filtered_df = sim_df[sim_df['Introgression Type'] == 'N to E']
migrating_nead_to_euro = filtered_df[['Start Site', 'Stop Site']].to_numpy()

filtered_df = sim_df[sim_df['Introgression Type'] == 'E to N']
migrating_euro_to_nean = filtered_df[['Start Site', 'Stop Site']].to_numpy()

filtered_df = sim_df[sim_df['Introgression Type'] == 'Recip']
recip_introgression = filtered_df[['Start Site', 'Stop Site']].to_numpy()

win_df['Average_Site'] = (win_df['Window_Start_Site'] + win_df['Window_Stop_Site']) / 2


#Create a plot of introgression tracts
#Locations of vert lines




fig = plt.figure(figsize=(10.0,5.0))

plt.ylim(-1,6)

### plot the introgressed tracts N -> E
plt.hlines(
    [2] * len(migrating_nead_to_euro), migrating_nead_to_euro[:,0], migrating_nead_to_euro[:,1], color="C0", lw=10, label="N -> E introgression")

### plot the introgressed tracts E -> N
plt.hlines(
    [3] * len(migrating_euro_to_nean), migrating_euro_to_nean[:,0], migrating_euro_to_nean[:,1], color="C1", lw=10, label="E -> N introgression")


### plot the recip introgressed tracts
plt.hlines(
    [4] * len(recip_introgression), recip_introgression[:,0], recip_introgression[:,1], color="C2", lw=10, label="Reciprocal introgression")

#plt.axvline(x=first_quart, color='b', linestyle='-')
#plt.axvline(x=halfway, color='b', linestyle='-')
#plt.axvline(x=third_quart, color='b', linestyle='-')

#Locations of vert lines
       
        
       # plot=df.plot(title = 'D-stat plot', x = 'Window_Number', y = 'D-Statistic', ylim = (-1.1,1.1), figsize=(10,5))
       # plt.xlabel('Window Number')
       # plt.ylabel('D-Statistic')
       # plt.axhline(y=0, color='r', linestyle='-')
       # plt.axvline(x=first_quart, color='b', linestyle='-')
       # plt.axvline(x=halfway, color='b', linestyle='-')
       # plt.axvline(x=third_quart, color='b', linestyle='-')
        #plt.savefig('savedfig1.png')

#ax2 = plt.twinx()
#ax2.set_ylim(-1,1)
#plt.plot(sim_df['Start Site'], sim_df['Average_Dstat_for_windows_in_tract'], '-', color = 'blue', linewidth = 1)

plt.plot(win_df['Window_Start_Site'], win_df['D-Statistic'], '-', color = 'blue', linewidth = 1)     
plt.axhline(y=0, color='r', linestyle='-')

#Format plot
plt.title(f"Introgressed tracks")
plt.xlabel("Genomic position")
plt.text(-0.025, 0.57, 'Tracts', transform=plt.gca().transAxes, rotation=90, va='center')
plt.text(-0.025, 0.15, 'D-Statistic', transform=plt.gca().transAxes, rotation=90, va='center')

plt.yticks([])
xmin, xmax = plt.xlim()
quarter_points = [xmin + (xmax - xmin) * i / 4 for i in range(1, 4)]
for x in quarter_points:
    plt.axvline(x=x, color='b', linestyle='-')
plt.legend()
#plt.show()
fig_file = os.path.splitext(sim_file)[0] + '_figsandstats.pdf'

fig.savefig(fig_file)

plt.close()
plt.rcdefaults()
#Violin Plot
sns.violinplot(x = "Introgression Type", y = "Average_Dstat_for_windows_in_tract", data = sim_df, split = True)
plt.ylim(-1,1)

violin_file = os.path.splitext(sim_file)[0] + '_violin.pdf'

plt.savefig(violin_file)
