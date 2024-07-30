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

#Create a parser
parser = argparse.ArgumentParser(description= "Figs and Stats")

#Add to the parser
parser.add_argument('-s', '--sim_file', type = str, required = True, help = 'Full path to the Int_sim output csv')
parser.add_argument('-w', '--win_file', type = str, required = True, help = 'Full path to the Sliding_window  output csv')

#Get the arguments from the parser and store as variables
args = parser.parse_args()
sim_file = args.sim_file
win_file = args.win_file


#Read in CSV files as dataframe
sim_df = pd.read_csv(sim_file)
win_df = pd.read_csv(win_file)

# Create blank lists
Davglist = []
Intwindows = []


## Add "no int" tracts to the sim df
num_rows_sim_df_before = len(sim_df)
num_rows_win_df_before = len(win_df)
print(f"simdf before: {num_rows_sim_df_before}")
print(f"windf before: {num_rows_win_df_before}")

#Sort the whole Data Frame by Start Site 	
sim_df = sim_df.sort_values(by = 'Start_Site')

#Create a counter
current_site = 0

#Loop through the rows of the sim dataframe
for index, row in sim_df.iterrows():
    
    #print current site for debugging
    print(f"Current site: {current_site}, Start_Site: {row['Start_Site']}")
	#ask if the tract start site is greater than the counter
    if row['Start_Site'] > current_site:
    	# create a new dataframe with one row
        new_row = pd.DataFrame({'Introgression_Type': ['No_Int'], 'Start_Site': [current_site], 'Stop_Site': [row['Start_Site'] -1]})
                
        # Concatenate  the new row onto the sim dataframe
        sim_df = pd.concat([sim_df, new_row], ignore_index = True)

    
    #Advance the counter so that it is one greater than stop site
    if row['Stop_Site'] > current_site:
        current_site = row['Stop_Site'] + 1
    

#Get the last row of the df
last_row = win_df.iloc[-1]

#Add a "no int" tract to the very end 
if last_row['Window_Stop_Site'] > current_site:
        new_row = pd.DataFrame({'Introgression_Type': ['No_Int'], 'Start_Site': [current_site], 'Stop_Site': [last_row['Window_Stop_Site']]})
        sim_df = pd.concat([sim_df, new_row], ignore_index = True)

##Sort again to make sure the df is ordered
sim_df = sim_df.sort_values(by = 'Start_Site')


### ASK WHETHER SIMDF HAS MORE ROWS THAN WINDOWS DF??
num_rows_sim_df = len(sim_df)
num_rows_win_df = len(win_df)

print(f"sim df: {num_rows_sim_df}")
print(f"win_df: {num_rows_win_df}")

if num_rows_sim_df > num_rows_win_df:
    print("ERROR: # of tracts exceeded # of windows, which could result in unexpected behavior. Revise parameters used in Int_sim.py and/or sliding_window.py. Stopping...")
    sys.exit()

#Loop through rows in dataframe
for ind, row in sim_df.iterrows():
    #Get the start and stop site
	temp_block_start = sim_df['Start_Site'][ind]
	temp_block_stop = sim_df['Stop_Site'][ind]

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
    
    #Get all windows falling within blocks
	filteredtemp = win_df[(win_df['Window_Start_Site'] >= temp_block_start) & (win_df['Window_Stop_Site'] <= temp_block_stop)]
	
	#print(filteredtemp)
	
	#Get Dstat values for windows within the block
	Dstatslist = filteredtemp[['D_Statistic']].to_numpy()

	
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
	
	#Check if there are items in the list and take average if so
	if len(Dstatslist) < 1:
	    Daverage = 0
	else:
	    Daverage = np.mean(Dstatslist)
    
    #Add D average values to list
	Davglist.append(Daverage)
	
	
	for i in list(filteredtemp['Window_Number']):
	    if i not in Intwindows:
	        Intwindows.append(i)
	
#Add Int windows to list
#print("Intwindows:")
#print(Intwindows)

NoIntWindows=[]

for i in list(win_df['Window_Number']):
    if i not in Intwindows:
        NoIntWindows.append(i)

#print("NoIntWindows:")
#print(NoIntWindows)


#Add the Dstat averages for the windows in each tract
sim_df['Average_Dstat_for_windows_in_tract'] = Davglist
#print("THIS IS THE SIM_DF:")
#print(sim_df)




#Added 'No Introgression Tract'

#no_introgression_df = win_df[(win_df['D_Statistic'].abs() < threshold)]
#average_d_stat_no_introgression = no_introgression_df['D_Statistic'].mean()

#print(sim_df.tail())

output_file = os.path.splitext(sim_file)[0] + '_figsandstats.csv'

#output_file = 'figs_and_stats_output_file.csv' 
sim_df.to_csv(output_file, index=False)


#Get the start and stop sites of each type
migrating_pop3_to_pop2 = sim_df[sim_df['Introgression_Type'] == 'pop3 to pop2'][['Start_Site', 'Stop_Site']].to_numpy()
migrating_pop2_to_pop3 = sim_df[sim_df['Introgression_Type'] == 'pop2 to pop3'][['Start_Site', 'Stop_Site']].to_numpy()
recip_introgression = sim_df[sim_df['Introgression_Type'] == 'Recip'][['Start_Site', 'Stop_Site']].to_numpy()
no_introgression = sim_df[sim_df['Introgression_Type'] == 'No_Int'][['Start_Site', 'Stop_Site']].to_numpy()


#filtered_df = sim_df[sim_df['Introgression_Type'] == 'pop3 to pop2']
#migrating_pop3_to_pop2 = filtered_df[['Start_Site', 'Stop_Site']].to_numpy()

#filtered_df = sim_df[sim_df['Introgression_Type'] == 'pop2 to pop3']
#migrating_pop2_to_pop3 = filtered_df[['Start_Site', 'Stop_Site']].to_numpy()

#filtered_df = sim_df[sim_df['Introgression_Type'] == 'Recip']
#recip_introgression = filtered_df[['Start_Site', 'Stop_Site']].to_numpy()

#print(no_introgression_df)
#no_introgression = no_introgression_df[['Window_Start_Site', 'Window_Stop_Site']].to_numpy()


#Create a new column to store the middle point of a window
win_df['Average_Site'] = (win_df['Window_Start_Site'] + win_df['Window_Stop_Site']) / 2


#Create a plot of introgression tracts
#Locations of vert lines

fig = plt.figure(figsize=(10.0,5.0))

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

plt.plot(win_df['Window_Start_Site'], win_df['D_Statistic'], '-', color = 'blue', linewidth = 1)     
plt.axhline(y=0, color='r', linestyle='-')

#Format plot
plt.title(f"Introgressed tracks")
plt.xlabel("Genomic position")
plt.text(-0.025, 0.57, 'Tracts', transform=plt.gca().transAxes, rotation=90, va='center')
plt.text(-0.025, 0.15, 'D-Statistic', transform=plt.gca().transAxes, rotation=90, va='center')

plt.yticks([])
'''
xmin, xmax = plt.xlim()
quarter_points = [xmin + (xmax - xmin) * i / 4 for i in range(1, 4)]
for x in quarter_points:
    plt.axvline(x=x, color='b', linestyle='-')
'''
plt.legend()
#plt.show()

fig_file = os.path.splitext(sim_file)[0] + '_figsandstats.pdf'

fig.savefig(fig_file)

plt.close()
plt.rcdefaults()

###################
### Violin Plot ###
###################
# Print the df for testing
#print("THIS IS SIM_DF")
#print(sim_df)

#Generate the plot
sns.violinplot(x = "Introgression_Type", y = "Average_Dstat_for_windows_in_tract", data = sim_df, split = False)
# Set ylimits
plt.ylim(-1,1)

plt.axhline(y=0, color='r', linestyle = '--')

#Create a file handle
violin_file = os.path.splitext(sim_file)[0] + '_fviolin.pdf'
#Sive the figure
plt.savefig(violin_file)


##########################
### Save summary stats ###
##########################

# Create a file handle for a csv file (create in append mode, so that each 'run' of script create a new line in the file)
stats_handle = open(os.path.splitext(sim_file)[0] + '_runstats.csv', 'a')
# Calculate summary stats by doing the following

#print(sim_df.loc[sim_df[‘Introgression_Type’] == ‘No_Int’][“Average_Dstat_for_windows_in_tract”])
#print(sim_df[‘Introgression_Type’ == ‘No_Int’])
#print(sim_df.loc[sim_df['Introgression_Type'] == 'foo'])

avg_noint=np.average(list(sim_df.loc[sim_df['Introgression_Type'] == 'No_Int']['Average_Dstat_for_windows_in_tract']))
avg_32=np.average(list(sim_df.loc[sim_df['Introgression_Type'] == 'pop3 to pop2']['Average_Dstat_for_windows_in_tract']))
avg_23=np.average(list(sim_df.loc[sim_df['Introgression_Type'] == 'pop2 to pop3']['Average_Dstat_for_windows_in_tract']))
avg_recip=np.average(list(sim_df.loc[sim_df['Introgression_Type'] == 'Recip']['Average_Dstat_for_windows_in_tract']))

#Add a line to the file
stats_handle.write(','.join([os.path.splitext(sim_file)[0], str(avg_noint), str(avg_32), str(avg_23), str(avg_recip)]))


# Slice sim_df to get a df that only contains no_int values
# Get the average D-stat from those windows and save that value as a variable
# Repeat for the other three types of int
# After you the four different averages, append them into a list
# Write that list to the csv file (which will add a row to the csv file)
stats_handle.close()
