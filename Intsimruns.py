import os
import sys


run_name = 'devin_test20240515_'
threshold_increment = 0.05

for i in range(1,10):
	print('Running Iteration ' + str(i))
	command1 ='python Int_sim.py -j Intsimruns/'+ run_name + str(i) + ' -s 10000000 -p 0.20 -m 0.00001 -r 0.0000000002 -n 10000'
	command2 = 'python sliding_window.py Intsimruns/' + run_name + str(i) +'.fa -w 10000 -o Outgroup -1 Pop1 -2 Pop2 -3 Pop3'
	command3 = 'python figsandstats.py -w Intsimruns/'+ run_name + str(i) + '_slidingwindow.csv -s Intsimruns/' + run_name + str(i) + '.csv -t '
	#os.system(command1)
	#os.system(command2)
	
	for j in range(1,10):
	    print(command3 + str(threshold_increment * j))
	    sys.exit()
	    os.system(command3 + str(threshold_increment * j))
	    os.rename('Intsimruns/' + run_name + str(i) + '_figsandstats.csv', 'Intsimruns/' + run_name + str(i) + '-' + str(j) + '_figsandstats.csv')
	    os.rename('Intsimruns/' + run_name + str(i) + '_figsandstats.pdf', 'Intsimruns/' + run_name + str(i) + '-' + str(j) + '_figsandstats.pdf')
	    os.rename('Intsimruns/' + run_name + str(i) + '_fviolin.pdf', 'Intsimruns/' + run_name + str(i) + '-' + str(j) + '_fviolin.pdf')
	    os.rename('Intsimruns/' + run_name + str(i) + '_runstats.csv', 'Intsimruns/' + run_name + str(i) + '-' + str(j) + '_runstats.csv')
	#print(command1)
	#print(command2)
	#print(command3)


## j = 1, threshold = 0.05
## j = 2, threshold = 0.10
## j = 3, threshold = 0.15

