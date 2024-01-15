import os


for i in range(2,51):
	print('Running Iteration ' + str(i))
	command1 ='python Int_sim.py -j Intsimruns/tenmb'+ str(i) + ' -s 10000000 -p 0.20 -m 0.00001 -r 0.0000000002 -n 10000'
	command2 = 'python sliding_window.py Intsimruns/tenmb' + str(i) +'.fa -w 10000 -o n0 -1 n1 -2 n2 -3 n3'
	command3 = 'python figsandstats.py -w Intsimruns/tenmb'+ str(i) + '_slidingwindow.csv -s Intsimruns/tenmb' + str(i) + '.csv' 
	os.system(command1)
	os.system(command2)
	os.system(command3)
	#print(command1)
	#print(command2)
	#print(command3)
