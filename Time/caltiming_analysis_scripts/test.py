import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import cPickle as pickle 
import time

import DataSet
import Event

if __name__ == "__main__":

	print "Loading file...", 
	sys.stdout.flush()
	data = pickle.load(open("../pickles/pions/noB/5GeV_1000evts.p", 'rb'))
	print "Done."
	sys.stdout.flush()

	#data.timeReco(plotting = True)
	#data.events[7].energyDisplay(True)
	#data.events[100].timeDisplay(True)
	#data.events[100].printEvent()
	#data.events[10].plotTimeHist(40)
	#data.events[10].algo_linearFirstTimeByLayer()

	#ti = time.time()
	#data.smearAndSave(.001, 0.01, "../pickles/electrons/10GeV_smeared_1ps_1perc.p")
	#tf = time.time()
	#print "took " + str(tf - ti) + " seconds to smear" 
	#smearedData = pickle.load(open("../pickles/electrons/10GeV_smeared_1ps_1perc.p", 'rb'))
	#smearedData.events[10].plotTimeHist(40)
	#data.avTimeHist(1000, 0, 1)
	#d, t = data.events[9].timeVsDepth()
	#plt.plot(d, t, 'ko')
	#plt.show()

	#data.plotAllDvsT()
	"""
	for i in range(len(data.events)):
		d, t = data.events[i].timeVsDepth()
		plt.plot(d, t, 'ko')
		#plt.xlabel("Depth into cal. (mm)")
		#plt.ylabel("Time of hit (ns)")
		plt.show()
		if(raw_input(">") == 'q'):
			break
	"""

	"""
	hits, bin_edges = data.events[100].timeHist(10)
	times = []
	for i in range(0, len(bin_edges)-1):
		times.append((bin_edges[i]+bin_edges[i+1])/2.0)

	plt.bar(times, hits, width = times[1]-times[0], color = 'blue')
	plt.show()
	"""
