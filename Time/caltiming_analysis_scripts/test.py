import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import cPickle as pickle 

import DataSet
import Event

if __name__ == "__main__":
	data = pickle.load(open("../pickles/electrons/AnaHit_electron_1GeV_2501.p", 'rb'))
	#data.events[100].energyDisplay()
	#data.events[100].timeDisplay()
	#data.events[100].printEvent()
	#data.events[101].plotTimeHist(1000)
	#data.avTimeHist(1000, 0, 1)
	
	data.plotAllDvsT()
	d, t = data.events[110].timeVsDepth()
	plt.plot(d, t, 'ko')
	plt.xlabel("Depth into cal. (cm)")
	plt.ylabel("Time of hit (ns)")
	plt.show()
	"""
	hits, bin_edges = data.events[100].timeHist(10)
	times = []
	for i in range(0, len(bin_edges)-1):
		times.append((bin_edges[i]+bin_edges[i+1])/2.0)

	plt.bar(times, hits, width = times[1]-times[0], color = 'blue')
	plt.show()
	"""
