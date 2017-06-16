import sys
import Event
import matplotlib.pyplot as plt

class DataSet:
	def __init__(self, events = None):
		self.events = events

	def avTimeHist(self, numBins):
		timeHist = [0 for _ in range(numBins)]
		for event in self.events:
			timeHist = [x+y for x, y in zip(timeHist, event.timeHist(numBins)[0])]
		
		bin_edges = self.events[0].timeHist(numBins)[1]	
		times = []
		for i in range(0, len(bin_edges)-1):
			times.append((bin_edges[i]+bin_edges[i+1])/2.0)

		plt.bar(times, timeHist, width = times[1]-times[0], color = 'blue')
		plt.show()
		
