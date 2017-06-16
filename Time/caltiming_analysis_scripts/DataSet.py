import sys
import Event
import matplotlib.pyplot as plt

class DataSet:
	def __init__(self, events = None, tSmear=0, eSmear=0):
		self.events = events

	
	# Plots the time-average histogram of all the events
	def avTimeHist(self, numBins, rangeMin, rangeMax):
		timeHist = [0 for _ in range(numBins)]
		for event in self.events:
			timeHist = [x+y for x, y in zip(timeHist, event.timeHist(numBins, rangeMin, rangeMax)[0])]
		
		binCenters = self.events[0].timeHist(numBins)[1]	

		plt.bar(binCenters, timeHist, width = binCenters[1]-binCenters[0], color = 'blue')
		plt.xlabel("Time (nanoseconds)")
		plt.ylabel("Total shower energy/bin (GeV/ns)")
		plt.show()

	def plotAllDvsT(self):
		d = []
		t = []
		for event in self.events:
			(d0, t0) = event.timeVsDepth()
			d += d0
			t += t0

		plt.plot(d, t, 'k.')
		plt.xlabel("Depth into Cal. (mm)")
		plt.ylabel("Time of Hit (ns)")
		plt.show()


