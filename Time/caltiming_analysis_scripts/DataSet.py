import sys
import Event
import cPickle as pickle 
import matplotlib.pyplot as plt

class DataSet:
	def __init__(self, events = None):
		self.events = events
		self.tSmear = tSmear 	#smear 1sigma in ns
		#fractional energy resolution
		#takes values between [0.0, 1.0]
		self.eSmear = eSmear

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

	def smearAndSave(self, tsm, esm, outfilename):
		#don't smear this event if it already has been
		#smeared. Then the full smear is the addition of the
		#smears in quadrature. Just start from the original file
		#if(self.tSmear != 0 or self.eSmear != 0):
		#	print "This data set has already been smeared"
		#	print "Find the original set and smear that instead"
		#	return


		#make a new list of events
		smearedEvents = []

		for ev in self.events:
			newEvent = ev.getSmearedEvent(tsm, esm)
			smearedEvents.append(newEvent)

		#create new data set
		smearedData = DataSet(smearedEvents, tsm, esm)
		#save the dataset to file
		pickle.dump(smearedData, open(outfilename, 'wb'))
		return
