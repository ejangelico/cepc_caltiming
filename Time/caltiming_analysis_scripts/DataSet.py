import sys
import Event
import cPickle as pickle 
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import Helper

def FWHM(xVals, yVals):
		xVals = xVals.tolist()
		yVals = yVals.tolist()

		maxIndex = yVals.index(max(yVals))
		peakCount      = yVals[maxIndex]
		peakCountx  = xVals[maxIndex]

		for i in range(maxIndex, len(xVals)):
			if yVals[i] < peakCount/2.0:
				break
		x0 = xVals[i-1]
		y0 = yVals[i-1]
		x1 = xVals[i]
		y1 = yVals[i]
		tHalfUp = x0 + (0.5 *peakCount - y0) *(x1-x0)/(y1-y0)

		for i in range(0, maxIndex):
			if yVals[i] > peakCount/2.0:
				break
		x0 = xVals[i-1]
		y0 = yVals[i-1]
		x1 = xVals[i]
		y1 = yVals[i]
		tHalfDown = x0 + (0.5 *peakCount - y0) *(x1-x0)/(y1-y0)

		return tHalfUp-tHalfDown


class DataSet:
	def __init__(self, events = None):
		self.events = events
		self.tSmear = tSmear 	#smear 1sigma in ns
		#fractional energy resolution
		#takes values between [0.0, 1.0]
		self.eSmear = eSmear

		self.layerBins = None

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

	#histograms the longitudinal shower depth
	#as defined in the CALICE paper as "Z0"
	#in units of the pion interaction length in tungsten
	def histLongitudinalDepth(self):
		Z0 = []
		for ev in self.events:
			Z0.append(ev.getShowerDepth())

		fig, ax = plt.subplots()
		ax.hist(Z0)
		ax.set_xlabel("Longitudinal shower depth (lambda_i)")
		plt.show()



	def timeReco(self, algo = 0, plotting = False):
		tEstList = []
		tTruList = []
		tDiffList = []
		for event in self.events:
			if algo == 0:
				tEst, tTru = event.algo_linearFirstTimeByLayer()
			else:
				print "Please specify the time reconstruction algorithm"
				sys.exit()
			tEstList.append(tEst)
			tTruList.append(tTru)
			tDiffList.append(tEst-tTru)

		tDiffCounts, tDiffBins = np.histogram(tDiffList, 50)
		tAv = np.average(tDiffList)
		tFWHM = FWHM(tDiffBins, tDiffCounts)
		tStd = np.std(tDiffList)
		tMed = np.median(tDiffList)
		tSkewness = scipy.stats.skew(tDiffList)
		tSkewTest = scipy.stats.skewtest(tDiffList)[0] # z-score of skewness test

		print "Mean difference:   ", round(1000*tAv, 4), "ps"
		print "Median difference: ", round(1000*tMed, 4), "ps"
		print "FHWM:              ", round(1000*tFWHM, 4), "ps"
		print "Standard Deviation:", round(1000*tStd, 4), "ps"
		print "Skewness:          ", round(tSkewness, 4) 
		print "Skewtest z-score:  ", round(tSkewTest, 4)

		if plotting:
			tDiffBins = [x*1000 for x in tDiffBins]
			plt.bar(tDiffBins[:-1], tDiffCounts, width = tDiffBins[1]-tDiffBins[0], color = 'b')
			plt.xlabel("$t_{reco} - t_{true}$" + " (ps) ", fontsize = 20)
			plt.ylabel("Counts/bin for 1 GeV electrons", fontsize = 20)
			plt.show()



