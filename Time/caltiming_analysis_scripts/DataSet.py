import sys
import Event
import cPickle as pickle 
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import Helper

# Given x-values and y-values of a curve with one peak, computes the full width at half max, in units of x
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
	def __init__(self, events = None, tSmear = None, eSmear = None):
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

	# Plots the depth (radius) vs time of every hit point
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

	# Makes a new dataset with smeared time/energy
	def smear(self, tsm, esm):
		print "Smearing data: tsm =", tsm, "ns; esm =", esm
		#make a new list of events
		smearedEvents = []

		for ev in self.events:
			newEvent = ev.getSmearedEvent(tsm, esm)
			smearedEvents.append(newEvent)

		#create new data set
		return DataSet(smearedEvents, tsm, esm)

	def smearAndSave(self, tsm, esm, outfilename):
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



	# Given the index of the event algorithm to use, reconstructs the reco-truth times and does stats
	def timeReco(self, algo = 0, plotting = False):
		tEstList = []
		tTruList = []
		tDiffList = []
		print "Reconstructing time for each event...",
		sys.stdout.flush()
		for event in self.events:
			if algo == 0:
				tEst, tTru = event.algo_linearFirstTimeByLayer()
			else:
				print "Please specify the time reconstruction algorithm"
				sys.exit()
			tEstList.append(tEst)
			tTruList.append(tTru)
			tDiffList.append(tEst-tTru)
		print "Done."
		sys.stdout.flush()

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

		return [tAv, tMed, tFWHM, tStd, tSkewness, tSkewTest]

	def timeRecoSmearing(self, smearTimeList = None):
		if smearTimeList is None:
			smearTimeList = [0.0, 0.001, 0.01, 0.1]
		AvList = []
		MedList = []
		FWHMList = []
		StdList = []
		for deltaT in smearTimeList:
			smearedData = self.smear(tsm = deltaT, esm = 0)
			print "Time smearing:", deltaT, "ns"
			recoStats = smearedData.timeReco()
			AvList.append(1000*recoStats[0])
			MedList.append(1000*recoStats[1])
			FWHMList.append(1000*recoStats[2])
			StdList.append(1000*recoStats[3])

		smearTimeList = [1000*x for x in smearTimeList]
		plt.plot(smearTimeList, AvList, 'k')
		plt.plot(smearTimeList, MedList, 'r')
		plt.ylabel("Difference from Truth (ps)")
		plt.xlabel("Smear Time (ps)")
		plt.legend(["Mean", "Median"])
		plt.show()

		plt.plot(smearTimeList, FWHMList, 'k')
		plt.plot(smearTimeList, StdList, 'r')
		plt.ylabel("Uncertainty on Reconstructed Time (ps)")
		plt.xlabel("Pixel Smear Time (ps)")
		plt.legend(["FWHM", "Std. Dev."])
		plt.show()

