import sys
import Event
import cPickle as pickle 
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import Helper

class DataSet:
	def __init__(self, events = None, tSmear = None, eSmear = None):
		self.events = events

		self.tSmear = tSmear 	#smear 1sigma in ns

		#fractional energy resolution
		#takes values between [0.0, 1.0]
		self.eSmear = eSmear

		self.layerBins = None

	def getAxis(self):
		pass

	# Plots the time-average histogram of all the events
	def avTimeHist(self, numBins, rangeMin, rangeMax):
		timeHist = [0 for _ in range(numBins)]
		for event in self.events:
			timeHist = [x+y for x, y in zip(timeHist, event.timeHist(numBins, rangeMin, rangeMax)[0])]
		
		binCenters = self.events[0].timeHist(numBins)[1]	

		fig, ax = plt.subplots()
		ax.bar(binCenters, timeHist, width = binCenters[1]-binCenters[0], color = 'blue')
		ax.set_xlabel("Time (nanoseconds)")
		ax.set_ylabel("Total shower energy/bin (GeV/ns)")
		plt.show()

	# Plots the depth (radius) vs time of every hit point
	def plotAllDvsT(self):
		d = []
		t = []
		for event in self.events:
			(d0, t0) = event.timeVsDepth()
			d += d0
			t += t0

		fig, ax = plt.subplots()
		ax.plot(d, t, 'k.')
		ax.set_xlabel("Depth into Cal. (mm)")
		ax.set_ylabel("Time of Hit (ns)")
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
		tFWHM = Helper.FWHM(tDiffBins, tDiffCounts)
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
			fig, ax = plt.subplots(figsize=(25, 13))
			ax.set_title("10 GeV e-, 10ps pixel resolution", fontsize=23)
			Helper.resize(fig, ax)
			tDiffBins = [x*1000 for x in tDiffBins]
			ax.bar(tDiffBins[:-1], tDiffCounts, width = tDiffBins[1]-tDiffBins[0], color = 'b')
			ax.set_xlabel("$t_{reco} - t_{true}$" + " (ps) ", fontsize = 20)
			ax.set_ylabel("Counts/bin for 1 GeV electrons", fontsize = 20)
			
			#plt.show()
			plt.savefig("../../../midterm_report/720plots/timereco_10Gev_10ps.png", bbox_inches='tight')

		return [tAv, tMed, tFWHM, tStd, tSkewness, tSkewTest]

	def timeRecoSmearing(self, smearTimeList = None, plotting = False):
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

		if plotting:
			fig, ax = plt.subplots(figsize=(25, 13))
			Helper.resize(fig, ax)
			smearTimeList = [1000*x for x in smearTimeList]
			ax.plot(smearTimeList, AvList, 'b', linewidth=3)
			ax.plot(smearTimeList, MedList, 'r', linewidth=3)
			ax.set_ylabel("Difference from Truth (ps)")
			ax.set_xlabel("Pixel Smear Time (ps)")
			ax.grid(True)
			ax.legend(["Mean", "Median"], fontsize=30)
			plt.savefig("../../../midterm_report/720plots/manypoints_means.png", bbox_inches='tight')
			p

			fig, ax = plt.subplots(figsize=(25, 13))
			Helper.resize(fig, ax)
			ax.plot(smearTimeList, FWHMList, 'b', linewidth=3)
			ax.plot(smearTimeList, StdList, 'r', linewidth=3)
			ax.set_ylabel("Uncertainty on Reconstructed Time (ps)")
			ax.set_xlabel("Pixel Smear Time (ps)")
			ax.legend(["FWHM", "Std. Dev."], fontsize=30)
			ax.grid(True)
			plt.savefig("../../../midterm_report/720plots/manypoints_stds.png", bbox_inches='tight')

		StdList = [x/1000 for x in StdList]
		return StdList
