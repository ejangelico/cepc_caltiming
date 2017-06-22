import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
import cPickle as pickle 
from scipy.stats import norm
import time
import Point

import DataSet
import Event


def plotTimeHistogram(data1, data2, fig, ax):

	recoTimes1, eff1 = data1.listReconstructedTimes(algo=1)
	recoTimes2, eff2 = data2.listReconstructedTimes(algo=1)

	#best gaus fit
	(mu1, sig1) = norm.fit(recoTimes1)
	(mu2, sig2) = norm.fit(recoTimes2)
	n, bins1, _ = ax.hist(recoTimes1, 80, normed=1)
	y1 = mlab.normpdf(bins1, mu1, sig1)
	n, bins2, _ = ax.hist(recoTimes2, 80, normed=1)
	y2 = mlab.normpdf(bins2, mu2, sig2)
	ax.plot(bins1, y1)
	ax.plot(bins2, y2)


def getMisidentification(data1, data2):
	recoTimes1, eff1 = data1.listReconstructedTimes(algo=1)
	recoTimes2, eff2 = data2.listReconstructedTimes(algo=1)


	#make a cut down the middle in between the means of 
	#the two distributions. 
	mean1 = np.mean(recoTimes1)
	mean2 = np.mean(recoTimes2)

	print abs(mean1 - mean2)

	N1 = len(recoTimes1)
	N2 = len(recoTimes2)
	misidentified1 = 0 	#the number of particles type 1 that are on the particle 2 half
	misidentified2 = 0 	#the number of particles of type 2 that are on the particle 1 half

	if(mean1 < mean2):
		cutvalue = mean1 + abs(mean1 - mean2)/2.0 #time that splits the middle of the distributions
		for t in recoTimes2:
			if(t <= cutvalue):
				misidentified2 += 1
		for t in recoTimes1:
			if(t >= cutvalue):
				misidentified1 += 1

	elif(mean1 > mean2):
		cutvalue = mean2 + abs(mean1 - mean2)/2.0 #time that splits the middle of the distributions
		for t in recoTimes2:
			if(t >= cutvalue):
				misidentified2 += 1
		for t in recoTimes1:
			if(t <= cutvalue):
				misidentified1 += 1

	else:
		#the two means are equal, one cannot separate at all
		misidentified1 = N1
		misidentified2 = N2


	correct1 = N1 - misidentified1
	correct2 = N2 - misidentified2

	return (correct1, misidentified1, correct2, misidentified2)


def getOverlap(data1, data2):
	recoTimes1, eff1 = data1.listReconstructedTimes(algo=1)
	recoTimes2, eff2 = data2.listReconstructedTimes(algo=1)


	#This is the chunk of code for doing a binned separation
	#and counting area overlap by looking at overlapping bins

	#make a list with all of the times
	#this helps binning both distributions at once
	collectiveTimes = [_ for _ in recoTimes1]
	for t in recoTimes2:
		collectiveTimes.append(t)

	#this binning is super critical, figure out how to do it
	numBins = [20, 50, 100, 200, 400]

	for nb in numBins:
		tothist, bin_edges = np.histogram(collectiveTimes, nb)
		#bin both distributions the same
		hist1, bin_edges = np.histogram(recoTimes1, bin_edges, normed = 1)
		hist2, bin_edges = np.histogram(recoTimes2, bin_edges, normed = 1)

		binwidth = abs(bin_edges[0] - bin_edges[1])
		overlapArea = 0
		for i in range(len(bin_edges) - 1):
			if(hist1[i] > 0 and hist2[i] > 0):
				overlapArea += min(hist1[i], hist2[i])*binwidth
			else:
				overlapArea += 0

		print overlapArea




		

	




if __name__ == "__main__":

	#--Start data loading--#

	kaonpath = "../../../data/pickles/kaons/noB/"
	pionpath = "../../../data/pickles/pions/noB/"

	kfilenames = ["AnaHit_Simu_kaon-_2.5GeV_E30L_E10mm_H40L_H10mm.p"]
	pfilenames = ["AnaHit_Simu_pi-_2.5GeV_E30L_E10mm_H40L_H10mm.p"]

	kMomenta = [2.5]
	pMomenta = [2.5]


	kcorrect = []
	kmis = []
	pcorrect = []
	pmis = []

	smears = [0.001, 0.003, 0.009, 0.018, 0.036, 0.072, 0.144, 0.288, 0.576]
	for sm in smears:
		kdata = []
		pdata = []
		print "Loading kaon data"
		for kfn in kfilenames:
			print "Loading file " + str(kfn) + "...", 
			sys.stdout.flush()
			data = pickle.load(open(kaonpath+kfn, 'rb'))
			smdata = data.smear(sm, 0)
			print "Done."
			sys.stdout.flush()
			kdata.append(smdata)
		print "Done loading kaons"
		print "Loading pion data"
		for pfn in pfilenames:
			print "Loading file " + str(pfn) + "...", 
			sys.stdout.flush()
			data = pickle.load(open(pionpath+pfn, 'rb'))
			smdata = data.smear(sm, 0)
			print "Done."
			sys.stdout.flush()
			pdata.append(smdata)
		print "Done loading pions"

		for i in range(len(kMomenta)):
			kdata[i].setMomentum(kMomenta[i])
		for i in range(len(pMomenta)):
			pdata[i].setMomentum(pMomenta[i])

		#align these data lists in momentum
		kdata = sorted(kdata, key=lambda x: x.getMomentum())
		pdata = sorted(pdata, key=lambda x: x.getMomentum())

		#---End data loading---#

		a, b, c, d = getMisidentification(kdata[0], pdata[0])
		kcorrect.append(a/1000.0)
		kmis.append(b/(a + b))
		pcorrect.append(c/1000.0)
		pmis.append(d/(c + d))


	#rescale smears
	smears = [_*1000 for _ in smears]

	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols = 2, nrows = 2, figsize=(12, 7))
	ax1.plot(smears, kcorrect, 'mo-', label="Correctly identified Kaons")
	ax2.plot(smears, kmis, 'm^-', label="Misidentification fraction kaons")
	ax3.plot(smears, pcorrect, 'ro-', label="Correctly identified pions")
	ax4.plot(smears, pmis, 'r^-', label="Misidentification fraction pions")
	ax1.set_xlabel("Time smear (ps)")
	ax2.set_xlabel("Time smear (ps)")
	ax3.set_xlabel("Time smear (ps)")
	ax4.set_xlabel("Time smear (ps)")
	ax1.set_ylabel("Fraction of total particle generated")
	ax2.set_ylabel("Fraction of total particle generated")
	ax3.set_ylabel("Fraction of total particle generated")
	ax4.set_ylabel("Fraction of total particle generated")
	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')
	ax4.set_xscale('log')
	ax1.set_ylim([0, 1])
	ax2.set_ylim([0, 1])
	ax3.set_ylim([0, 1])
	ax4.set_ylim([0, 1])
	ax1.grid(True)
	ax2.grid(True)
	ax3.grid(True)
	ax4.grid(True)
	plt.savefig("misidentPlot_2.5GeV_700pscut.png", bbox_inches='tight')

	








	


	
		

