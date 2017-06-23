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



def getMisidentification(data1, data2, alg):
	recoTimes1, eff1 = data1.listReconstructedTimes(algo=alg)
	recoTimes2, eff2 = data2.listReconstructedTimes(algo=alg)


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

	return (N1, correct1, misidentified1, N2, correct2, misidentified2)



		
def pickleFullDataSet():
	kaonpath = "../../../data/pickles/kaons/noB/"
	pionpath = "../../../data/pickles/pions/noB/"

	kfilenames = ["AnaHit_Simu_kaon-_1GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_2.5GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_5GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_7.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_kaon-_10GeV_E30L_E10mm_H40L_H10mm.p"]
	pfilenames = ["AnaHit_Simu_pi-_1GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_2.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_7.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_10GeV_E30L_E10mm_H40L_H10mm.p"]

	momenta = [1, 2.5, 5, 7.5, 10]


	#structure of final data:
	#[[[data type 1 for time smear 1 as a function of energy], [data type 2 for time smear 1 ...], ...], [time smear 2]] 

	smears = [0.001, 0.003, 0.009, 0.018, 0.036, 0.072, 0.144, 0.288, 0.576]
	separationData = {}

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

		#---End data loading---#
		separationData[str(sm)] = []
		for i in range(len(momenta)):
			print "--On momentum " + str(momenta[i]) + "GeV"
			kdata[i].setMomentum(momenta[i])
			pdata[i].setMomentum(momenta[i])

			N1, correct1, mis1, N2, correct2, mis2 = getMisidentification(kdata[i], pdata[i])

			momentumPointData = [N1, correct1, mis1, N2, correct2, mis2]
			separationData[str(sm)].append(momentumPointData)


	pickle.dump([momenta, smears, separationData], open("snakeperformance_fulldata_fixedTimeCut_15mmRod.p", 'wb'))

			
def plotFullData(sepdatafile):
	momenta, smears, separationData = pickle.load(open(sepdatafile, 'rb'))

	fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(ncols = 2, nrows = 2, figsize=(13, 7))

	for sm in smears:
		energyIndexedData = separationData[str(sm)]




	




if __name__ == "__main__":


	#pickleFullDataSet()
	#sys.exit()

	#--Start data loading--#

	kaonpath = "../../../data/pickles/kaons/noB/"
	pionpath = "../../../data/pickles/pions/noB/"

	kfilenames = ["AnaHit_Simu_kaon-_2.5GeV_E30L_E10mm_H40L_H10mm.p"]
	pfilenames = ["AnaHit_Simu_pi-_2.5GeV_E30L_E10mm_H40L_H10mm.p"]

	kMomenta = [2.5]
	pMomenta = [2.5]

	#smears = [0.010]

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

		#---End data loading---#

		t0pi = pdata[0].simpleReco()
		t0k = kdata[0].simpleReco()

		collective = t0pi
		for t in t0k:
			collective.append(t)
		fig, ax = plt.subplots()
		#ax.hist(t0pi, 300, alpha=0.5)
		#ax.hist(t0k, 300, alpha=0.5)
		ax.hist(collective, 300)
		plt.show()









	


	
		

