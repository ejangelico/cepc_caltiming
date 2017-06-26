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
	median1 = np.median(recoTimes1)
	median2 = np.median(recoTimes2)


	N1 = len(recoTimes1)
	N2 = len(recoTimes2)
	misidentified1 = 0 	#the number of particles type 1 that are on the particle 2 half
	misidentified2 = 0 	#the number of particles of type 2 that are on the particle 1 half

	if(median1 < median2):
		cutvalue = median1 + abs(median1 - median2)/2.0 #time that splits the middle of the distributions
		for t in recoTimes2:
			if(t <= cutvalue):
				misidentified2 += 1
		for t in recoTimes1:
			if(t >= cutvalue):
				misidentified1 += 1

	elif(median1 > median2):
		cutvalue = median2 + abs(median1 - median2)/2.0 #time that splits the middle of the distributions
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


def gaussianFullDataSet(alg):
	kaonpath = "../../../data/pickles/kaons/noB/"
	pionpath = "../../../data/pickles/pions/noB/"

	kfilenames = ["AnaHit_Simu_kaon-_1GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_2.5GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_5GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_7.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_kaon-_10GeV_E30L_E10mm_H40L_H10mm.p"]
	pfilenames = ["AnaHit_Simu_pi-_1GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_2.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_7.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_10GeV_E30L_E10mm_H40L_H10mm.p"]

	momenta = [1, 2.5, 5, 7.5, 10]


	#structure of final data:
	#[[[data type 1 for time smear 1 as a function of energy], [data type 2 for time smear 1 ...], ...], [time smear 2]] 

	smears = [0.001, 0.003, 0.009, 0.018, 0.036, 0.072, 0.144, 0.288, 0.576]
	separationData = {}

	kdata = []
	pdata = []

	for kfn in kfilenames:
		print "Loading file " + str(kfn) + "...", 
		sys.stdout.flush()
		data = pickle.load(open(kaonpath+kfn, 'rb'))
		print "Done."
		sys.stdout.flush()
		kdata.append(data)
	print "Done loading kaons"
	print "Loading pion data"
	for pfn in pfilenames:
		print "Loading file " + str(pfn) + "...", 
		sys.stdout.flush()
		data = pickle.load(open(pionpath+pfn, 'rb'))
		print "Done."
		sys.stdout.flush()
		pdata.append(data)
	print "Done loading pions"

	for sm in smears:
		print "On timesmear " + str(sm)
		ksmdata = []
		psmdata = []
		for kd in kdata:
			ksmdata.append(kd.smear(sm, 0))
		for pd in pdata:
			psmdata.append(pd.smear(sm, 0))


		#---End data loading---#
		separationData[str(sm)] = []
		for i in range(len(momenta)):
			print "--On momentum " + str(momenta[i]) + "GeV"
			ksmdata[i].setMomentum(momenta[i])
			psmdata[i].setMomentum(momenta[i])

			t1, ef1 = ksmdata[i].listReconstructedTimes(alg)
			t2, ef2 = psmdata[i].listReconstructedTimes(alg)
			(mu1, sig1) = norm.fit(t1)
			(mu2, sig2) = norm.fit(t2)
			fig, ax = plt.subplots()
			n1, bin1, patches1 = ax.hist(t1, 100, normed=1, alpha=0.5)
			n2, bin2, patches2 = ax.hist(t2, 100, normed=1, alpha=0.5)
			y1 = mlab.normpdf(bin1, mu1, sig1)
			y2 = mlab.normpdf(bin2, mu2, sig2)
			ax.plot(bin1, y1, linewidth=2)
			ax.plot(bin2, y2, linewidth=2)
			plt.show()



		
def pickleFullDataSet(alg):
	kaonpath = "../../../data/pickles/kaons/noB/"
	pionpath = "../../../data/pickles/pions/noB/"

	kfilenames = ["AnaHit_Simu_kaon-_1GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_2.5GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_5GeV_E30L_E10mm_H40L_H10mm.p", "AnaHit_Simu_kaon-_7.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_kaon-_10GeV_E30L_E10mm_H40L_H10mm.p"]
	pfilenames = ["AnaHit_Simu_pi-_1GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_2.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_7.5GeV_E30L_E10mm_H40L_H10mm.p","AnaHit_Simu_pi-_10GeV_E30L_E10mm_H40L_H10mm.p"]

	momenta = [1, 2.5, 5, 7.5, 10]


	#structure of final data:
	#[[[data type 1 for time smear 1 as a function of energy], [data type 2 for time smear 1 ...], ...], [time smear 2]] 

	smears = [0.001, 0.003, 0.009, 0.018, 0.036, 0.072, 0.144, 0.288, 0.576]
	separationData = {}

	kdata = []
	pdata = []

	for kfn in kfilenames:
		print "Loading file " + str(kfn) + "...", 
		sys.stdout.flush()
		data = pickle.load(open(kaonpath+kfn, 'rb'))
		print "Done."
		sys.stdout.flush()
		kdata.append(data)
	print "Done loading kaons"
	print "Loading pion data"
	for pfn in pfilenames:
		print "Loading file " + str(pfn) + "...", 
		sys.stdout.flush()
		data = pickle.load(open(pionpath+pfn, 'rb'))
		print "Done."
		sys.stdout.flush()
		pdata.append(data)
	print "Done loading pions"

	for sm in smears:
		print "On timesmear " + str(sm)
		ksmdata = []
		psmdata = []
		for kd in kdata:
			ksmdata.append(kd.smear(sm, 0))
		for pd in pdata:
			psmdata.append(pd.smear(sm, 0))


		#---End data loading---#
		separationData[str(sm)] = []
		for i in range(len(momenta)):
			print "--On momentum " + str(momenta[i]) + "GeV"
			ksmdata[i].setMomentum(momenta[i])
			psmdata[i].setMomentum(momenta[i])

			N1, correct1, mis1, N2, correct2, mis2 = getMisidentification(ksmdata[i], psmdata[i], alg)

			momentumPointData = [N1, correct1, mis1, N2, correct2, mis2]
			separationData[str(sm)].append(momentumPointData)


	pickle.dump([momenta, smears, separationData], open("actualsnake_final.p", 'wb'))

			
def plotFullData(sepdatafile):
	momenta, smears, separationData = pickle.load(open(sepdatafile, 'rb'))

	fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(ncols = 2, nrows = 2, figsize=(40, 22))

	#heatmap calc
	def rgb(minimum, maximum, value):
	    minimum, maximum = float(minimum), float(maximum)
	    ratio = 2 * (value-minimum) / (maximum - minimum)
	    b = int(max(0, 255*(1 - ratio)))
	    r = int(max(0, 255*(ratio - 1)))
	    g = 255 - b - r
	    return r/255., g/255., b/255.

	vals = np.linspace(0, 255, len(smears))
	minrgb = 0
	maxrgb = 255
	heats = [rgb(minrgb, maxrgb, _) for _ in vals]
	heat_count = -1
	for sm in smears:
		heat_count += 1
		energyIndexedData = separationData[str(sm)]
		kCorrect = []
		pCorrect = []
		kmis = []
		pmis = []		
		for i in range(len(energyIndexedData)):
			N1, cor1, mis1, N2, cor2, mis2 = energyIndexedData[i]
			kCorrect.append(cor1/1000.0)
			kmis.append(mis1/1000.0)
			pCorrect.append(cor2/1000.0)
			pmis.append(mis1/1000.0)

		ax1.plot(momenta, kCorrect, marker='o', color=heats[heat_count], label=r"$\sigma_t = $" + str(sm*1000) + "ps", linewidth=3, markersize=10)
		ax2.plot(momenta, kmis, marker='o', color=heats[heat_count], label=r"$\sigma_t = $" + str(sm*1000) + "ps", linewidth=3, markersize=10)
		ax3.plot(momenta, pmis, marker='o', color=heats[heat_count], label=r"$\sigma_t = $" + str(sm*1000) + "ps", linewidth=3, markersize=10)
		ax4.plot(momenta, pCorrect, marker='o', color=heats[heat_count], label=r"$\sigma_t = $" + str(sm*1000) + "ps", linewidth=3, markersize=10)



	ax1.grid(True)
	ax2.grid(True)
	ax3.grid(True)
	ax4.grid(True)


	ax1.set_ylim([0, 1])
	ax2.set_ylim([0, 1])
	ax3.set_ylim([0, 1])
	ax4.set_ylim([0, 1])

	ax1.set_ylabel(r"K$_{true}$ . K$_{true}$", fontsize=22)
	ax2.set_ylabel(r"K$_{true}$ . K$_{false}$", fontsize=22)
	ax3.set_ylabel(r"K$_{false}$ . K$_{true}$", fontsize=22)
	ax4.set_ylabel(r"K$_{false}$ . K$_{false}$", fontsize=22)
	ax1.set_xlabel("momentum", fontsize=22)
	ax2.set_xlabel("momentum", fontsize=22)
	ax3.set_xlabel("momentum", fontsize=22)
	ax4.set_xlabel("momentum", fontsize=22)
	ax1.set_title(r"K$_{true}$ . K$_{true}$", fontsize=22)
	ax2.set_title(r"K$_{true}$ . K$_{false}$", fontsize=22)
	ax3.set_title(r"K$_{false}$ . K$_{true}$", fontsize=22)
	ax4.set_title(r"K$_{false}$ . K$_{false}$", fontsize=22)

	ax1.legend(fontsize=18)
	ax2.legend(fontsize=18)
	ax3.legend(fontsize=18)
	ax4.legend(fontsize=18)

	majorLocatorY = MultipleLocator(0.2)
	minorLocatorY = MultipleLocator(0.05)
	majorLocatorX = MultipleLocator(2)
	minorLocatorX = MultipleLocator(1)
	ax1.get_xaxis().set_major_locator(majorLocatorX)
	ax1.get_xaxis().set_minor_locator(minorLocatorX)
	ax1.get_yaxis().set_major_locator(majorLocatorY)
	ax1.get_yaxis().set_minor_locator(minorLocatorY)
	ax1.get_xaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax1.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax1.get_yaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax1.get_yaxis().set_tick_params(length=10, width=2, which='minor')

	ax2.get_xaxis().set_major_locator(majorLocatorX)
	ax2.get_xaxis().set_minor_locator(minorLocatorX)
	ax2.get_yaxis().set_major_locator(majorLocatorY)
	ax2.get_yaxis().set_minor_locator(minorLocatorY)
	ax2.get_xaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax2.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax2.get_yaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax2.get_yaxis().set_tick_params(length=10, width=2, which='minor')

	ax3.get_xaxis().set_major_locator(majorLocatorX)
	ax3.get_xaxis().set_minor_locator(minorLocatorX)
	ax3.get_yaxis().set_major_locator(majorLocatorY)
	ax3.get_yaxis().set_minor_locator(minorLocatorY)
	ax3.get_xaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax3.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax3.get_yaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax3.get_yaxis().set_tick_params(length=10, width=2, which='minor')

	ax4.get_xaxis().set_major_locator(majorLocatorX)
	ax4.get_xaxis().set_minor_locator(minorLocatorX)
	ax4.get_yaxis().set_major_locator(majorLocatorY)
	ax4.get_yaxis().set_minor_locator(minorLocatorY)
	ax4.get_xaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax4.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax4.get_yaxis().set_tick_params(labelsize=22, length=20, width=2, which='major')
	ax4.get_yaxis().set_tick_params(length=10, width=2, which='minor')


	plt.savefig("test.png", bbox_inches='tight')





	




if __name__ == "__main__":


	#plotFullData("highway_firstTry.p")
	#pickleFullDataSet(1)
	#gaussianFullDataSet(2)
	#sys.exit()

	#--Start data loading--#

	kaonpath = "../../../data/pickles/kaons/noB/"
	pionpath = "../../../data/pickles/pions/noB/"

	kfilenames = ["AnaHit_Simu_kaon-_2.5GeV_E30L_E10mm_H40L_H10mm.p"]
	pfilenames = ["AnaHit_Simu_pi-_2.5GeV_E30L_E10mm_H40L_H10mm.p"]

	kMomenta = [2.5]
	pMomenta = [2.5]

	#smears = [0.010]

	kCorrect = []
	pCorrect = []
	kmis = []
	pmis = []

	smears = [0.1]
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

		t1, ef1 = kdata[0].listReconstructedTimes(algo=2)
		t2, ef2 = pdata[0].listReconstructedTimes(algo=2)

		mm = np.median(t2) + 0.5*abs(np.median(t1) - np.median(t2))
		fig, ax = plt.subplots(figsize=(13, 7))
		ax.hist(t1, 1000, alpha=.5, label="kaons")
		ax.hist(t2, 1000, alpha=.5, label="pions")
		ax.plot([mm, mm], [0, 100])
		ax.plot([np.median(t1), np.median(t1)], [0, 100])
		ax.plot([np.median(t2), np.median(t2)], [0, 100])
		ax.set_xlim([-0.3, 0.05])
		ax.legend()
		ax.set_title("Simple 2.5GeV k/pi 100ps smear")
		plt.show()
		sys.exit()


		N1, cor1, mis1, N2, cor2, mis2 = getMisidentification(kdata[0], pdata[0], 2)
		kCorrect.append(cor1/1000.0)
		kmis.append(float(mis1)/(cor2 + mis1))
		pCorrect.append(cor2/1000.0)
		pmis.append(float(mis2)/(cor1 + mis2))
		treco1, _ = kdata[0].listReconstructedTimes(algo=2)
		treco2, _ = pdata[0].listReconstructedTimes(algo=2)

		fig, ax = plt.subplots()
		mean1 = np.mean(treco1)
		mean2 = np.mean(treco2)
		mm = abs(mean1 + mean2)/2.0
		ax.hist(treco1, 500, alpha=0.4)
		ax.hist(treco2, 500, alpha=0.4)
		ax.plot([mm + mean2, mm + mean2], [0, 700])
		plt.show()
		sys.exit()



	fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(13, 7))
	ax1.plot(smears, kCorrect, 'mo-', label='Correctly identified kaons')
	ax2.plot(smears, kmis, 'md-', label="Misidentified fraction kaons")
	ax3.plot(smears, pCorrect, 'ro-', label="Correctly identified pion")
	ax4.plot(smears, pmis, 'rd-', label="Misidentified fraction pions")

	ax1.grid(True)
	ax2.grid(True)
	ax3.grid(True)
	ax4.grid(True)

	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')
	ax4.set_xscale('log')

	ax1.set_ylim([0, 1])
	ax2.set_ylim([0, 1])
	ax3.set_ylim([0, 1])
	ax4.set_ylim([0, 1])

	ax1.set_ylabel("Fraction of total events generated")
	ax2.set_ylabel("Misidentified/(cut + misidentified")
	ax3.set_ylabel("Fraction of total events generated")
	ax4.set_ylabel("Misidentified/(cut + misidentified")

	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()

	plt.savefig("algosimple_misidentified_timesmearpscut_75GeV.png", bbox_inches='tight')
		









	


	
		

