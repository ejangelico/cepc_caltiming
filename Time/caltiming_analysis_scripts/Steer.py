import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
import cPickle as pickle 
import time
import Point

import DataSet
import Event

if __name__ == "__main__":

	print "Loading file...", 
	sys.stdout.flush()
	data = pickle.load(open("../../../data/pickles/kaons/noB/AnaHit_Simu_kaon-_5GeV_E30L_E10mm_H40L_H10mm.p", 'rb'))
	print "Done."
	sys.stdout.flush()


	#testing of algorithms given a line axis
	#data.events[12].projectionDisplay()
	#data.events[900].algo_rodLinearWithDepth(0.01)
	smeardata = data.smear(0.1, 0)
	#smeardata.setMomentum(10)
	#smeardata.events[12].algo_rodLinearWithDepth(0.01)
	#smeardata.events[50].algo_Snake(0)
	#smeardata.testSnake()
	#smeardata.timeReco(1, True)
	#smeardata.simpleReco()
	t0, eff = smeardata.listReconstructedTimes(algo=1)
	fig, ax = plt.subplots(figsize=(13, 7))
	ax.hist(t0, 300)
	ax.set_xlim([5.9, 6.3])
	ax.set_xlabel("reconstructed time (ns)")
	ax.set_title("5 GeV Kaon time reconstruction, 100ps smear")
	plt.show()




	sys.exit()

	#visualizing display loop
	while(True):
		i = np.random.randint(0, len(data.events))
		print "on " + str(i)
		smeared = data.smear(0.01, 0)
		smeared.events[i].algo_Snake(0.01)
		

