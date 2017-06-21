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
	data = pickle.load(open("../../../data/pickles/pions/noB/AnaHit_Simu_pi-_10GeV_E30L_E10mm_H40L_H10mm.p", 'rb'))
	print "Done."
	sys.stdout.flush()


	#testing of algorithms given a line axis
	#data.events[12].projectionDisplay()
	#data.events[628].algo_rodLinearWithDepth()
	#smeardata = data.smear(0.01, 0)
	#smeardata.timeReco(algo=1, plotting=True)



	sys.exit()

	#visualizing display loop
	while(raw_input("> ") != 'q'):
		i = np.random.randint(0, len(data.events))
		print "on " + str(i)
		smeared = data.smear(0.01, 0)
		smeared.events[i].algo_rodLinearWithDepth()
		

