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
	data = pickle.load(open("../../../data/pickles/kaons/noB/AnaHit_Simu_kaon-_1GeV_E30L_E10mm_H40L_H10mm.p", 'rb'))
	print "Done."
	sys.stdout.flush()


	#testing of algorithms given a line axis
	#data.events[14].projectionDisplay()



	#sys.exit()

	#visualizing display loop
	for i in range(len(data.events)):
		print "on " + str(i)
		smeared = data.smear(0.0, 0)
		smeared.events[i].algo_rodLinearWithDepth()
		if(raw_input(">  ") == 'q'):
			break

