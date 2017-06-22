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
	data = pickle.load(open("../pickles/pions/noB/5GeV_100.p", 'rb'))
	print "Done."
	sys.stdout.flush()

	event = data.events[11].hadronicNoiseCut()
	dotArray = []
	yhat = Point.Point(0, 1, 0, cart = True)
	w0List = np.linspace(0.5, 100, 1000)
	for w0 in w0List:
		print round(w0, 4)
		ShowerAx = event.getShowerAxisWeighted(w0)
		dotArray.append(yhat*ShowerAx[1])
	#event.projectionDisplay(line = ShowerAx)
	plt.plot(w0List, dotArray, 'k')
	plt.xlabel("w0")
	plt.ylabel("Dot product from y-axis")
	plt.show()

	sys.exit()

	#testing of algorithms given a line axis
	#data.events[12].projectionDisplay()
	#data.events[5].algo_rodLinearWithDepth()



	#sys.exit()

	#visualizing display loop
	while(raw_input("> ") != 'q'):
		i = np.random.randint(0, len(data.events))
		print "on " + str(i)
		smeared = data.smear(0.01, 0)
		smeared.events[i].algo_rodLinearWithDepth()
		

