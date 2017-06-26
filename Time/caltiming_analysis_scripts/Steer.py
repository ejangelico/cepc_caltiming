import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
import cPickle as pickle 
import time
import Point
import Octagon
import DataSet
import Event

if __name__ == "__main__":
	Oct = Octagon.Octagon(R = 1847.4)

	print "Loading file...", 
	sys.stdout.flush()
	data = pickle.load(open("../pickles/kaons/noB/AnaHit_Simu_kaon-_5GeV_E30L_E10mm_H40L_H10mm.p", 'rb'))
	print "Done."
	sys.stdout.flush()

	data = data.smear(0.05, 0)
	
	for i in range(0, len(data.events)):
		print "Event number:", i
		print data.events[i].algo_Highway(rodRadius = 15, showerAxis = data.getAxis(isB = False), plotting = True)
	sys.exit()

	for i in range(0, len(data.events)):
		event = data.events[i].hadronicNoiseCut()
		event.projectionDisplay(line = data.getAxis(isB = True))

	"""
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
	"""
	sys.exit()

	#testing of algorithms given a line axis
	#data.events[12].projectionDisplay()
	#data.events[900].algo_rodLinearWithDepth(0.01)
	#smeardata = data.smear(0.0, 0)
	#smeardata.events[12].algo_rodLinearWithDepth(0.01)
	#smeardata.timeReco(algo=1, plotting=True)
	data.energyDepthPlot()



	sys.exit()

	#visualizing display loop
	while(raw_input("> ") != 'q'):
		i = np.random.randint(0, len(data.events))
		print "on " + str(i)
		smeared = data.smear(0.01, 0)
		smeared.events[i].algo_rodLinearWithDepth()
		

