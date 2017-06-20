import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
import cPickle as pickle 
import time

import DataSet
import Event

if __name__ == "__main__":

	print "Loading file...", 
	sys.stdout.flush()
	data = pickle.load(open("../pickles/pions/noB/200Gev_20evts.p", 'rb'))
	print "Done."
	sys.stdout.flush()

	#data.events[10].timeVsDepth(True)
	#data.events[10].energyDisplay(False)
	#data.events[12].algo_linearFirstTimeByLayer(plotting=True)
	#smdata = data.smear(0.01, 0)
	#data.timeReco(0, False)
	#data.timeRecoSmearing(np.linspace(0, 0.5, 30))

	layerWidth = 10	#mm
	layers = data.events[10].makeLayersWithRadii(layerWidth)
	layers = sorted(layers, key=lambda x: x.getCenter())
	firstHit = layers[0].getFirstPoint()

	#cuts
	phiwidth = np.pi/3.0
	zwidth = 20 	#mm
	zrange = [firstHit.getZ() - zwidth, firstHit.getZ() + zwidth]
	phirange = [firstHit.getPhi() - phiwidth, firstHit.getPhi() + phiwidth]

	#layer coloring
	count = 0
	newHPs = []

	#centroid testing
	centroids = []
	for l in layers:
		count += 1
		l.cutHitPoints(rhorange=None, phirange=phirange, zrange=zrange, trange=None)
		if(l.getSpaceTimeCentroid() == None):
			pass
		else:
			centroids.append(l.getSpaceTimeCentroid())

		for hp in l.hitPoints:
			hp.setE(count)
			newHPs.append(hp)

	trimmedEvent = Event.Event(newHPs, 0)
	cx = [_.getX() for _ in centroids]
	cy = [_.getY() for _ in centroids]
	cz = [_.getZ() for _ in centroids]
	#trimmedEvent.energyDisplay(True)
	fig = plt.figure(figsize=(10,7))
	ax = fig.gca(projection='3d')
	ax.scatter(cx, cy, cz, s=5)
	plt.show()





