import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import DataSet

class Event:
	def __init__(self, hitPoints=None, evNum=None):
		self.hitPoints = hitPoints  #array of hit points corresponding to pixel positions in the ECAL	
		self.evNum = evNum 	#integer id for the event

	def printEvent(self):
		print "Event Number:", self.evNum
		print
		print "Event Position:"
		for hit in self.hitPoints:
			print hit
		print
		print "Event Energy:", self.hitEn

	# Returns an array of layers of all the points in the event
	def makeLayers(self):
		pass

	# Produces 3D event display of the pixels, where the color is the energy deposition
	def energyDisplay(self):
		x = []
		y = []
		z = []
		E = []
		for hit in self.hitPoints:
			x.append(hit.getX())
			y.append(hit.getY())
			z.append(hit.getZ())
			E.append(hit.getE()) 
		
		cm = plt.get_cmap('jet')
		cNorm = matplotlib.colors.Normalize(vmin=min(E), vmax=max(E))
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
		fig = plt.figure()
		ax = Axes3D(fig)
		ax.scatter(x, y, z, c=scalarMap.to_rgba(E))
		scalarMap.set_array(E)
		fig.colorbar(scalarMap)
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_zlabel("z")
		ax.set_xlim(-35, 35)
		ax.set_ylim(1840, 1910)
	        ax.set_zlim(-35, 35)
		plt.show()	

	# Produces 3D event display of the pixels, where the color is the time of the event
	def timeDisplay(self):
		x = []
		y = []
		z = []
		t = []
		for hit in self.hitPoints:
			x.append(hit.getX())
			y.append(hit.getY())
			z.append(hit.getZ())
			t.append(hit.getT())	
	
		cm = plt.get_cmap('jet')
		cNorm = matplotlib.colors.Normalize(vmin=min(t), vmax=max(t))
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
		fig = plt.figure()
		ax = Axes3D(fig)
		ax.scatter(x, y, z, c=scalarMap.to_rgba(t))
		scalarMap.set_array(t)
		fig.colorbar(scalarMap)
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_zlabel("z")
		ax.set_xlim(-35, 35)
		ax.set_ylim(1840, 1910)
	        ax.set_zlim(-35, 35)
                plt.show()	

	# Histograms the times of each pixel, weighted by the energy deposited.
	# The time is relative to the first hit
	# Returns (energyPerBin, binCenters)
	def timeHist(self, numBins, rangeMin = None, rangeMax = None):
		t = []
		E = []
		for hit in self.hitPoints:
			t.append(hit.getT())
			E.append(hit.getE())
		t = [x - min(t) for x in t]

		if rangeMin is None:
			rangeMin = min(t)
		if rangeMax is None:
			rangeMax = max(t)
		hist, bin_edges = np.histogram(t, numBins, (rangeMin, rangeMax), weights = E)

		binCenters = []
		for i in range(0, len(bin_edges)-1):
			binCenters.append((bin_edges[i]+bin_edges[i+1])/2.0)

		return  (hist, binCenters)

	# Plots the time histogram of the event
	def plotTimeHist(self, numBins):
		hits, binCenters = self.timeHist(numBins)
		plt.bar(binCenters, hits, width = binCenters[1]-binCenters[0], color = 'blue')
		plt.show()

	# Returns two arrays of the depth and time of each hit
	def timeVsDepth(self):
		d = []
		t = []
		for hit in self.hitPoints:
			d.append(hit.getRho())
			t.append(hit.getT())
		return (d, t)

	# Plot the time vs depth of every hit in this event
	def plotTvsD(self):
		d, t = data.events[110].timeVsDepth()
		plt.plot(d, t, 'ko')
		plt.xlabel("Depth (rho) into cal. (mm)")
		plt.ylabel("Time of hit (ns)")
		plt.show()
		
	# Does a linear fit to the first time of arrival vs depth in each layer
	def algo_linearFirstTimeByLayer(self):
		layers = self.makeLayers()
		tList = []
		dList = []
		for layer in layers:
			tList.append(layer.getFirstTime())
			dList.append(layer.d0)
