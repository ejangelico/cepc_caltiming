import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import DataSet

class Event:
	def __init__(self, hitPoints=None, hitEn=None, evNum=None):
		self.hitPoints = hitPoints  #array of hit points corresponding to pixel positions in the ECAL
		self.hitEn = hitEn  #array of energies that lines up with the hit points	
		self.evNum = evNum 	#integer id for the event

	def printEvent(self):
		print "Event Number:", self.evNum
		print
		print "Event Position:"
		for hit in self.hitPoints:
			print hit
		print
		print "Event Energy:", self.hitEn

	# Produces 3D event display of the pixels, where the color is the energy deposition
	def energyDisplay(self):
		x = []
		y = []
		z = []
		E = []
		for hit in self.hitPoints:
			x.append(hit[0])
			y.append(hit[1])
			z.append(hit[2])
		E = self.hitEn
		
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
			x.append(hit[0])
			y.append(hit[1])
			z.append(hit[2])
			t.append(hit[3])	
	
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
		for hit in self.hitPoints:
			t.append(hit[3])
		t = [x - min(t) for x in t]

		if rangeMin is None:
			rangeMin = min(t)
		if rangeMax is None:
			rangeMax = max(t)
		hist, bin_edges = np.histogram(t, numBins, (rangeMin, rangeMax), weights = self.hitEn)

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
			d.append(hit[1])
			t.append(hit[3])
		return (d, t)
