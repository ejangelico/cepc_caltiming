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

	def timeHist(self, numBins):
		t = []
		for hit in self.hitPoints:
			t.append(hit[3])
		t = [x - min(t) for x in t]
		#startIndex = t.index(0) #Remove 0 to see if initial peak goes away
		#t.pop(startIndex)
		#if startIndex == 0: E = self.hitEn[1:]
		#else: E = self.hitEn[:startIndex-1] + self.hitEn[startIndex:]
		return np.histogram(t, numBins, weights = self.hitEn)

	def plotTimeHist(self, numBins):
		hits, bin_edges = self.timeHist(numBins)
		times = []
		for i in range(0, len(bin_edges)-1):
			times.append((bin_edges[i]+bin_edges[i+1])/2.0)
		plt.bar(times, hits, width = times[1]-times[0], color = 'blue')
		plt.show()
