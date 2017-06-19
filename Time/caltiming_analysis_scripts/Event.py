import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.stats import linregress
import numpy as np
import time
import sys
import DataSet
import HitPoint
import Layer

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
	def makeLayers(self, width):
		layerList = []
		for hitPoint in self.hitPoints:
			hitAdded = False
			for layer in layerList:
				if layer.addPoint(hitPoint):
					hitAdded = True
					break
			if not hitAdded:
				layerList.append(Layer.Layer())
				layerList[-1].initializeWithPoint(hitPoint, width)
		return layerList

	# Smear all the hit point times and energies by Gaussians with width tsm, esm, respectively
	def getSmearedEvent(self, tsm, esm):
		smearedHitPoints = []
		for hp in self.hitPoints:
			#smear with gaussian of sigma = tsm or esm
			if tsm > 0:
				newt = np.random.normal(hp.getT(), tsm)
			else:
				newt = hp.getT()
			
			if esm > 0:
				newe = np.random.normal(hp.getE(), esm*hp.getE())
			else:
				newe = hp.getE()

			#keep positions unsmeared
			x, y, z = hp.getXYZ()
			newHP = HitPoint.HitPoint(x, y, z, newt, newe, 1)
			smearedHitPoints.append(newHP)

		#return a new event with these hitpoints and
		#the same event number
		return Event(smearedHitPoints, self.evNum)

	#draws the detector as a set of cylinders
	#with hard coded radii. See
	#"GearOutput.xml" files for geometry values
	def drawDetector(self, ax):
		#L is half-z length
		def drawCylinder(ax, rad, L, color):
			#cylinder mesh
			x = np.linspace(-rad, rad, 100)
			z = np.linspace(-L, L)
			Xc, Zc = np.meshgrid(x, z)
			Yc = np.sqrt(rad**2 - Xc**2)

			#grid parameters
			rstride = 20
			cstride = 10
			ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride, color=color)
			ax.plot_surface(Xc, -Yc, -Zc, alpha=0.2, rstride=rstride, cstride=cstride, color=color)

		#half zs
		tpc_ecal_hcal_z = 2350

		#radii
		tpc_rin = 329
		tpc_rout = 1808
		ecal_rin = 1847.4
		hcal_rin = 2058
		hcal_rout = 3385.5

		drawCylinder(ax, tpc_rin, tpc_ecal_hcal_z, 'b')
		drawCylinder(ax, tpc_rout, tpc_ecal_hcal_z, 'b')
		drawCylinder(ax, ecal_rin, tpc_ecal_hcal_z, 'y')
		drawCylinder(ax, hcal_rin, tpc_ecal_hcal_z, 'y')
		drawCylinder(ax, hcal_rout, tpc_ecal_hcal_z, 'r')

	# Produces 3D event display of the pixels, where the color is the energy deposition
	def energyDisplay(self, drawDetector=False):
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
		if(drawDetector == True):
			self.drawDetector(ax)

		plt.show()	

	# Produces 3D event display of the pixels, where the color is the time of the event
	def timeDisplay(self, drawDetector=False):
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
		if(drawDetector == True):
			self.drawDetector(ax)
        
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
			d.append(hit.getY())
			t.append(hit.getT())
		return (d, t)

	# Plot the time vs depth of every hit in this event
	def plotTvsD(self):
		d, t = data.events[110].timeVsDepth()
		plt.plot(d, t, 'ko')
		plt.xlabel("Depth into cal. (mm)")
		plt.ylabel("Time of hit (ns)")
		plt.show()
		
	# Does a linear fit to the first time of arrival vs depth in each layer
	# layerWidth = width around center point of each layer, mm
	# timeCutoffLo = first hit time accepted by algo, ns
	# timeCutoffHi = last hit time accepted, ns
	# dCutoffLo = lowest layer position accepted, mm
	# dCutoffHi = highest layer position accepted, mm
	def algo_linearFirstTimeByLayer(self, layerWidth = 1.0, timeCutoffLo = 5.5, timeCutoffHi = 6.6, dCutoffLo = 1800, dCutoffHi = 1050, plotting = False):
		layerWidth = 1.0 # mm
		timeCutoffLo = 5.5 # ns
		timeCutoffHi = 6.6 # ns 
		dCutoffLo = 1800 # mm
		dCutoffHi = 2050 # mm

		layers = self.makeLayers(layerWidth)

		dList = []
		tList = []
		for layer in layers:
			dList.append(layer.d0)
			tList.append(layer.getFirstTime())
	
		tList_new = []
		dList_new = []
		for i in range(0, len(tList)):
			if (timeCutoffLo < tList[i] < timeCutoffHi) and (dCutoffLo < dList[i] < dCutoffHi):
				tList_new.append(tList[i])
				dList_new.append(dList[i])
		tList = tList_new
		dList = dList_new

		fitParams = linregress(dList, tList)
		
		tEst = fitParams[0]*min(dList)+fitParams[1]

		if plotting:
			print "Slope:", fitParams[0], "ns/mm"
			print "1/Slope:", 1/fitParams[0], "mm/ns"
			print "y-int:", fitParams[1], "mm"
			print "Estimate of shower start time:", tEst, "ns"
			print "Truth value:", min(tList), "ns"
			print "Difference:", np.abs(tEst - min(tList)), "ns"

			linFitFunc = np.poly1d(fitParams[:2])	
			x = np.linspace(min(dList), max(dList), 10)
			y = [linFitFunc(z) for z in x]	

			plt.plot(x, y, 'r')
			plt.plot(dList, tList, 'ko')
			plt.xlabel("Depth (mm)")
			plt.ylabel("Time (ns)")
			plt.show()

		return tEst, min(tList)
