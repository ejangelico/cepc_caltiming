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
import Helper
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
		#ax.set_xlim([-70, 70])
		#ax.set_zlim([-70, 70])
		#ax.set_ylim([1848, 3350])
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
		ax.set_xlim([-30, 30])
		ax.set_zlim([-30, 30])
		ax.set_ylim([1848, 2050])
        
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
		fig, ax = plt.subplots()
		hits, binCenters = self.timeHist(numBins)
		ax.bar(binCenters, hits, width = binCenters[1]-binCenters[0], color = 'blue')
		plt.show()

	# Returns two arrays of the depth and time of each hit
	def timeVsDepth(self, plotting = False):
		d = []
		t = []
		for hit in self.hitPoints:
			d.append(hit.getRho())
			t.append(hit.getT())

		if(plotting == True):
			fig, ax = plt.subplots()
			ax.plot(d, t, 'ko')
			ax.set_xlabel("Depth into cal. (mm)")
			ax.set_ylabel("Time of hit (ns)")
			plt.show()

		return (d, t)


	#function that calculates the shower depth
	#of an event based on the definition of 
	#"Z0" from the CALICE paper 2014
	def getShowerDepth(self):
		rho_start = 1847.3 #the front face of the e-cal in rho (mm)
		timeCutoffLo = 6 # ns
		timeCutoffHi = 8 # ns 
		dCutoffLo = 1847.3 # mm
		dCutoffHi = 3385 # mm

		Z0 = 0
		esum = 0
		for hit in self.hitPoints:
			if(timeCutoffLo < hit.getT() < timeCutoffHi) and (dCutoffLo < hit.getRho() < dCutoffHi):
				e = hit.getE()
				esum += e
				rho = hit.getRho()
				Z0 += e*(rho - rho_start)

		Z0 = Z0/esum
		#divide by interaction length for pions in 
		#tungsten ~11.33cm
		Z0 = Z0/113.3

		return Z0

		

		
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
		dCutoffLo = 1847.3 # mm
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
			#string to hold the fit info
			fitinfo = 'Slope:' + str(fitParams[0]) + 'ns/mm\n'
			fitinfo += "1/Slope:" + str(1/fitParams[0]) + "mm/ns\n"
			fitinfo += "y-int:" + str(fitParams[1]) + "mm\n"
			fitinfo += "Estimate of shower start time:" + str(tEst) + "ns\n"
			fitinfo += "Truth value:" + str(min(tList)) + "ns\n"
			fitinfo += "Difference:" + str(np.abs(tEst - min(tList))) +  "ns\n"
			print fitinfo


			linFitFunc = np.poly1d(fitParams[:2])	
			x = np.linspace(min(dList), max(dList), 10)
			y = [linFitFunc(z) for z in x]	

			fig, ax = plt.subplots()
			Helper.resize(fig, ax)
			ax.plot(x, y, 'r', linewidth=2)
			ax.plot(dList, tList, 'ko', markersize=15)
			ax.set_xlabel("Depth (mm)")
			ax.set_ylabel("Time (ns)")
			plt.show()

		return tEst, min(tList)
