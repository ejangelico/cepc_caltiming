import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.stats import linregress
from scipy.optimize import least_squares
import numpy as np
import time
import sys
import Point
import DataSet
import HitPoint
import Point
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
	def makeLayersWithPoints(self, width):
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

	#makes layers of constant width at various
	#rho values throughout the ecal and hcal
	def makeLayersWithRadii(self, width):
		#make empty layers first
		layerList = []
		rad = 1847.3 + width/2.0 	#mm start of ecal + 1 width
		hcalrad = 3385
		while (rad < hcalrad):
			tempLay = Layer.Layer()
			tempLay.initializeWithRadius(rad, width)
			layerList.append(tempLay)
			rad += width

		#fill layers with the hitpoints
		for hitPoint in self.hitPoints:
			for layer in layerList:
				if layer.addPoint(hitPoint):
					break

		#remove layers with no hitpoints
		removes = []
		for l in layerList:
			if(l.hitPoints == None):
				removes.append(l)

		for rm in removes:
			layerList.remove(rm)


		return layerList

	# removes hit points in the event that do not
	# pass a set of cuts. Passing is constituted by 
	# the point lying in the range, for example trange
	# that is a two element list, trange = [floor, ceiling]
	#
	#remove = False/True option decides whether to actually remove
	#elements of the self.hitPoints list, or just to preserve them and 
	#return a new list of hits that passed
	def cutHitPoints(self, rhorange=None, phirange=None, zrange=None, trange=None, remove=False):
		removes = [] #list of hitpoints to remove
		passes = []
		for hp in self.hitPoints:
			#cumbersome if structure is so that
			#all combinations of ranges can be
			#allowed to be None or non-None
			if(rhorange == None):
				pass
			elif(min(rhorange) > hp.getRho() or max(rhorange) < hp.getRho()):
				removes.append(hp)
				continue

			if(phirange == None):
				pass
			elif(min(phirange) > hp.getPhi() or max(phirange) < hp.getPhi()):
				removes.append(hp)
				continue

			if(zrange == None):
				pass
			elif(min(zrange) > hp.getZ() or max(zrange) < hp.getZ()):
				removes.append(hp)
				continue

			if(trange == None):
				pass
			elif(min(trange) > hp.getT() or max(trange) < hp.getT()):
				removes.append(hp)
				continue

			#if it made it this far
			passes.append(hp)

		#remove the points that didn't pass
		if(remove==True):
			for rm in removes:
				self.hitPoints.remove(rm)

		return passes


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
		ecal_rin = 1847.3
		hcal_rin = 2058
		hcal_rout = 3385.5

		drawCylinder(ax, tpc_rin, tpc_ecal_hcal_z, 'b')
		drawCylinder(ax, tpc_rout, tpc_ecal_hcal_z, 'b')
		drawCylinder(ax, ecal_rin, tpc_ecal_hcal_z, 'y')
		drawCylinder(ax, hcal_rin, tpc_ecal_hcal_z, 'y')
		drawCylinder(ax, hcal_rout, tpc_ecal_hcal_z, 'r')

	#make 2 subplots that are 2D projections of eachother
	#line will plot a line
	#point will plot a big red point
	def projectionDisplay(self, line=None, point=None):
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
		cNorm = matplotlib.colors.Normalize(vmin=min(t), vmax=min(t)+6)
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
		fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(13, 6))


		ax1.scatter(x, y, c=scalarMap.to_rgba(t))
		ax1.set_xlim([-1000, 1000])
		ax1.set_ylim([1800, 2100])
		ax1.set_xlabel("x-projected axis")
		ax1.set_ylabel("y-projected axis")
		#hcal boundary
		ax1.plot([-1000, 1000], [2058, 2058], 'k-')

		ax2.scatter(z, y, c=scalarMap.to_rgba(t))
		ax2.set_xlim([-1000, 1000])
		ax2.set_ylim([1800, 2100])
		ax2.set_xlabel("z-projected axis")
		ax2.set_ylabel("y-projected axis")
		ax2.plot([-1000, 1000], [2058, 2058], 'k-')

		#plot a line on both plots as well
		if(line==None):
			pass
		else:
			#line is a point and a unit vector
			#line = [Point, Point] from the Point class

			#i use this to generate two points that I draw
			#a line between
			t1 = -10000
			t2 = 10000
			d0 = line[0] #the intercepts
			v = line[1]  #the slope parameter
			xx = [d0.getX() + v.getX()*t1, d0.getX() + v.getX()*t2]
			yy = [d0.getY() + v.getY()*t1, d0.getY() + v.getY()*t2]
			zz = [d0.getZ() + v.getZ()*t1, d0.getZ() + v.getZ()*t2]

			ax1.plot(xx, yy, 'r-')
			ax2.plot(zz, yy, 'r-')

		if(point == None):
			pass
		else:
			ax1.plot(point.getX(), point.getY(), 'kd', markersize=9)
			ax2.plot(point.getZ(), point.getY(), 'kd', markersize=9)




		scalarMap.set_array(t)
		plt.show()


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
		return (fig, ax)	

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
		rho_start = 1847.4 #the front face of the e-cal in rho (mm)
		timeCutoffLo = 6 # ns
		timeCutoffHi = 8 # ns 
		dCutoffLo = 1847.4 # mm
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


	#this is the first pass cut for a general
	#hadronic shower in the e-cal. It trims noise
	#and neutron fat to leave a nice shower core
	#Tested on 7.5 GeV pions with no B field
	def hadronicNoiseCut(self):
		#make a very rough noise cut

		#cut out late time first
		tlist = [_.getT() for _ in self.hitPoints]
		tmin = min(tlist)
		twidth = 6	
		trange = [tmin, tmin + twidth]
		cutpoints = self.cutHitPoints(rhorange=None, phirange=None, zrange=None, trange=trange, remove=False)

		#use those cut points to calculate spacial cuts
		philist = [_.getPhi() for _ in cutpoints]
		zlist = [_.getZ() for _ in cutpoints]
		meanphi = np.mean(philist)
		meanz = np.mean(zlist)
		#halfwidths
		zwidth = 250 #mm
		phiwidth = np.pi/8.0 
		#cutranges
		zrange = [meanz - zwidth, meanz + zwidth]
		phirange = [meanphi - phiwidth, meanphi + phiwidth]
		trimmedEvent = Event(cutpoints, 0)
		cutpoints = trimmedEvent.cutHitPoints(rhorange=None, phirange=phirange, zrange=zrange, trange=None, remove=False)
		trimmedEvent = Event(cutpoints, self.evNum)

		return trimmedEvent


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
		dCutoffLo = 1847.4 # mm
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

	def algo_rodLinearWithDepth(self):


		#perform initial rough cuts
		cutEvent = self.hadronicNoiseCut()
		if(len(cutEvent.hitPoints) == 0):
			print "did pass rough cut, adapt the cutting algorithm later"
			print "you need to be able to make a general rough cut on noise"
			return

		#here input the algorithm that finds
		#the shower axis
		showerAxis = [Point.Point(0,0,0,1), Point.Point(0,1,0,1).normalize()]
		#find the intersection point of this showerAxis with
		#the cylinder of the e-cal radius. Two points satisfy
		#equation
		point_plus, point_minus = Helper.getCylinderIntersection(1847.4, showerAxis)

		#find which point is closest to the earliest
		#point in the remaining trimmed hit points. This
		#should always be the correct intersection with the cylinder
		#because the intersections are on opposite poles. 
		firstHit = Helper.getFirstHit(cutEvent.hitPoints)
		firstPoint = Point.Point(firstHit.getX(), firstHit.getY(), firstHit.getZ(), 1)
		ecalIntersect = None
		if((firstPoint - point_plus).getMag() < (firstPoint - point_minus).getMag()):
			ecalIntersect = point_plus
		else:
			ecalIntersect = point_minus

	

		#here, find the best rod radius to use
		rodRadius = 15 	#mm


		passed = []		#hit points that are inside the rod
		rodDepths = []	#hit depths relative to the shower axis intersection with ecal
		rodTimes = []	#global hit time but for events that pass rod cut

		#two points on the line
		x1 = ecalIntersect
		x2 = x1 - showerAxis[1]
		for hp in cutEvent.hitPoints:
			x0 = Point.Point(hp.getX(), hp.getY(), hp.getZ(), 1)
			#find the distance of the perpendicular
			#to the axis line from hp
			d = ((x0 - x1).cross((x0 - x2))).getMag()/(x2 - x1).getMag()
			#if this distance is inside the rod radius
			#and the point is in the calorimeter "1847.4"mm
			if(d <= rodRadius and hp.getRho() >= 1847.4):
				passed.append(hp)
				#distance along axis from point x1 to 
				#the perpendicular intersection of axis
				#with hit point
				h = (x1 - x0).getMag()
				D = np.sqrt(h*h - d*d)
				#because x1 is the ecal Intersection,
				#D is the depth from the intersection
				rodDepths.append(D)
				rodTimes.append(hp.getT())
				
		if(len(passed) == 0):
			print "Rod captured no hits"
			print "Either the rod is too small or the showerAxis is poorly fit"
			return None


		#calculate fraction of total energy
		#contained in the rod--uses trimmed event so
		#as to ignore noise/etc
		etot = np.sum([_.getE() for _ in cutEvent.hitPoints])
		erod = np.sum([_.getE() for _ in passed])
		efrac = erod/etot


		#fit the depth vs time data with a line
		#that has slope of 1/speed of light. One
		#parameter fit
		"""
		#line with p as the y intercept
		c = 299.792458
		fitfunc = lambda p, x: (1.0/c)*x + p[0]
		errfunc = lambda p, x, y: fitfunc(p, x) - y
		pguess = [4]	#ns
		result = least_squares(errfunc, pguess, args=(np.array(rodDepths), np.array(rodTimes)))
		cept = result.x[0]
		"""

		






