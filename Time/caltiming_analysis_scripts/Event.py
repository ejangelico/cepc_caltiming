import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from matplotlib.colors import LogNorm
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
import Point
import Octagon


global ecalRIN 
global speedOfLight

ecalRIN = 1847.4 #mm
speedOfLight = 299.792458 #mm/ns


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
		rad = ecalRIN + width/2.0 	#mm start of ecal + 1 width
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
		ecal_rin = ecalRIN
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
			t.append(hit.getE())
		
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
		Oct = Octagon.Octagon(1847.4)
		Oct.plot(ax1)

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
		rho_start = ecalRIN #the front face of the e-cal in rho (mm)
		timeCutoffLo = 6 # ns
		timeCutoffHi = 8 # ns 
		dCutoffLo = ecalRIN # mm
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


	# Returns [Point(x0, y0, z0), Point(vx, vy, vz)]
	def getShowerAxis(self):
		points = []
		for hitPoint in self.hitPoints:
			points.append([hitPoint.getX(), hitPoint.getY(), hitPoint.getZ()])
		return Helper.LinFit3D(points)
	# Returns [Point(x0, y0, z0), Point(vx, vy, vz)]
	def getShowerAxisWeighted(self, w0 = 2):
		points = []
		for hitPoint in self.hitPoints:
			points.append([hitPoint.getX(), hitPoint.getY(), hitPoint.getZ(), hitPoint.getE()])
		return Helper.LinFit3DWeighted(points, w0)

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

	#does a rod cut on all hit points in the event
	#given a shower axis. 
	#Returns: a list of hit points that passed,
	#a list of depths relative to the ecal Intersection point
	#and the shower axis, and their times
	def rodFilter(self, radius, showerAxis):


		#start the rod filtering
		rodRadius = radius 	#mm
		passed = []		#hit points that are inside the rod
		rodDepths = []	#hit depths relative to the shower axis intersection with ecal
		rodTimes = []	#global hit time but for events that pass rod cut

		#two points on the line
		x1 = showerAxis[0] 	#intersection point with e-cal
		x2 = x1 - showerAxis[1]
		for hp in self.hitPoints:
			x0 = Point.Point(hp.getX(), hp.getY(), hp.getZ(), 1)
			#find the distance of the perpendicular
			#to the axis line from hp
			d = ((x0 - x1).cross((x0 - x2))).getMag()/(x2 - x1).getMag()
			#if this distance is inside the rod radius
			#and the point is in the calorimeter "1847.4"mm
			if(d <= rodRadius and hp.getRho() >= ecalRIN):
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


		return (passed, rodDepths, rodTimes)



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
		dCutoffLo = ecalRIN # mm
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

	def algo_Snake(self, timesmear, radius, showerAxis):

		#perform initial rough cuts
		cutEvent = self.hadronicNoiseCut()
		if(len(cutEvent.hitPoints) <= 7):
			#print "did pass rough cut, adapt the cutting algorithm later"
			#print "you need to be able to make a general rough cut on noise"
			return (None, 1)

		passed, rodDepths, rodTimes = cutEvent.rodFilter(radius, showerAxis)

		if(len(passed) <= 8):
			return (None, 2)

		rodEvent = Event(passed, 0)
		#rodEvent.projectionDisplay(showerAxis, ecalIntersect)

		#---BEGIN N-hit iteration fitting---#
		#order the passed hits based on time
		timebank, depthbank = (list(t) for t in zip(*sorted(zip(rodTimes, rodDepths), key=lambda x: x[0])))
		tcept = []
		tcept_av = []
		tcept_std = []
		nit = []
		fitvels = []
		differentialTCept = []
		costs = []

		#delete this
		ddd = []
		ttt = []

		#iterate through the points
		#and fit at each iteration
		n = 0
		while n < (len(timebank) - 2):
			fittimes = [timebank[i] for i in range(n + 3)]
			fitdepths = [depthbank[i] for i in range(n + 3)]

			#fit these points to a line with floating slope
			fitfunc = lambda p, x: p[1]*x + p[0]
			errfunc = lambda p, x, y: fitfunc(p, x) - y
			pguess = [4, 1.0/299.0]	#ns, ns/mm
			result = least_squares(errfunc, pguess, args=(np.array(fitdepths), np.array(fittimes)))
			cept = result.x[0]
			fitvel = result.x[1]

			#the meat of the algorithm, 
			#rejecting additional points 
			#based on a criteria
			if(n > 0):
				if(timesmear <= 0.05):
					pscut = 0.05/(n**(0.5))
				else:
					pscut = 1.5*timesmear/(n**(0.5))
				if(abs(cept - tcept[-1]) > pscut):
					#skip this point by removing
					#from the bank
					timebank.remove(fittimes[-1])
					depthbank.remove(fitdepths[-1])
					continue


			#push the results
			nit.append(n)
			tcept.append(cept)
			tcept_av.append(np.mean(tcept))
			tcept_std.append(np.std(tcept))
			fitvels.append(fitvel)
			costs.append(result.cost)
			if(n > 0):
				differentialTCept.append(tcept[-1] - tcept[-2])

			n += 1


		#----------------------#

		#---Final outlier cuts---#

		#if you haven't used at least 7
		#points to do the fit
		if(len(timebank) < 7):
			return (None, 3)


		#for testing cutting efficiencies
		#return (len(passed), len(timebank))
		#-----------------------#

		#ttrue = 6.590 #ns for 1GeV charged kaon
		"""
		
		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(13,7))
		ax1.plot(nit, tcept, 'ro-')
		ax1.set_xlabel("number of iterations")
		ax1.set_ylabel("time intercept of linear fit")
		ax2.plot(rodDepths, rodTimes, 'ko')
		ax2.plot(rodDepths, fitfunc([tcept[-1], fitvels[-1]], np.array(rodDepths)), 'b-', label="N used = " + str(len(timebank)))
		ax2.plot(depthbank, timebank, 'ro')
		ax2.legend()
		ax2.set_xlabel("rod depth (mm)")
		ax2.set_ylabel("hit time (ns)")
		ax3.plot(nit[1:], differentialTCept, 'go-')
		ax3.set_ylabel("differenetial time intercept")
		ax4.plot(nit, tcept_std, 'mo-')
		ax4.set_ylabel("running standard deviation")
		plt.show()

		cutEvent.projectionDisplay()

		print self.evNum
		"""

		return (tcept[-1], 0)


	def algo_Highway(self, rodRadius, showerAxis, plotting = False):
		cLight = speedOfLight # mm/ns

		# Hit cut and event cut
		eventCut = self.hadronicNoiseCut()
		passedHits, rodDepths, rodTimes = eventCut.rodFilter(radius = 15, showerAxis = showerAxis)
		if len(passedHits) < 3:
			return None

		# The main highway algorithm

		# Parameters of the algorithm
		dCut = 1.0    # Initial width of the band
		dStep = 0.01 # How much to step dCut by
		numIter = int(dCut/dStep)+1 # Maximum number of iterations
		cutPointCountLength = 10 # Stores the number of points cut over the last cutPointCountLength iterations
		derivThreshold = 5 # If the sum of the number of points cut over the last cutPointCountLength iterations
				   # is greater than derivThreshold, then stop iterating and return the time from
				   # cutPointCountLength iterations ago

		# Sort hits by time
		times, depths = (list(t) for t in zip(*sorted(zip(rodTimes, rodDepths), key=lambda x: x[0])))
		
		# Fit functions for the linear fit
		lineFunc = lambda p, x: p[1]*x+p[0]
		resFunc  = lambda p, x, y: lineFunc(p, x) - y

		line = [times[0]-depths[0]/cLight, 1/cLight] # Initial Guess, [yint, m] for line params

		# The variables that stores the points left after band cut
		timesRemaining = times
		depthsRemaining = depths
		
		# The quantities to measure at each iteration, not including the first guess
		t0List = []
		cutPointCount = []

		if plotting:
			numPointsList = []
			derivList = []
			pointsList = []
			lineList = []
			dCutList = []

		for n in range(0, numIter):
			# The length before this iteration's band cut
			lengthPrev = len(timesRemaining)

			# The list of times/depths after the band cut
			timesRemaining = []
			depthsRemaining = []

			# Do the band cut
			for i in range(0, len(times)):
				# This is the perpendicular distance of a point from the line
				if np.abs((times[i] - (line[0] + line[1]*depths[i])))/np.sqrt(1+line[1]**2) < dCut:
					timesRemaining.append(times[i])
					depthsRemaining.append(depths[i])

			# If there are only two points left, stop the algorithm, it failed
			if len(timesRemaining) <= 2:
				#print "Highway failed."
				return None

			# If there are more than two points, do a least squares fit
			result = least_squares(resFunc, [line[0], line[1]], args = (np.array(depthsRemaining), np.array(timesRemaining)))
			t0, m = result.x

			# Update the line parameters and the band width
			line[0] = t0
			line[1] = m
			dCut -= dStep

			# Update the t0List
			t0List.append(t0)

			# Check if losing points too quickly and return t0 if so
			numPointsCut = lengthPrev - len(timesRemaining)
			cutPointCount.append(numPointsCut)

			# Update all the parameters to plot
			if plotting:
				numPointsList.append(len(timesRemaining))
				derivList.append(numPointsCut)
				pointsList.append([depthsRemaining, timesRemaining])
				lineList.append([t0, m])
				dCutList.append(dCut+dStep)

			if len(cutPointCount) > cutPointCountLength:
				del cutPointCount[0]
			if sum(cutPointCount) >= derivThreshold:
				if not plotting:
					if(len(t0List) < cutPointCountLength):
						return t0List[0]
					else:
						return t0List[-cutPointCountLength]
				else:
					break

		if plotting:
			iterCount = range(0, len(t0List))
		
			fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(10,7))
		
			ax1.plot(iterCount, t0List, 'k')
			ax1.set_xlabel("Iteration")
			ax1.set_ylabel("$t_0$ (ns)")
		
			ax2.plot(iterCount, numPointsList, 'k')
			ax2.set_xlabel("Iteration")
			ax2.set_ylabel("Number of Points")

			ax3.plot(iterCount, derivList, 'k')
			ax3.set_xlabel("Iteration")
			ax3.set_ylabel("Number of Points Cut")

			t0 = lineList[-cutPointCountLength][0]
			m  = lineList[-cutPointCountLength][1]
			xPoints = pointsList[-cutPointCountLength][0]
			yPoints = pointsList[-cutPointCountLength][1]
			dCut = dCutList[-cutPointCountLength]

			x = np.linspace(min(xPoints), max(xPoints))
			y = [m*z +t0 for z in x]
			yHi = [m*z + t0 + dCut*np.sqrt(1+m**2) for z in x]
			yLo = [m*z + t0 - dCut*np.sqrt(1+m**2) for z in x]

			ax4.plot(depths, times, 'ro')
			ax4.plot(xPoints, yPoints, 'ko')
			ax4.plot(x, y, 'k')
			ax4.plot(x, yLo, 'k--')
			ax4.plot(x, yHi, 'k--')
			ax4.set_xlabel("Depth along Cylinder (mm)")
			ax4.set_ylabel("Time (ns)")

			plt.show()

			return t0List[-cutPointCountLength]

		# If you make it outside the for loop, there has been some strange failure
		return None
		

	def algo_Huff(self, timesmear, showerAxis):
		#perform initial rough cuts
		cutEvent = self.hadronicNoiseCut()
		if(len(cutEvent.hitPoints) <= 5):
			#print "did pass rough cut, adapt the cutting algorithm later"
			#print "you need to be able to make a general rough cut on noise"
			return (None, 1)

		#here input the algorithm that finds
		#the shower axis.

		#cutEvent.projectionDisplay(showerAxis)
		radius = 15 #mm
		passed, rodDepths, rodTimes = cutEvent.rodFilter(radius, showerAxis)

		if(len(passed) <= 5):
			return (None, 2)

		rodEvent = Event(passed, 0)
		#rodEvent.projectionDisplay(showerAxis, ecalIntersect)

		#---BEGIN Hough transform fitting---#
		
		#parametrizes a point (d, t)
		#returns an array of rho's
		#and thetas
		def hough(d, t):
			thetas = np.linspace(0, np.pi, 1000, endpoint=False)
			rhos = [d*np.cos(x) + t*np.sin(x) for x in thetas]
			return (thetas, rhos)

		#take in a th and rho from hough
		#transform and return b, m
		#for y = mx + b
		def hough_inv(th, rho):
			return (rho/np.sin(th), -1*np.cos(th)/np.sin(th))

		#line function mx + b
		def line(x, b, m):
			return (m*x + b)

		rhos = []
		thetas = []
		for i in range(len(rodDepths)):
			tt, rr = hough(rodDepths[i], rodTimes[i])
			for t in tt:
				thetas.append(t)
			for r in rr:
				rhos.append(r)


		H, thedges, rhoedges = np.histogram2d(thetas, rhos, bins=1500)
		#find the max x and y
		idx = list(H.flatten()).index(H.max())
		th_i, rho_i = idx / H.shape[1], idx % H.shape[1] 	#copied from stack exchange
		b, m = hough_inv(thedges[th_i + 1], rhoedges[rho_i + 1])


		fig, ax = plt.subplots(figsize=(10,7))
		ax.plot(rodDepths, rodTimes, 'ko')
		ax.plot(rodDepths, line(np.array(rodDepths), b, m), 'b-', label="treco - ttrue = ")
		plt.show()



	#an algorithm that histograms time of arrivals in 
	#the event, then forms to hypotheses and decides
	#two arrival times based on the average of points
	#contained in the fast component of the time distribution
	def algo_Simple(self, momentum, timesmear):	

		cutEvent = self.hadronicNoiseCut()
		if(len(cutEvent.hitPoints) <= 15):
			#print "did pass rough cut, adapt the cutting algorithm later"
			#print "you need to be able to make a general rough cut on noise"
			return (None, 1)

		#calculations of beta based on
		#a hypothesis of what the particle is
		mkaon = 0.493677 #GeV
		ekaon = np.sqrt(mkaon**2 + momentum**2)
		beta_kaon = momentum/ekaon

		T0_kaonHyp = []
		for hp in cutEvent.hitPoints:
			rho = np.sqrt(hp.getX()**2 + hp.getY()**2 + hp.getZ()**2)
			t = hp.getT()
			T0_kaonHyp.append(t - rho/(speedOfLight*beta_kaon))


		if(timesmear <= 0.007):
			pscut = 0.007
		else:
			pscut = timesmear

		passed = []
		earliest = min(T0_kaonHyp)
		for t0 in T0_kaonHyp:
			if(t0 <= pscut + earliest):
				passed.append(t0)

		if(len(passed) == 0):
			return (None, None)


		#return the average of those times
		#and the number of hits kept
		return (np.mean(passed), len(passed))


	#function that returns all T0 = t - r/v
	#quantities assuming kaon velocity
	def getKaonTime(self, momentum):
		#calculations of beta based on
		#a hypothesis of what the particle is
		mkaon = 0.493677 #GeV
		ekaon = np.sqrt(mkaon**2 + momentum**2)
		beta_kaon = momentum/ekaon

		T0_kaonHyp = []
		for hp in self.hitPoints:
			rho = np.sqrt(hp.getX()**2 + hp.getY()**2 + hp.getZ()**2)
			t = hp.getT()
			if(t > 20):
				continue
			T0_kaonHyp.append(t - rho/(speedOfLight*beta_kaon))

		return T0_kaonHyp

	#function that returns all T0 = t - r/v
	#quantities assuming pion velocity
	def getPionTime(self, momentum):
		#calculations of beta based on
		#a hypothesis of what the particle is
		mpion = 0.139570 #GeV
		epion = np.sqrt(mpion**2 + momentum**2)
		beta_pion = momentum/epion

		T0_pionHyp = []
		for hp in self.hitPoints:
			rho = np.sqrt(hp.getX()**2 + hp.getY()**2 + hp.getZ()**2)
			t = hp.getT()
			if(t > 20):
				continue
			T0_pionHyp.append(t - rho/(speedOfLight*beta_pion))

		return T0_pionHyp
