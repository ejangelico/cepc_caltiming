import Point

# Contains all the hits in a given cylindrical layer in the calorimeter
class Layer:
	def __init__(self):
		self.d0   = None # The center of the layer
		self.dmin = None # The min depth that is still in the layer (closer to IP)
		self.dmax = None # The max depth that is still in the layer

		# Use average depth of all the hits as the depth of the layer!!! This can be d0

		self.hitPoints = [] # An array of the points that are in the layer


	def getCenter(self):
		return self.d0

	# Given a point, will initialize the layer
	def initializeWithPoint(self, hitPoint, width = None, (dmin, dmax) = (None, None)):
		if (width is None) and (dmin is None):
			print "Must initialize layer size"
			return
		self.d0 = hitPoint.getRho()
		if width is not None:
			self.dmin = self.d0-width/2.0
			self.dmax = self.d0+width/2.0
		else:
			self.dmin = dmin
			self.dmax = dmax
		self.hitPoints.append(hitPoint)

	#initialize a layer given a center rho coordinate
	def initializeWithRadius(self, rad, width):
		if(rad == None or width == None):
			return
		self.d0 = rad
		self.dmin = self.d0-width/2.0
		self.dmax = self.d0+width/2.0
		return


	# Test whether the point is in the layer
	def containsPoint(self, hitPoint):
		if self.dmin <= hitPoint.getRho() < self.dmax:
			return True
		else:
			return False

	# Tests if a point is in the layer, and if it is, add it
	def addPoint(self, hitPoint):
		if self.containsPoint(hitPoint):
			self.hitPoints.append(hitPoint)
			return True
		else:
			return False

	# Returns the earliest time in the layer
	def getFirstTime(self):
		tList = []
		for hitPoint in self.hitPoints:
			tList.append(hitPoint.getT())
		return min(tList)

	# Returns the earliest point in the layer
	def getFirstPoint(self):
		if(self.hitPoints == None):
			return None
		if(len(self.hitPoints) == 0):
			return None

		hpBest = self.hitPoints[0]
		tBest = hpBest.getT()
		for hp in self.hitPoints:
			if(hp.getT() < tBest):
				hpBest = hp
				tBest = hp.getT()

		return hpBest

	# removes hit points in the layer that do not
	# pass a set of cuts. Passing is constituted by 
	# the point lying in the range, for example trange
	# that is a two element list, trange = [floor, ceiling]
	def cutHitPoints(self, rhorange=None, phirange=None, zrange=None, trange=None):
		removes = [] #list of hitpoints to remove
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

		#remove the points that didn't pass
		for rm in removes:
			self.hitPoints.remove(rm)

		return


	# Returns the plain geometric
	# centroid of the hit points in the layer
	def getSpaceTimeCentroid(self):
		if(self.hitPoints == None):
			return None
		if(len(self.hitPoints) == 0):
			return None

		#define points from the hitpoints
		points = []
		for hp in self.hitPoints:
			points.append(Point.Point(hp.getX(), hp.getY(), hp.getZ(), 1))

		#geometric 3D centroid
		spaceCentroid = Point.Point(0,0,0,1)
		for p in points:
			spaceCentroid = spaceCentroid + p

		spaceCentroid = spaceCentroid.scale(float(1.0/len(points)))
		#time centroid is just average
		timeCentroid = 0
		for hp in self.hitPoints:
			timeCentroid += hp.getT()

		timeCentroid = timeCentroid/float(len(self.hitPoints))
		#put it together into a point
		spaceCentroid.setT(timeCentroid)
		return spaceCentroid



