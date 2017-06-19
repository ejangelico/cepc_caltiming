# Contains all the hits in a given cylindrical layer in the calorimeter
class Layer:
	def __init__(self):
		self.d0   = None # The center of the layer
		self.dmin = None # The min depth that is still in the layer (closer to IP)
		self.dmax = None # The max depth that is still in the layer

		# Use average depth of all the hits as the depth of the layer!!! This can be d0

		self.hitPoints = [] # An array of the points that are in the layer

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

	# Returns the energy-weighted center position of all the hits in the layer
	def getCentroid(self):
		pass
