
# Contains all the hits in a given cylindrical layer in the calorimeter
class Layer:
	def __init__(self):
		self.d0   = None # The center of the layer
		self.dmin = None # The min depth that is still in the layer (closer to IP)
		self.dmax = None # The max depth that is still in the layer

		self.hitPoints = [] # An array of the points that are in the layer

	# Given a point, will initialize the layer
	def intializeWithPoint(self, HitPoint, width = None, (dmin, dmax) = (None, None)):
		if (width is None) and (dmin is None):
			print "Must initialize layer size"
			return

	# Test whether the point is in the layer
	def containsPoint(self, HitPoint):
		pass #return True or False

	# Tests if a point is in the layer, and if it is, add it
	def addPointe(self, HitPoint):
		pass

	# Returns the earliest time in the layer
	def getFirstTime(self):
		pass

	# Returns the energy-weighted center position of all the hits in the layer
	def getCentroid(self):
		pass
