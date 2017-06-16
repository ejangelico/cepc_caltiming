import numpy as np
import sys


#these set of points defined
#with the origin at the center of the 
#barrel
class HitPoint:
	def __init__(self, x, y, z, t, E, coord=1):
		self.E = E	#GeV
		self.t = t	#ns

		#input coordinates given
		#in CARTESIAN
		if(coord == 1):
			self.x = x
			self.y = y
			self.z = z
			self.rho = np.sqrt(x*x + y*y)
			self.phi = np.arctan2(y, x)


		#initialization coordinates
		#given in CYLINDRICAL
		#(i.e. x ~ rho, y ~ phi, z = z)
		elif(coord == 0):
			self.rho = x
			self.phi = y
			self.z = z
			self.x = (self.rho)*np.cos(self.phi)
			self.y = (self.rho)*np.sin(self.phi)


	#---OPERATIONS---#

	def __mul__(self, b):
		return ((self.getX() * b.getX()) + (self.getY() * b.getY()) + (self.getZ() * b.getZ()))

	#scalar multiplication
	def scale(self, b):
		return Point(b*self.getX(), b*self.getY(), b*self.getZ(), self.t, self.E, 1)

	def __str__(self):
		return "(" + str(self.getX()) + ", " + str(self.getY()) + ", " + str(self.getZ()) + "), (t = " + str(self.getT()) + ", E = " + str(self.getE()) + ")\n"

	def getMag(self):
		return np.sqrt(self.getX()**2 + self.getY()**2 + self.getZ()**2)

	#get the distance between this point
	#and b, returns in mm
	def getDistance(self, b):
		return np.sqrt((self.getX() - b.getX())**2 + (self.getY() - b.getY())**2 + (self.getZ() - b.getZ())**2) 


	#---OPERATIONS END---#



	#---SETS AND GETS---#

	def getXYZ(self):
		return (self.x, self.y, self.z)

	def getRhoPhiZ(self):
		return (self.rho, self.phi, self.z)

	def getX(self):
		return self.x

	def getY(self):
		return self.y

	def getZ(self): 
		return self.z

	def getRho(self):
		return self.rho

	def getPhi(self):
		return self.phi

	def getE(self):
		return self.E

	def getT(self):
		return self.t

	def setE(self, e):
		self.E = e

	def setT(self, t):
		self.t = t

	#---SETS AND GETS END---#



