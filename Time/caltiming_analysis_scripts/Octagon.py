import numpy as np
import Point
import matplotlib.pyplot as plt

class Octagon:
	def __init__(self, R):
		self.R = R # shortest distance between origin and the edge
		self.corners = None

	def getCorners(self):
		if self.corners is not None:
			return self.corners
		corners = []
		outerR = self.R/np.cos(np.pi/8.0)
		for i in range(0, 8):
			x = outerR*np.cos(np.pi/8.0 * (2*i+1))
			y = outerR*np.sin(np.pi/8.0 * (2*i+1))
			corners.append([x, y])
		self.corners = corners
		return self.corners

	def plot(self, ax = None):
		self.getCorners()
		N = range(-1, 7)
		if ax is None:
			for i in N:
				x = [self.corners[i][0], self.corners[i+1][0]]
				y = [self.corners[i][1], self.corners[i+1][1]]
				plt.plot(x, y, 'k')
			plt.show()
		else:
			for i in N:
				x = [self.corners[i][0], self.corners[i+1][0]]
				y = [self.corners[i][1], self.corners[i+1][1]]
				ax.plot(x, y, 'k')
			return ax
	
	# Given the angle, gets the radius of the octagon at that point
	# t = theta
	def getRadius(self, t):
		while t > 2*np.pi:
			t -= 2*np.pi
		while t < 0:
			t += 2*np.pi
	
		while t > np.pi/8:
			t -= np.pi/4
		return self.R/np.cos(t)
		
	# Checks whether a point [x0, y0] is in the octagon
	def isInside(self, x0, y0):
		# Find theta and r of the point. Find radius of octagon at that angle
		if self.getRadius(np.arctan(y0/x0)) > np.sqrt(x0**2+y0**2):
			return True
		else:
			return False
		
	# Given the radius of a circle, find intersection of circle and octagon
	# Assumes circle is centered at x = -rho, y = 0
	# Returns first point of intersect, starting from origin
	# Return None if no intersection
	# Returns [x0, y0] intersection point otherwise
	# tol = tolerance in units of mm
	def circleIntersect(self, rho, tol = 1e-4):
		thetaMin = 0.01
		thetaMax = np.pi-0.01
		foundIntersect = False
		while rho*(thetaMax-thetaMin) > tol:
			thetaList = np.linspace(thetaMin, thetaMax, 100)
			for i in range(0, len(thetaList)):
				t = thetaList[i]
				r = rho*np.sqrt(2*(1-np.cos(t)))
				phi = np.arctan(np.sin(t)/(np.cos(t)-1))
				if r > self.getRadius(phi):
					thetaMin = thetaList[i-1]
					thetaMax = thetaList[i]
					foundIntersect = True
					break
			if not foundIntersect:
				return None
		thetaInt = 0.5*(thetaMax+thetaMin)
		x0 = rho*(np.cos(thetaInt)-1)
		y0 = rho*np.sin(thetaInt)
		return [x0, y0]
