import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import scipy.optimize as opt
import Point


def resize(fig, ax):
	ax.yaxis.label.set_size(21)
	ax.xaxis.label.set_size(21)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=20, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	return (fig, ax)


# Given x-values and y-values of a curve with one peak, computes the full width at half max, in units of x
def FWHM(xVals, yVals):
		xVals = xVals.tolist()
		yVals = yVals.tolist()

		maxIndex = yVals.index(max(yVals))
		peakCount      = yVals[maxIndex]
		peakCountx  = xVals[maxIndex]

		for i in range(maxIndex, len(xVals)):
			if yVals[i] < peakCount/2.0:
				break
		x0 = xVals[i-1]
		y0 = yVals[i-1]
		x1 = xVals[i]
		y1 = yVals[i]
		tHalfUp = x0 + (0.5 *peakCount - y0) *(x1-x0)/(y1-y0)

		for i in range(0, maxIndex):
			if yVals[i] > peakCount/2.0:
				break
		x0 = xVals[i-1]
		y0 = yVals[i-1]
		x1 = xVals[i]
		y1 = yVals[i]
		tHalfDown = x0 + (0.5 *peakCount - y0) *(x1-x0)/(y1-y0)

		return tHalfUp-tHalfDown


# LineParams = [x0, y0, z0, theta, phi] <= 5 params of line in 3D
# Params = [w0, [[x1, y1, z1, E1], [x2, y2, z2, E2], ... [xN, yN, zN, EN]]] list of all points to fit to
# w0 = weight function parameter
def LinChiSq3DWeighted(LineParams, Params):
	x1 = np.array(LineParams[:3])
	theta = LineParams[3]
	phi = LineParams[4]

	w0 = Params[0]
	points = Params[1]

	x2 = x1 + np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
	Etot = np.sum(points, axis=0)[3]

	D = 0
	for p in points:
		w = w0 + np.log(p[3]/Etot)
		if w < 0: 
			continue
		x0 = np.array(p[:3])
		d = np.linalg.norm(np.cross(x0-x1, x0-x2))
		D += (d*w)**2
	return D

# points = [[x1, y1, z1, E1], [x2, y2, z2, E2], ... [xN, yN, zN, EN]]
# w0 = Energy weighting function parameter
# Returns a list of two points: [LinePoint, LineVector]
def LinFit3D(points, w0):
	points = np.array(points)
	x0, y0, z0, E0 = points[0]
	theta = np.pi/2.0
	phi = np.pi/2.0
	output = opt.minimize(LinChiSq3DWeighted, [x0, y0, z0, theta, phi], args = [w0, points], bounds = [[None, None], [None, None], [None, None], [0, np.pi], [0, 2*np.pi]])
	x0fit, y0fit, z0fit, thetafit, phifit = output.x
	mxfit, myfit, mzfit = np.cos(phifit)*np.sin(thetafit), np.sin(phifit)*np.sin(thetafit), np.cos(thetafit)
	
	return [Point.Point(x0fit, y0fit, z0fit, None), Point.Point(mxfit, myfit, mzfit, None)]
