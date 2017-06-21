import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
		peakCount  = yVals[maxIndex]
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

#get the first hit in a set of hitpoints
def getFirstHit(hitpoints):
	tsort = sorted(hitpoints, key=lambda x: x.getT())
	return tsort[0]


#return the two intersection points of an infinite
#cylinder with a line, assuming the cylinder has 
#its principle axis as the z axis and "line"
#is [Point, Point] , i.e. [shift point, vector]
def getCylinderIntersection(radius, line):
	#At^2 + Bt + C = 0 from mrl.nyu.edu/~dzorin/rend05/lecture2.pdf
	cylinderAxis = Point.Point(0,0,1,1) #z axis
	v = line[1]
	x = line[0]
	A = (v - cylinderAxis.scale(v*cylinderAxis))*(v - cylinderAxis.scale(v*cylinderAxis))
	B = 2*(v - cylinderAxis.scale(v*cylinderAxis))*(x - cylinderAxis.scale(x*cylinderAxis))
	C = (x - cylinderAxis.scale(x*cylinderAxis))*(x - cylinderAxis.scale(x*cylinderAxis)) - radius**2
	#two solution from quadratic equation
	t_plus = (1.0/(2*A))*(-1*B + np.sqrt(B*B - 4*A*C))
	t_minus = (1.0/(2*A))*(-1*B - np.sqrt(B*B - 4*A*C))
	point_plus = Point.Point(t_plus*v.getX() + x.getX(), t_plus*v.getY() + x.getY(), t_plus*v.getZ() + x.getZ(), 1)
	point_minus = Point.Point(t_minus*v.getX() + x.getX(), t_minus*v.getY() + x.getY(), t_minus*v.getZ() + x.getZ(), 1)
	return (point_plus, point_minus)


