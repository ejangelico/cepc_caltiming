import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


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

