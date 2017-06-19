import DataSet
import Helper
import numpy as np
import cPickle as pickle
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import sys

mg  = 0
me  = 511.00
mPi = 139.57
mK  = 493.68
mPi0= 134.98
mK0 = 497.65
mp  = 938.27
mn  = 939.57

def tof(p, m, q):
	if q == 0:
		return 6.004/1000.0 * np.sqrt(p**2 + m**2)/p
	else:
		return 3.177/1000.0 * np.sqrt(p**2 + m**2) *np.arccos(1 - 0.5 * (1890.0/p)**2)

def normal(x, mu, sigma):
	return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-mu)**2/(2*sigma**2))

def overlap(x, mu1, sig1, mu2, sig2):
	return min([normal(x, mu1, sig1), normal(x, mu2, sig2)])

def plotGaussianOverlap():
	p = 4000.0
	tofPi = tof(p, mPi, q = 1)
	tofK  = tof(p, mK,  q = 1)
	deltaT = 0.0255 # Corresponding to dtpix = 50 ps for 1 GeV electrons

	print
	print
	print "tofPi:", tofPi
	print "tofK: ", tofK
	print "deltaT:", deltaT

	fig, ax = plt.subplots()
	x = np.linspace(tofPi - 5*deltaT, tofPi + 5*deltaT, 1000)
	y = [normal(i, tofPi, deltaT) for i in x]
	ax.plot(x, y, color = 'k', linewidth = 2.0)

	x = np.linspace(tofK - 5*deltaT, tofK + 5*deltaT, 1000)
	y = [normal(i, tofK, deltaT) for i in x]
	ax.plot(x, y, color = 'k', linewidth = 2.0)

	av = 0.5*(tofK + tofPi)
	x = np.linspace(av - 8*deltaT, av + 8*deltaT, 1000)
	y = [overlap(i, tofPi, deltaT, tofK, deltaT) for i in x]
	plt.fill(x, y, fill = False, hatch = "\\")
	plt.plot(x, y, 'r', linewidth = 2.0)
	plt.title("Pixel Time Resolution: 50 ps; Momentum: 4 GeV", fontsize = 20)
	plt.xlabel("Time of Arrival (ns)")
	plt.ylabel("Prob. Density")
	Helper.resize(fig, ax)
	ax.xaxis.set_major_locator(MultipleLocator(0.05))
	ax.xaxis.set_minor_locator(MultipleLocator(0.01))
	ax.yaxis.set_major_locator(MultipleLocator(2))
	ax.yaxis.set_minor_locator(MultipleLocator(0.5))
	plt.show()

print "Loading file...", 
sys.stdout.flush()
data = pickle.load(open("../pickles/electrons/AnaHit_electron_10GeV_3760.p", 'rb'))
print "Done."
sys.stdout.flush()

deltaTPixList = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
deltaTList = data.timeRecoSmearing(deltaTPixList)

print "Doing integrations."
sys.stdout.flush()
pList = np.linspace(1000, 10000, 100)
pListGeV = [p/1000.0 for p in pList]
overlapArea = []
fig, ax = plt.subplots()
for deltaT in deltaTList:
	overlapArea.append([])
	for p in pList:
		tofPi = tof(p, mPi, q = 1)	
		tofK  = tof(p, mK,  q = 1)
		av = 0.5*(tofPi + tofK)
		overlapArea[-1].append(scipy.integrate.quad(overlap, av-8*deltaT, av+8*deltaT, args = (tofPi, deltaT, tofK, deltaT))[0])
	ax.plot(pListGeV, overlapArea[-1], color = 'k')
plt.xlabel("Momentum (GeV)")
plt.ylabel("Overlap Coefficient")
Helper.resize(fig, ax)
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
plt.show()
