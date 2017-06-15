from rootpy.tree import Tree
from rootpy.io import root_open
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import cPickle as pickle 
"""
global hbar
global cc
hbar = 6.582119514e-16 #ev - s
cc = 0.299792458 #m/ns


#input is VACUUM WAVELENGTH in nm
#this is the Sellmeier equation
#with parameters taken from Daimon and Masumura
#for HPLC DI water (fancy 16-18MohM)
#at 21.5C
def getIndex(lam):
	A1 = 5.689083732e-1
	A2 = 1.719708856e-1
	A3 = 2.062501582e-2
	A4 = 1.123965424e-1
	l12 = 5.110301794e-3
	l22 = 1.825180155e-2
	l32 = 2.624158904e-2
	l42 = 1.067505178e1

	lam = lam*1e-3 	#put lambda in microns

	n2 = 1 + (A1*lam*lam)/(lam*lam - l12) \
	+ (A2*lam*lam)/(lam*lam - l22) \
	+ (A3*lam*lam)/(lam*lam - l32) \
	+ (A4*lam*lam)/(lam*lam - l42)

	return np.sqrt(n2)

#the derivative of above 
def dndl(lam):
	A1 = 5.689083732e-1
	A2 = 1.719708856e-1
	A3 = 2.062501582e-2
	A4 = 1.123965424e-1
	l12 = 5.110301794e-3
	l22 = 1.825180155e-2
	l32 = 2.624158904e-2
	l42 = 1.067505178e1
	nlam = getIndex(lam) #lam still in nm

	lam = lam*1e-3 	#put lambda in microns

	bigfactor = (-2*A1*l12*lam)/((lam*lam - l12)**2) +\
	(-2*A2*l22*lam)/((lam*lam - l22)**2) +\
	(-2*A3*l32*lam)/((lam*lam - l32)**2) +\
	(-2*A4*l42*lam)/((lam*lam - l42)**2)

	dndl = 0.5*nlam*bigfactor
	return dndl*1e-3 #back to units of inverse nm




#finds the group velocity
#for a given photon energy
def getVelocity(E):
	omega = (E/hbar)*1e-9 #1/ns
	lam_vacuum = (2*np.pi*cc/omega)*1e9 #nm
	n = getIndex(lam_vacuum)
	dn = dndl(lam_vacuum)
	vg = cc/(n - lam_vacuum*dn)
	return vg
	

#get energy given vacuum wavelength
#in nm
def getEnergy(lam):
	return (2*np.pi*cc/(lam*1e-9))*hbar*1e9


def velocityPlot():
	lam = np.linspace(200, 1000, 100)
	vs = []
	for l in lam:
		omega = 2*np.pi*cc/(l*1e-9) #1/ns
		omega = omega*1e9 #1/s
		E = hbar*omega
		vs.append(getVelocity(E))


	fig, ax = plt.subplots(figsize=(25, 10))
	ax.plot(lam, vs, 'g', linewidth=3)
	ax.set_xlabel("vacuum wavelength (nm)", fontsize=26)
	ax.set_ylabel("group velocity (m/ns)", fontsize=26)
	ax.set_title("Group velocity in HPLC water", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	plt.savefig("plots/physics/groupv.png", bbox_inches='tight')
	#plt.show()
	sys.exit()

#plots both group and phase index
def plotIndex():
	wl = np.linspace(200, 1000, 100)
	gns = []
	pns = []
	for w in wl:
		gns.append(getIndex(w) - w*dndl(w))
		pns.append(getIndex(w))

	fig, ax = plt.subplots(figsize=(25, 10))
	ax.plot(wl, gns, 'r--', linewidth=3, label="Group index")
	ax.plot(wl, pns, 'r-', linewidth=3, label="Phase index")
	ax.set_xlabel("vacuum wavelength (nm)", fontsize=26)
	ax.set_ylabel("refractive index", fontsize=26)
	ax.set_title("Group and phase refractive index in HPLC DI water", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	ax.legend(fontsize=26)
	plt.savefig("plots/physics/refIndex.png", bbox_inches='tight')
	#plt.show()
	sys.exit()


#plots time of arrival differences
#for different wavelengths, relative
#to 400nm
def plotTimeSpread():
	fig, ax = plt.subplots(figsize=(10, 10))
	lam_vac = np.arange(200, 900, 50)
	pl = np.linspace(0, 10, 1000) #meters of path length
	v400 = getVelocity(getEnergy(400))
	for lam in lam_vac:
		timeDiff = []
		v = getVelocity(getEnergy(lam))
		for d in pl:
			timeDiff.append((d*v - d*v400)/(v*v400))

		rgb = wlrgb.wavelength_to_rgb(lam)
		ax.plot(pl, timeDiff, color=rgb, label=r"$\lambda = $" + str(lam))

	ax.legend()
	ax.set_xlabel("path length (m)", fontsize=26)
	ax.set_ylabel("time difference \nrelative to 400nm (ns)", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	majorLocatorY = MultipleLocator(1)
	minorLocatorY = MultipleLocator(0.25)
	majorLocatorX = MultipleLocator(1)
	minorLocatorX = MultipleLocator(0.5)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.set_ylim([-4, 4])
	ax.grid(color='r', linestyle='--', which='major')
	ax.grid(color='r', linestyle='--', which='minor', alpha=0.5)

	#plt.savefig("plots/physics/timespread.png", bbox_inches='tight')
	plt.show()






if __name__ == "__main__":

	for fn in sys.argv[1:]:
		f = root_open(fn, "read")
		print "On file " + fn
		if(int(fn.split('_')[1].split('.')[0][0]) < 3f):
			continue
		phtree = f.PhTree
		eltree = f.ElTree
		rays = []
		ens = []
		evcount = 0
		tref = None 	#time reference to t=0
		for event in phtree:
			print "On event " + str(evcount)
			for i in range(len(event.X)):
				pinit = Point.Point(event.X[i], event.Y[i], event.Z[i], 1)
				vinit = Point.Point(event.Px[i], event.Py[i], event.Pz[i], 1)
				energy = event.Energy[i]
				tinit = event.Tg[i]
				if(energy == 0):
					continue
				if(tref == None):
					tref = tinit
				elif(tref > tinit):
					#find first photon and use that as t=0 reference
					tref = tinit
				else:
					omeg = (energy*1e6/hbar)*1e-9
					mat_lambda = (2*np.pi*getVelocity(energy*1e6)/omeg)*1e9
					rays.append(Ray.Ray(pinit.scale(1e-3), vinit, getVelocity(energy*1e6), tinit, mat_lambda))


			for r in rays:
				r.setStartTime(r.getStartTime() - tref)
			pickle.dump(rays, open(fn[:-5] + ".p", "wb"))



	
"""
