import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import cPickle as pickle 
import time

import DataSet
import Event

if __name__ == "__main__":

	print "Loading file...", 
	sys.stdout.flush()
	data = pickle.load(open("../pickles/pions/noB/5Gev_1000evts.p", 'rb'))
	print "Done."
	sys.stdout.flush()

	#data.events[10].timeVsDepth(True)
	data.events[10].energyDisplay(False)
	#data.events[12].algo_linearFirstTimeByLayer(plotting=True)
	#smdata = data.smear(0.01, 0)
	#data.timeReco(0, False)
	#data.timeRecoSmearing(np.linspace(0, 0.5, 30))