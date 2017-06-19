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
	data = pickle.load(open("../pickles/electrons/AnaHit_electron_10GeV_3760.p", 'rb'))
	print "Done."
	sys.stdout.flush()

	data.timeRecoSmearing(np.linspace(0, 0.1, 15))
