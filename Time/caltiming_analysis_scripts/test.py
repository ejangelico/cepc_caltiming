import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import cPickle as pickle 

import DataSet
import Event

if __name__ == "__main__":
	data = pickle.load(open("../pickles/electrons/AnaHit_electron_1GeV_2501.p", 'rb'))
	print data.events[124].evNum