import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def resize(fig, ax):
	fig.set_size_inches(25, 15)
	ax.yaxis.label.set_size(20)
	ax.xaxis.label.set_size(20)
	ax.get_xaxis().set_tick_params(labelsize=20, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=20, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	return (fig, ax)


