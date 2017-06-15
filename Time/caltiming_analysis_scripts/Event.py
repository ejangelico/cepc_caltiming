import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys

class Event:
	def __init__(self, hitPoints=None, hitEn=None, evNum=None):
		self.hitPoints = hitPoints  #array of hit points corresponding to pixel positions in the ECAL
		self.hitEn = hitEn  #array of energies that lines up with the hit points	
		self.evNum = evNum 	#integer id for the event