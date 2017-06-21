from rootpy.tree import Tree
from rootpy.io import root_open
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
import cPickle as pickle 
import Event
import DataSet
import HitPoint

#USAGE:
#python rootConverter.py <root file paths>
#
#the program will loop through all root files given in arguments
#so you can use, for example, "python rootConverter.py ./roots/*"



if __name__ == "__main__":


	for fn in sys.argv[1:]:
		f = root_open(fn, "read")
		print "On file " + fn

		#put trees into variables
		#syntax: file."TreeName"
		mcpTree = f.MCP
		evtTree = f.Evt
		hitTree = f.DigiHit
			

		#python objects that will be saved to output file
		eventList = []	#list of event class objects that store info for each event

		curEvent = None
		hitPoints = [] #list of HitPoints. See HitPoint class

		#tree has all of the events smushed into
		#one list, so I call "el" one of the list elements
		#and I split the smushed list into many lists indexed
		#by events
		for el in hitTree:
			#very first event only
			if(curEvent == None):
				curEvent = el.EventNr

			thisEvent = el.EventNr

			#we are now moving on to a new
			#event number, so store all of the 
			#active lists into an Event class object
			if(curEvent != thisEvent):
				eventList.append(Event.Event(hitPoints, curEvent))
				hitPoints = [] 

				#print to let you know that you're now on the next event
				curEvent = thisEvent
				print "On event " + str(thisEvent)

			

			#in any case, append this list element to our 
			#lists that are quantities of interest
			hitPoints.append(HitPoint.HitPoint(el.HitX, el.HitY, el.HitZ, el.Time, el.HitEn, 1))

		
		data = DataSet.DataSet(eventList)
		
		#save the file as a pickle file
		#the name of the saved file is the root file name
		#just with a ".p" instead of ".root"
		filename = fn[:-5]+".p"
		pickle.dump(data, open(filename, "wb"))

		#to load this file, do 
		#eventList = pickle.load(open(filname.p, 'rb'))

