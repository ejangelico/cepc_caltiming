/*
 * used to readout the FD, FD_10First, FD_20Later, ... and other quantities of bush
 * objective: pattern tagging for ECAL as: penetrating mip, deep interaction CH, EM Cluster; 
 * 				  HCAL as: penetrating mip, early interaction H, later interaction H, etc.
 *	//Used as SP Analysis Code
 * */

#include <TimeAna.hh>
#include "ArborTool.hh"
#include "ArborToolLCIO.hh"
#include "DetectorPos.hh"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <UTIL/CellIDDecoder.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TVector3.h>

using namespace std;

const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";

TimeAna aTimeAna ;
TimeAna::TimeAna()
	: Processor("TimeAna"),
	_output(0)
{
	_description = "Measure Bush Quantities" ;

	_treeFileName="TimeAna.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);


	_treeName="Neutral";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);

	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);


}

void TimeAna::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB


	_outputTree->Branch("EventNr", 	&_eventNr, 	"EventNr/I");
	_outputTree->Branch("NClu", 	&_NClu, 	"NClu/I");

	_outputTree->Branch("NeuAveTime", 	&_NeuAveTime, 	"NeuAveTime/F");
	_outputTree->Branch("NeuPeakTime", 	&_NeuPeakTime, 	"NeuPeakTime/F");
	_outputTree->Branch("NeuNCount", 	&_NeuNCount, 	"NeuNCount/F");
	_outputTree->Branch("NeuAngle", 	&_NeuAngle, 	"NeuAngle/F");
	_outputTree->Branch("NeuEnergy", 	&_NeuEnergy, 	"NeuEnergy/F");
	_outputTree->Branch("NeuPos", 	&_NeuPos, 	"NeuPos/F");
	_outputTree->Branch("NeuNHits", 	&_NeuNHits, 	"NeuNHits/F");
	
	_outputTree->Branch("FragAveTime", 	&_FragAveTime, 	"FragAveTime/F");
	_outputTree->Branch("FragPeakTime", 	&_FragPeakTime, "FragPeakTime/F");
	_outputTree->Branch("FragNCount", 	&_FragNCount, 	"FragNCount/F");
	_outputTree->Branch("FragAngle", 	&_FragAngle, 	"FragAngle/F");
	_outputTree->Branch("FragEnergy", 	&_FragEnergy, 	"FragEnergy/F");
	
	}

void TimeAna::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		_eventNr=evtP->getEventNumber();
		//if( _eventNr % 100 == 0)
			std::cout<<_eventNr<<" evts processed"<<std::endl;

		_NeuAveTime=0;
		_NeuPeakTime=0;
		_NeuAngle=0;
		_NeuNCount=0;
		_NeuEnergy=0;
		_FragAveTime=0;
		_FragPeakTime=0;
		_FragAngle=0;
		_FragNCount=0;
		_FragEnergy=0;
		std::vector<Cluster*> mergedcluster;
		mergedcluster.clear();

		TVector3 NeuCluTimeVec;
		TVector3 FragCluTimeVec;

		try{
			LCCollection * colNeutral = evtP->getCollection("ArborNeutral");
			LCCollection * colNeutral1 = evtP->getCollection("ClusterNeutralCore");

			_NClu = colNeutral->getNumberOfElements();
			cout<<"cluster number "<<_NClu<<" : "<<colNeutral1->getNumberOfElements()<<endl;
			if(_NClu>0){
			Cluster * Neutral = dynamic_cast<Cluster*>(colNeutral->getElementAt(0));
			_NeuAveTime=ClusterTime(Neutral)[2];
			_NeuPeakTime=ClusterTime(Neutral)[3];
			_NeuNCount=ClusterTime(Neutral)[4];
			_NeuAngle=ClusterTime(Neutral)[5];
			_NeuEnergy=Neutral->getEnergy();
			_NeuNHits=Neutral->getCalorimeterHits().size();
			TVector3 cluPos=Neutral->getPosition();
			_NeuPos=cluPos.Mag();
			cout<<"peak "<<_NeuPeakTime<<" ncount "<<_NeuNCount<<endl;
			for(int i=1;i<_NClu;i++){

				Cluster * a_Frag = dynamic_cast<Cluster*>(colNeutral->getElementAt(i));
				mergedcluster.push_back(a_Frag);
			}
			ClusterImpl * Frag = NaiveMergeClu(mergedcluster);

			_FragAveTime=ClusterTime(Frag)[0];
			_FragPeakTime=ClusterTime(Frag)[1];
			_FragNCount=ClusterTime(Frag)[2];
			_FragAngle=ClusterTime(Frag)[3];
			_FragEnergy=Frag->getEnergy();
			_outputTree->Fill();
			}
		}
		catch (lcio::DataNotAvailableException err) { }

	}  	

}	

void TimeAna::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


