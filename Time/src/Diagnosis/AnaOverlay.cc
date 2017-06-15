#include <AnaOverlay.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
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
#include <TMath.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TVector3.h>

using namespace std;

AnaOverlay aAnaOverlay ;
AnaOverlay::AnaOverlay()
	: Processor("AnaOverlay"),
	_output(0)
{
	_description = "CLuster Analysis" ;

	_treeFileName="AnaOverlay.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);
	_Merge=1;
	registerProcessorParameter( "Merge" , 
			"If zero No Merge" ,
			_Merge ,
			_Merge);
}

void AnaOverlay::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree("Evt","Evt");
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB


	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("Num", &num,"Num/I");
	_outputTree->Branch("NClu", &_nclu,"NClu/I");
	_outputTree->Branch("THEn",&_THEn,"THEn/F");
	_outputTree->Branch("TCEn",&_TCEn,"TCEn/F");
	_outputTree->Branch("LCEn",&_LCEn,"LCEn/F");
	_outputTree->Branch("SCEn",&_SCEn,"SCEn/F");
	_outputTree->Branch("LCEnRSCEn",&_LCEnRSCEn,"LCEnRSCEn/F");
	_outputTree->Branch("LCEnPSCEn",&_LCEnPSCEn,"LCEnPSCEn/F");

	num = 0;
}

void AnaOverlay::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{ 

			_eventNr=evtP->getEventNumber();
			cout << "Event:" << _eventNr << endl;
			LCCollection * MCP = evtP->getCollection("MCParticle");
			_nclu = 0;	


			LCCollection * col_Clu;
			if(_Merge==0)  col_Clu = evtP->getCollection("EHBushes");
			//if(_Merge==1)  col_Clu = evtP->getCollection("ArborNeutral");
			if(_Merge==1)  col_Clu = evtP->getCollection("CluAB_1st");

			_THEn = 0;
			_TCEn = 0;
			_LCEn = 0;


			for(int c0 = 0; c0 < col_Clu->getNumberOfElements(); c0++)
			{
				Cluster *a_clu = dynamic_cast<EVENT::Cluster *>(col_Clu->getElementAt(c0));
				_nclu++;
				float aCluEn = a_clu->getEnergy();
				_TCEn += aCluEn;

				if(aCluEn > _SCEn) _SCEn = aCluEn;
				if(aCluEn > _LCEn){
					_SCEn = _LCEn;
					_LCEn = aCluEn;
				}

			}
			_LCEnRSCEn=_LCEn/_SCEn;
			_LCEnPSCEn=_LCEn+_SCEn;


			//----------------------------DigiHit------------------------------
			LCCollection * col = evtP->getCollection("DigiSiHit");
			int numElements = col->getNumberOfElements();

			if (col->getTypeName() == LCIO::CALORIMETERHIT)
			{

				for (int j(0); j < numElements; ++j) 
				{
					CalorimeterHit *a_DigiHit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;

					_HitE=a_DigiHit->getEnergy();
					_THEn += _HitE;			

				}
			}

			//--------------------------------------------------------------------

			_outputTree->Fill();
			num++;
		}catch (lcio::DataNotAvailableException err) { }
	}
}

void AnaOverlay::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


