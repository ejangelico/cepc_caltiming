#include <AnaHit.hh>
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

#include <DetectorPos.hh>
#include <ArborTool.hh>
#include <ArborToolLCIO.hh>

using namespace std;

AnaHit aAnaHit ;
AnaHit::AnaHit()
	: Processor("AnaHit"),
	_output(0)
{
	_description = "CLuster Analysis" ;

	_treeFileName="AnaHit.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_overwrite=1;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);
}

void AnaHit::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputMCP = new TTree("MCP","MCP");
	_outputMCP->SetAutoSave(32*1024*1024);  // autosave every 32MB

	_outputTree = new TTree("Evt","Evt");
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB

	_outputDigiHit = new TTree("DigiHit", "DigiHit" );
	_outputDigiHit->SetAutoSave(32*1024*1024);  // autosave every 32MB

	_outputMCP->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputMCP->Branch("PDG", &_PDG, "PDG/I");
	_outputMCP->Branch("ParentPDG", &_ParentPDG, "ParentPDG/I");
	_outputMCP->Branch("Energy", &_MCPEn, "Energy/F");
	_outputMCP->Branch("Vertex", &_MCPVertex, "Vertex[3]/F");
	_outputMCP->Branch("Endpoint", &_MCPEndpoint, "Endpoint[3]/F");
	_outputMCP->Branch("MCPP", &_MCPP, "MCPP[3]/F");

	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("THEn",&_THEn,"THEn/F");

	_outputDigiHit->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputDigiHit->Branch("HitX",&_HitPosX,"HitX/F");
	_outputDigiHit->Branch("HitY",&_HitPosY,"HitY/F");
	_outputDigiHit->Branch("HitZ",&_HitPosZ,"HitZ/F");
	_outputDigiHit->Branch("HitEn",&_HitE,"HitEn/F");
	_outputDigiHit->Branch("Time",&_Time,"Time/F");
	_outputDigiHit->Branch("ID0",&_ID0,"ID0/I");
	_outputDigiHit->Branch("ID1",&_ID1,"ID1/I");
	_outputDigiHit->Branch("M",&_M,"M/I");
	_outputDigiHit->Branch("S",&_S,"S/I");
	_outputDigiHit->Branch("I",&_I,"I/I");
	_outputDigiHit->Branch("J",&_J,"J/I");
	_outputDigiHit->Branch("K",&_K,"K/I");
	_outputDigiHit->Branch("MCPID", &_MCPID, "MCPID/I");

	num = 0;
}

void AnaHit::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{ 
			_eventNr=evtP->getEventNumber();
			cout << "Event:" << _eventNr << endl;
			LCCollection * MCP = evtP->getCollection("MCParticle");
			for(int i0 = 0; i0 < MCP->getNumberOfElements(); i0++)
			{
				MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCP->getElementAt(i0));
				_PDG = a_MCP->getPDG();
				_MCPEn = a_MCP->getEnergy();                            
				_MCPP[0] = a_MCP->getMomentum()[0];
				_MCPP[1] = a_MCP->getMomentum()[1];
				_MCPP[2] = a_MCP->getMomentum()[2];

				_MCPVertex[0] = a_MCP->getVertex()[0];
				_MCPVertex[1] = a_MCP->getVertex()[1];
				_MCPVertex[2] = a_MCP->getVertex()[2];

				_MCPEndpoint[0] = a_MCP->getEndpoint()[0];
				_MCPEndpoint[1] = a_MCP->getEndpoint()[1];
				_MCPEndpoint[2] = a_MCP->getEndpoint()[2];

				if( a_MCP->getParents().size() == 0 ){
					_ParentPDG=0;
				}
				else if(a_MCP->getParents().size() != 0){
					_ParentPDG = a_MCP->getParent(0)->getPDG();
				}
				_outputMCP->Fill();
			}



			//LCCollection * col = evtP->getCollection("DigiSiHit");
			LCCollection * col = evtP->getCollection("SiCalCollection");
			int numElements = col->getNumberOfElements();

			_THEn = 0;
			if (col->getTypeName() == LCIO::CALORIMETERHIT)
			{
			const std::string ECALCellIDDecoder = "S-1:8,M:8,:I:16,K-1:12,GRZone:4,J:16";
			CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);
				for (int j(0); j < numElements; ++j) 
				{
					CalorimeterHit *a_DigiHit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;

					_HitPosX=a_DigiHit->getPosition()[0];
					_HitPosY=a_DigiHit->getPosition()[1];
					_HitPosZ=a_DigiHit->getPosition()[2];
					_HitE=a_DigiHit->getEnergy();
					_THEn += _HitE;			
					//HitEnError = a_DigiHit->getEnergyError();
					_ID0=a_DigiHit->getCellID0();
					_ID1=a_DigiHit->getCellID1();
					_Time=a_DigiHit->getTime();
					_MCPID = -999;

					_M = idDecoder(a_DigiHit)["M"];
					_S = idDecoder(a_DigiHit)["S-1"];
					_I = idDecoder(a_DigiHit)["I"];
					_J = idDecoder(a_DigiHit)["J"];
					_K = idDecoder(a_DigiHit)["K-1"];

					_outputDigiHit->Fill();
				}
			}

			else if (col->getTypeName() == LCIO::SIMCALORIMETERHIT)
			{
			const std::string ECALCellIDDecoder = "S-1:8,M:8,:I:16,K-1:12,GRZone:4,J:16";
			CellIDDecoder<SimCalorimeterHit> idDecoder(ECALCellIDDecoder);
				for (int j(0); j < numElements; ++j) 
				{
					SimCalorimeterHit *a_DigiHit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;

					_HitPosX=a_DigiHit->getPosition()[0];
					_HitPosY=a_DigiHit->getPosition()[1];
					_HitPosZ=a_DigiHit->getPosition()[2];
					_HitE=a_DigiHit->getEnergy();
					_THEn += _HitE;			
					//HitEnError = a_DigiHit->getEnergyError();
					_ID0=a_DigiHit->getCellID0();
					_ID1=a_DigiHit->getCellID1();
					_MCPID = -999;

					_M = idDecoder(a_DigiHit)["M"];
					_S = idDecoder(a_DigiHit)["S-1"];
					_I = idDecoder(a_DigiHit)["I"];
					_J = idDecoder(a_DigiHit)["J"];
					_K = idDecoder(a_DigiHit)["K-1"];

					//int NCont=a_DigiHit->getNMCParticles();
					int NCont=a_DigiHit->getNMCContributions();
					float mainCont=0;
					for(int n=0;n<NCont;n++){
						if(a_DigiHit->getEnergyCont(n)>mainCont) {
							mainCont=a_DigiHit->getEnergyCont(n);
							_Time=a_DigiHit->getTimeCont(n);
							_MCPID=a_DigiHit->getPDGCont(n);
						}
					}

					_outputDigiHit->Fill();
				}
			}
			_outputTree->Fill();

			num++;
		}catch (lcio::DataNotAvailableException err) { }
	}
}

void AnaHit::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


