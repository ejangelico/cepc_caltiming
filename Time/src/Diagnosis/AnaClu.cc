#include <AnaClu.hh>
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

AnaClu aAnaClu ;
AnaClu::AnaClu()
	: Processor("AnaClu"),
	_output(0)
{
	_description = "CLuster Analysis" ;

	_treeFileName="AnaClu.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);
}

void AnaClu::init() {

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

	_outputPFO = new TTree("Cluster", "Cluster" );
	_outputPFO->SetAutoSave(32*1024*1024);  // autosave every 32MB

	_outputArborHit = new TTree("ArborHit", "ArborHit" );
	_outputArborHit->SetAutoSave(32*1024*1024);  // autosave every 32MB

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
	_outputTree->Branch("Num", &num,"Num/I");
	_outputTree->Branch("NClu", &_nclu,"NClu/I");
	_outputTree->Branch("ThetaMC", &_thetaMC, "ThetaMC/F");
	_outputTree->Branch("THEn",&_THEn,"THEn/F");
	_outputTree->Branch("TCEn",&_TCEn,"TCEn/F");
	_outputTree->Branch("LCEn",&_LCEn,"LCEn/F");
	_outputTree->Branch("LCPos", _LCPos, "Pos[3]/F");

	_outputPFO->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputPFO->Branch("Num", &num,"Num/I");
	_outputPFO->Branch("Pos", _Pos, "Pos[3]/F");
	_outputPFO->Branch("Theta", &_Theta, "Theta/F");
	_outputPFO->Branch("Phi", &_Phi, "Phi/F");
	_outputPFO->Branch("EClu", &_EClu, "EClu/F");
	_outputPFO->Branch("CluDepth", &_CluDepth, "CluDepth/F");
	_outputPFO->Branch("CluFD", &_CluFD, "CluFD/F");
	_outputPFO->Branch("CluT0", &_CluT0, "CluT0/F");
	_outputPFO->Branch("CluCOG", &_CluCOG, "CluCOG[3]/F");
	_outputPFO->Branch("NHits", &_CluNHit,"NHits/I");

	_outputArborHit->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputArborHit->Branch("HitX",&_HitPosX,"HitX/F");
	_outputArborHit->Branch("HitY",&_HitPosY,"HitY/F");
	_outputArborHit->Branch("HitZ",&_HitPosZ,"HitZ/F");
	_outputArborHit->Branch("HitEn",&_HitE,"HitEn/F");
	_outputArborHit->Branch("Time",&_Time,"Time/F");
	_outputArborHit->Branch("ID0",&_ID0,"ID0/I");
	_outputArborHit->Branch("ID1",&_ID1,"ID1/I");
	_outputArborHit->Branch("M",&_M,"M/I");
	_outputArborHit->Branch("S",&_S,"S/I");
	_outputArborHit->Branch("I",&_I,"I/I");
	_outputArborHit->Branch("J",&_J,"J/I");
	_outputArborHit->Branch("K",&_K,"K/I");
	_outputArborHit->Branch("MCPID", &_MCPID, "MCPID/I");

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

void AnaClu::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{ 

			_eventNr=evtP->getEventNumber();
			cout << "Event:" << _eventNr << endl;
			LCCollection * MCP = evtP->getCollection("MCParticle");
			_nclu = 0;	
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
					_thetaMC = TMath::ACos(_MCPP[2]/sqrt(_MCPP[0]*_MCPP[0]+_MCPP[1]*_MCPP[1]+_MCPP[2]*_MCPP[2]));
				}
				else if(a_MCP->getParents().size() != 0){
					_ParentPDG = a_MCP->getParent(0)->getPDG();
				}
				_outputMCP->Fill();
			}


			//LCCollection * col_Clu = evtP->getCollection("EHBushes");
			//LCCollection * col_Clu = evtP->getCollection("ArborNeutral");
			LCCollection * col_Clu = evtP->getCollection("CluAB_1st");
			LCCollection * col_CaloH = evtP->getCollection("DigiSiHit");

			_THEn = 0;
			_TCEn = 0;
			_LCEn = 0;


			for(int c0 = 0; c0 < col_Clu->getNumberOfElements(); c0++)
			{
				Cluster *a_clu = dynamic_cast<EVENT::Cluster *>(col_Clu->getElementAt(c0));
				_nclu++;
				float aCluEn = a_clu->getEnergy();
				_TCEn += aCluEn;
				if(aCluEn > _LCEn) _LCEn = aCluEn;
				for(unsigned int h2= 0; h2 < a_clu->getCalorimeterHits().size(); h2++)
				{
					CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[h2];

					float HitEn = a_hit->getEnergy();
				}


				_EClu = aCluEn;
				TVector3 CluPos;
				CluPos = a_clu->getPosition();
				_Pos[0] = CluPos.X();
				_Pos[1] = CluPos.Y();
				_Pos[2] = CluPos.Z();
				_Theta = CluPos.Theta();
				_Phi = CluPos.Phi();
				_CluDepth = DisSeedSurfaceSimple(CluPos);
				_CluT0 = ClusterT0(a_clu);
				TVector3 Cog=ClusterCoG(a_clu);
				_CluCOG[0] = Cog.X();
				_CluCOG[1] = Cog.Y();
				_CluCOG[2] = Cog.Z();

				string ECALCellIDDecoder = col_CaloH->getParameters().getStringVal(LCIO::CellIDEncoding);
//cout<<ECALCellIDDecoder<<endl;
				_CluFD = FDV3(a_clu,ECALCellIDDecoder);
				CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);

				_EcalEn = 0;
				_HcalEn = 0;
				_EcalNHit = 0;
				_HcalNHit = 0;
				_CluNHit = 0;

				std::vector<CalorimeterHit*> Ecalhits;
				std::vector<CalorimeterHit*> Hcalhits;
				std::vector<CalorimeterHit*> allhits;

				allhits.clear();
				Ecalhits.clear();
				Hcalhits.clear();

				TVector3 HitPos;

				int currCluNHits = a_clu->getCalorimeterHits().size();
				if (currCluNHits == 0) continue;
				_CluNHit = currCluNHits;
				
				float totHitEn = 0;
				float totHitEnDis = 0;
				float HitEn;

				_outputPFO->Fill();

			}


			//------------------Arbor Leading Cluster Hit----------------------

				const std::string ECALCellIDDecoder = "S-1:8,M:8,:I:16,K-1:12,GRZone:4,J:16";
				CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);

			for(int c0 = 0; c0 < col_Clu->getNumberOfElements(); c0++)
			{
				Cluster *a_clu = dynamic_cast<EVENT::Cluster *>(col_Clu->getElementAt(c0));
				float aCluEn = a_clu->getEnergy();
				if(aCluEn!=_LCEn) continue;

				TVector3 LCluPos;
				LCluPos = a_clu->getPosition();
				_LCPos[0] = LCluPos.X();
				_LCPos[1] = LCluPos.Y();
				_LCPos[2] = LCluPos.Z();

				for(unsigned int h2= 0; h2 < a_clu->getCalorimeterHits().size(); h2++)
				{
					CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[h2];

					_HitPosX=a_hit->getPosition()[0];
					_HitPosY=a_hit->getPosition()[1];
					_HitPosZ=a_hit->getPosition()[2];
					_HitE=a_hit->getEnergy();
					//HitEnError = a_hit->getEnergyError();
					_ID0=a_hit->getCellID0();
					_ID1=a_hit->getCellID1();
					_Time=a_hit->getTime();
					_MCPID = -999;

					_M = idDecoder(a_hit)["M"];
					_S = idDecoder(a_hit)["S-1"];
					_I = idDecoder(a_hit)["I"];
					_J = idDecoder(a_hit)["J"];
					_K = idDecoder(a_hit)["K-1"];

					_outputArborHit->Fill();
				}

			}


			//----------------------------DigiHit------------------------------
			LCCollection * col = evtP->getCollection("DigiSiHit");
			int numElements = col->getNumberOfElements();

			if (col->getTypeName() == LCIO::CALORIMETERHIT)
			{
				//CellIDDecoder<CalorimeterHit> idDecoder0( col ) ;

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
			_outputTree->Fill();

			//--------------------------------------------------------------------

			num++;
		}catch (lcio::DataNotAvailableException err) { }
	}
}

void AnaClu::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


