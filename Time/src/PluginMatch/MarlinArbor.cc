#include <MarlinArbor.hh>
#include <ArborTool.hh>
#include <ArborToolLCIO.hh>
#include <ArborHit.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/LCGenericObject.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ClusterImpl.h>
#include "UTIL/CellIDDecoder.h"
#include "UTIL/Operators.h"

#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Rtypes.h>
#include <sstream>
#include <set>
#include <TVector3.h>
#include <vector>
#include <algorithm>

#include "DetectorPos.hh"

using namespace std;
using namespace lcio ;
using namespace marlin ;

extern linkcoll InitLinks; 
// extern linkcoll IterLinks_1; 
extern linkcoll IterLinks; 
extern linkcoll links_debug; 
extern branchcoll Trees; 
extern std::vector<int> IsoHitsIndex;

MarlinArbor aMarlinArbor ;
MarlinArbor::MarlinArbor()
        : Processor("MarlinArbor")
          ,
          _output(0)
{
        _description = " Tree Algorithm in Marlin" ;
}

void MarlinArbor::init() {

	printParameters();
}

void MarlinArbor::HitsPreparation()
{
	cout<<"Start to prepare Hits"<<endl;
}

void MarlinArbor::LinkVisulization( LCEvent * evtPP, std::string Name, std::vector<CalorimeterHit*> Hits, linkcoll inputLinks ) 
{

	LCCollectionVec *colllink = new LCCollectionVec(LCIO::LCRELATION);
	LCFlagImpl linkflag;
	linkflag.setBit(LCIO::CHBIT_LONG);
	colllink->setFlag(linkflag.getFlag());

	int NLink = inputLinks.size();
	std::pair<int, int> currlink; 
	for(int i0 = 0; i0 < NLink; i0++)
	{
		currlink = inputLinks[i0];

		CalorimeterHit* a_hit = Hits[currlink.first];
		CalorimeterHit* b_hit = Hits[currlink.second];

		LCRelationImpl *a_link = new LCRelationImpl(a_hit, b_hit);
		colllink->addElement( a_link );
	}

	evtPP->addCollection(colllink, Name);
}

void MarlinArbor::MakeIsoHits( LCEvent * evtPP, std::vector<CalorimeterHit*> inputCaloHits, std::string outputBushCollection )
{
	LCCollection * isohitcoll = new LCCollectionVec(LCIO::CALORIMETERHIT);

	string initString = "M:3,S-1:3,I:9,J:9,K-1:6";          //Need to verify
	isohitcoll->parameters().setValue(LCIO::CellIDEncoding, initString);

	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);                  
	flag.setBit(LCIO::CHBIT_ID1);
	flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);
	isohitcoll->setFlag(flag.getFlag());

	int nhit = inputCaloHits.size();

	for(int i = 0; i < nhit; i++)
	{
		CalorimeterHit* a_hit = inputCaloHits[i];
		CalorimeterHitImpl * collhit = new CalorimeterHitImpl();
		collhit->setPosition(a_hit->getPosition());
		collhit->setCellID0(a_hit->getCellID0());
		collhit->setCellID1(a_hit->getCellID1());
		collhit->setEnergy(a_hit->getEnergy());
		collhit->setTime(a_hit->getTime());
		isohitcoll->addElement(collhit);
	}

	evtPP->addCollection(isohitcoll, outputBushCollection);
}

void MarlinArbor::processEvent( LCEvent * evtP )
{
	if(evtP)
	{
		int EvtNr = evtP->getEventNumber();
		_eventNr = EvtNr; 

		std::cout<<EvtNr<<" events processed"<<std::endl;

		TVector3 currHitPos;
		std::vector< ArborHit > inputABHit; 
		std::vector< CalorimeterHit* > inputHits;  
		std::vector< std::vector<int> > Sequence; 
		int IsoHitSize = 0; 
		int LayerNum = 0; 
		float Depth = 0; 
		TVector3 hitPos;

		std::vector<CalorimeterHit*> IsoHits;

		try{
			LCCollection * CaloHitColl = evtP ->getCollection("DigiSiHit");
			int NHitsCurrCol = CaloHitColl->getNumberOfElements();
			CellIDDecoder<CalorimeterHit> idDecoder(CaloHitColl);
			for(int i2 = 0; i2 < NHitsCurrCol; i2++)
			{
				CalorimeterHit * a_hit = dynamic_cast<CalorimeterHit*>(CaloHitColl->getElementAt(i2));
				hitPos = a_hit->getPosition();
				if(a_hit->getEnergy() > 0.0 && a_hit->getTime() - hitPos.Mag()/300 < 25)
				{	
					currHitPos = a_hit->getPosition();
					Depth = DisSeedSurface(currHitPos);
					LayerNum = idDecoder(a_hit)["K-1"];
					ArborHit a_abhit(currHitPos, LayerNum, 0, Depth, 0, 1);
					inputABHit.push_back(a_abhit);
					inputHits.push_back(a_hit);
				}
			}
		}catch(lcio::DataNotAvailableException zero) { }

		Sequence = Arbor(inputABHit, 1, 7.15);   
		ClusterBuilding( evtP, "EHBushes", inputHits, Trees, 0 );
		LinkVisulization(evtP, "Links_init", inputHits, InitLinks);
		//		LinkVisulization(evtP, "Links_iter_1", inputHits, IterLinks_1);
		LinkVisulization(evtP, "Links_iter", inputHits, IterLinks);
		LinkVisulization(evtP, "Links_init_Debug", inputHits, links_debug);

		IsoHitSize = IsoHitsIndex.size();

		for(int i2 = 0; i2 < IsoHitSize; i2++)
		{
			CalorimeterHit* a_Isohit = inputHits[ IsoHitsIndex[i2] ];
			if(a_Isohit->getEnergy() > 0)	//Veto Trk End Hits
			{
				IsoHits.push_back(a_Isohit);
			}
		}
		MakeIsoHits(evtP, IsoHits, "AllIsolatedHits");
	}
}


void MarlinArbor::end()
{
	std::cout<<"Arbor Ends. Good luck"<<std::endl;
}
