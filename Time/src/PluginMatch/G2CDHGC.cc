#include <G2CDHGC.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include "UTIL/CellIDDecoder.h"

#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TMath.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;

G2CDHGC aG2CDHGC ;
G2CDHGC::G2CDHGC()
	: Processor("G2CDHGC"),
	_output(0)
{
	_description = "Naive HGC Digitizer";

	registerProcessorParameter("CalibCalo" ,
                             "Global Calibration Constant" ,
                             _CalibCalo,
                             float(1.0));
}

void G2CDHGC::init() {
	printParameters();
	flag.setBit(LCIO::CHBIT_LONG);                  //To set position & ID1
	flag.setBit(LCIO::CHBIT_ID1);
	flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);
}

void G2CDHGC::processEvent( LCEvent * evtP )
{
	if (evtP)
	{
		try{
			LCCollection *inputHGCcol = evtP->getCollection( "SiCalCollection" ) ;

			int NHGCSimHit = inputHGCcol->getNumberOfElements();
			float HitEn = 0; 
			float currTime = 0;
			float HitStepEn = 0;
			float EmaxStep = 0;

			LCCollectionVec *hgccol = new LCCollectionVec(LCIO::CALORIMETERHIT);
			string idString = inputHGCcol->getParameters().getStringVal(LCIO::CellIDEncoding);
			hgccol->parameters().setValue(LCIO::CellIDEncoding, idString);
			hgccol->setFlag(flag.getFlag());

			LCCollectionVec *relcol = new LCCollectionVec(LCIO::LCRELATION);
			relcol->setFlag(flag.getFlag());

			CellIDDecoder<SimCalorimeterHit> idDecoder(inputHGCcol);
			int LayerNum = 0; 

			for(int k1 = 0; k1 < NHGCSimHit; k1++)
			{       
				SimCalorimeterHit * SimEcalhit = dynamic_cast<SimCalorimeterHit*>( inputHGCcol->getElementAt( k1 ) ) ;
				HitEn = SimEcalhit->getEnergy();
				currTime = 0;
				EmaxStep = 0; 
				for(int k=0; k<SimEcalhit->getNMCContributions(); k++)
                                {
					HitStepEn = SimEcalhit->getEnergyCont(k);
					if(HitStepEn > EmaxStep)
                                        {
						EmaxStep = HitStepEn;
						currTime = SimEcalhit->getTimeCont(k);
					}
				}
				if( HitEn > 1.0E-9 )
				//if( HitEn > 1.47E-5*5 )
				{
					CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
					calhit->setCellID0( SimEcalhit->getCellID0() );        //Assume 100% efficiency
					calhit->setCellID1( SimEcalhit->getCellID1() );
					calhit->setPosition( SimEcalhit->getPosition() );
					LayerNum = idDecoder( SimEcalhit )["K-1"];
					if(LayerNum)
						calhit->setEnergy( SimEcalhit->getEnergy()*_CalibCalo );          //Charge
					else
						calhit->setEnergy( SimEcalhit->getEnergy());
					calhit->setTime(currTime);
					calhit->setEnergyError(0);
					hgccol->addElement(calhit);
					LCRelationImpl *rel = new LCRelationImpl(calhit, SimEcalhit, 1.0);    //only keep the leading contribution
					relcol->addElement(rel);
				}
			}
			evtP->addCollection(hgccol, "DigiSiHit");
			evtP->addCollection(relcol, "DigiToSim");
		}
		catch (lcio::DataNotAvailableException zero) { }
	}
}

void G2CDHGC::end()
{
	std::cout<<"General Gas Digitizer FINISHED"<<std::endl;
}
