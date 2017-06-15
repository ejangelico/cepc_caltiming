#ifndef _AnaOverlay_hh_
#define _AnaOverlay_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>

class TTree;

class AnaOverlay  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new AnaOverlay ; }

		AnaOverlay();

		~AnaOverlay() {};

		void init();

		void processEvent( LCEvent * evtP );
            
                
		void end();

	protected:
		std::string _treeFileName;
		std::ostream *_output;
		int _overwrite;
		int _Merge;
		TTree *_outputMCP;
		TTree *_outputTree;
		TTree *_outputPFO;
		TTree *_outputArborHit;
		TTree *_outputDigiHit;
                int _eventNr, num, _nclu;
		
		int  _CluNHit; 
                float _EClu,_cluDepth, _maxDepth, _minDepth;
		float _Pos[3];

		float _THEn, _TCEn, _LCEn,_SCEn;
		float  _LCEnRSCEn,_LCEnPSCEn;

                float _HitE, _HitPosX, _HitPosY, _HitPosZ, _Time;
                int _ID0, _ID1, _MCPID;
		int _M,_S,_I,_J,_K;

};

#endif




