#ifndef _AnaHit_hh_
#define _AnaHit_hh_

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

class AnaHit  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new AnaHit ; }

		AnaHit();

		~AnaHit() {};

		void init();

		void processEvent( LCEvent * evtP );
            
                
		void end();

	protected:
		std::string _treeFileName;
		std::ostream *_output;
		int _overwrite;
		TTree *_outputMCP;
		TTree *_outputTree;
		TTree *_outputDigiHit;
                int _eventNr, num;
		
		int _PDG,_ParentPDG;
                float _MCPEn, _MCPP[3],_MCPVertex[3],_MCPEndpoint[3];
		float _Pos[3];
		float _THEn;

                float _HitE, _HitPosX, _HitPosY, _HitPosZ, _Time;
                int _ID0, _ID1, _MCPID;
		int _M,_S,_I,_J,_K;

};

#endif




