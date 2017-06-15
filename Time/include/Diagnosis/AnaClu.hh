#ifndef _AnaClu_hh_
#define _AnaClu_hh_

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

class AnaClu  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new AnaClu ; }

		AnaClu();

		~AnaClu() {};

		void init();

		void processEvent( LCEvent * evtP );
            
                
		void end();

	protected:
		std::string _treeFileName;
		std::ostream *_output;
		int _overwrite;
		TTree *_outputMCP;
		TTree *_outputTree;
		TTree *_outputPFO;
		TTree *_outputArborHit;
		TTree *_outputDigiHit;
                int _eventNr, num, _nclu;
                float _thetaMC;
		
		int _PDG,_ParentPDG;
                float _MCPEn, _MCPP[3],_MCPVertex[3],_MCPEndpoint[3];
		int _EcalNHit, _HcalNHit, _CluNHit; 
                float _EcalEn, _HcalEn, _EClu,_cluDepth, _maxDepth, _minDepth;
                float _avEnDisHtoL;
		float _Pos[3];
		float _LCPos[3];
		float _Theta, _Phi;
		float _CluDepth,_CluFD,_CluT0;
		float _CluCOG[3];

		float _THEn, _TCEn, _LCEn;

                float _HitE, _HitPosX, _HitPosY, _HitPosZ, _Time;
                int _ID0, _ID1, _MCPID;
		int _M,_S,_I,_J,_K;

};

#endif




