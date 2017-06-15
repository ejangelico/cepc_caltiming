#ifndef _TimeAna_hh_
#define _TimeAna_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class TimeAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new TimeAna ; }

		TimeAna();

		~TimeAna() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::ostream *_output;

		TTree *_outputTree;
		int _overwrite;
		unsigned int _eventNr;


		float _NeuAveTime;
		float _NeuPeakTime;
		float _NeuAngle;
		float _NeuNCount;
		float _FragAveTime;
		float _FragPeakTime;
		float _FragAngle;
		float _FragNCount;
		float _NeuEnergy;
		float _FragEnergy;
		float _NeuPos;
		float _NeuNHits;
		int _NClu;

				//float _MCPEn, _TrkEn, _CluEn, _Energy, _Charge;
		//float _MCP[3], _TrkVtxP[3], _TrkEndP[3], _P[3], _EnSeg[10];
		//float _TotalEn, _TotalChEn, _TotalNeuEn; 
		//int _NArborC, _NArborN; 
		//float _Depth, _LeadCDepth, _LeadCEn, _FD, _SumDE, _LeadDE; 
		//int _LeadCHits; 
};

#endif


