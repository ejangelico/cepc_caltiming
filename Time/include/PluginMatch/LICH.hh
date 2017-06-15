#ifndef _LICH_hh_
#define _LICH_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/ParticleID.h>
#include <EVENT/Track.h>
#include <IMPL/LCFlagImpl.h>
#include <TNtuple.h>
#include <TObject.h>

#include <TTree.h>
#include <TFile.h>
#include <TH3.h>
#include <TVector3.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"



class TTree;

  class LICH  : public marlin::Processor
  {
  public:

    Processor*  newProcessor() { return new LICH ; }
 
    LICH();
    
    ~LICH() {};
    int SubDeFlag(TVector3 inputPos);
    int ActiveLayers(std::vector<CalorimeterHit*> clu, const std::string& encoder_str);
    int NHScaleV2( const std::string& encoder_str, std::vector<CalorimeterHit*> clu0, int RatioX, int RatioY, int RatioZ );
    float FDV2(std::vector<CalorimeterHit*> clu, const std::string& encoder_str);
    float DisSeedSurface( TVector3 SeedPos );
    TVector3 ClusterCoG(Cluster * inputCluster);


    void init();
    void gearPara();
    void NewClusterFlag(Cluster* a_tree, Track* a_trk);
    void SamplePro( LCEvent * evtP );
    void processEvent( LCEvent * evtP );   

    void end();

  protected:
    std::string _treeFileName;
    std::string _FileName;
    std::string _sampleEn;
    std::string _inputPFO;
    std::string _inputMCP;
    std::string _outputPFO;
    std::string weightDir;
    int _Training;

    LCFlagImpl Cluflag;
    std::ostream *_output;
    std::vector<float> inputLevel;
    std::vector<int> inputEn;
    std::vector<float> inputsubPos;
    std::vector<float> en_bin;
    int NEn;
    int NPos;

   std::vector<std::string> inputDetectorModules;
    std::vector<TMVA::Reader*> reader1;
    std::vector<std::vector<TMVA::Reader*> > reader_all;
    int eventNr, MCPDG, NPFO, EcalNHit, HcalNHit, CluNHit, NLEcal, NLHcal,AL_Ecal, AL_Hcal, NH_ECALF10, NH_ECALL20, FPFO, FMu, FE, NMCP;
    float _EcalNHit, _HcalNHit, _CluNHit, _NLEcal, _NLHcal,_AL_Ecal, _AL_Hcal, _NH_ECALF10, _NH_ECALL20, LCEn, TotalEn;
    float maxDisHtoL, minDisHtoL, avDisHtoL, avEnDisHtoL, EcalEn, HcalEn, EClu, graDepth, cluDepth, graAbsDepth, maxDepth, minDepth, MaxDisHel, MinDisHel, FD_all, FD_ECAL, FD_HCAL, crdis, EEClu_L10, EEClu_R, EEClu_r, EEClu_p, rms_Ecal, rms_Hcal, rms_Ecal2, rms_Hcal2, av_NHE, av_NHH, FD_ECALF10, FD_ECALL20, dEdx, cosTheta, Phi, LCcosTheta;
    float EE, E_10, E_R, E_r, minAngle,minEn;
    float mva_e, mva_mu, mva_pi;
	float rminBarrelEcal;
	float rmaxBarrelEcal;
	float zminBarrelEcal;
	float zmaxBarrelEcal;
	
	float rminEndCapEcal;
	float rmaxEndCapEcal;
	float zminEndCapEcal;
	float zmaxEndCapEcal;
	
	float rminBarrelHcal;
	float rmaxBarrelHcal;
	float zminBarrelHcal;
	float zmaxBarrelHcal;
	
	float rminEndCapHcal;
	float rmaxEndCapHcal;
	float zminEndCapHcal;
	float zmaxEndCapHcal;
	
	float BField;
	float TrkEn;
	int NMuon;
	int NElectron;
	float RecoEn;
	int PFOID,MCID,MCMotherID,muMotherID;
	float PFOcosTheta,closedPar,ELike,MuLike,PiLike,closedMCP,ImpactPara,LinkMCEn,trkD0,trkZ0,sigD0,sigZ0;
	int Num,IniMu,IniE;




	TTree *piTree;
	TTree *muTree;
	TTree *eTree;
	TTree *otherTree;
	TTree *evtTree;
	TTree *ArborPFO;
	TFile *tree_file;
	
	int ClusterID;
	int LCID;
	int NEnBin;
	std::vector<float> mvacut_pi;
	std::vector<float> mvacut_mu;
	std::vector<float> mvacut_e;
	float mvaVal_pi, mvaVal_mu, mvaVal_e;


  };

#endif


