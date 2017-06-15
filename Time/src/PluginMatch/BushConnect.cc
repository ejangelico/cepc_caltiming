#include <BushConnect.hh>
#include <ArborTool.hh>
#include <ArborToolLCIO.hh>
#include <DetectorPos.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/TrackImpl.h>
#include "UTIL/CellIDDecoder.h"
#include "HelixClass.hh"		
#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <sstream>
#include <TH1.h>
#include <TVector3.h>
//#include <ArborTrack.hh>

using namespace std;
using namespace TMath;

//const string ECALCellIDDecoder  = "M:3,S-1:3,I:9,J:9,K-1:6";
const string ECALCellIDDecoder  = "S-1:8,M:8,I:16,K-1:12,GRZone:4,J:16";

BushConnect aBushConnect ;
BushConnect::BushConnect()
	: Processor("BushConnect"),
        _output(0)
{

	_FLAG_DIAGNOSIS = 0;
	registerProcessorParameter( "FlagDiagnosis" ,
                             "If actived, Store intermediate Cluster collections into LCIO" ,
                                _FLAG_DIAGNOSIS ,
                                _FLAG_DIAGNOSIS);

	_FLAG_MCPMIMIC = 0;
	registerProcessorParameter( "MCPMIMIC" ,
                             "If actived, MIMIC MCParticle into Tracks" ,
                                _FLAG_MCPMIMIC ,
                                _FLAG_MCPMIMIC);

/*
	_description = "Track Cluster Linking. Track info represented by MCTruth at this moment" ;
	registerProcessorParameter( "TimeCut" ,
		             "The name of the file to which the pion ROOT tree will be written" ,
				_TimeCut ,
				_TimeCut);

	registerProcessorParameter( "HitAbsCut" ,
		             "The name of the file to which the pion ROOT tree will be written" ,
				_HitAbsCut ,
				_HitAbsCut);
*/
}

void BushConnect::init() {
	printParameters();
	Cluflag.setBit(LCIO::CHBIT_LONG);
}

void BushConnect::Clean(){
	MCPTrack_Type.clear();
	Track_Energy.clear();
	Track_Type.clear();
	Track_Theta.clear();
	Track_EndPoint.clear();
	Track_P3.clear();
	TrackStartPoint.clear();

	SortedTracks.clear();
	ChCoreID.clear();

	chargedclustercore.clear();
	selfmergedcluster.clear();
	non_chargedclustercore.clear();		//all clusters besides charged cluster core

	ClusterType_1stID.clear();
	CluFD.clear();
	CluT0.clear();
	CluCoG.clear();
	Clu_Depth.clear();
}

void BushConnect::MCPTrackSort(LCEvent* evtPP)
{
	LCCollection * InputMCPColl = evtPP->getCollection("MCParticle");
	int NMCP = InputMCPColl->getNumberOfElements();

	SortedTracks.clear();
	LCCollection *mcpmimictrack = new LCCollectionVec(LCIO::TRACK);

	std::vector<int> CHMCPIndex; 
	std::vector<float> CHMCPEn; 
	CHMCPIndex.clear();
	CHMCPEn.clear();
	Track_EndPoint.clear();
	Track_Energy.clear();
	Track_Type.clear();

	TVector3 VTX, EndP;
	float MCPEn = 0;
	float SumEn = 0;  
	float SumChEn = 0; 

	for(int s = 0; s < NMCP; s++)
	{
		MCParticle * a_MCP = dynamic_cast<MCParticle*>(InputMCPColl->getElementAt(s));
		VTX = a_MCP->getVertex();
		EndP = a_MCP->getEndpoint();
		if(   VTX.Mag() < 10 && EndP.Mag() > 10  )
		{
			MCPEn = a_MCP->getEnergy();	
			SumEn += MCPEn; 
			if(a_MCP->getCharge() != 0)
			{
				CHMCPEn.push_back(MCPEn);
				CHMCPIndex.push_back(s);
				SumChEn += MCPEn; 
			}
		}
	}

	std::vector<int> SortCHMCPIndex = SortMeasure(CHMCPEn, 1);
	TVector3 MCP, CheckMCP; 

	for(int i = 0; i < (int) SortCHMCPIndex.size(); i++)
	{
		MCParticle * a_MCP = dynamic_cast<MCParticle*>(InputMCPColl->getElementAt(CHMCPIndex[SortCHMCPIndex[i]]));
		// cout<<i<<"th Sorted Charged MCP En "<<a_MCP->getEnergy()<<endl; 
		MCP = a_MCP->getMomentum();
		TrackImpl * a_trk = new TrackImpl();
		a_trk->setD0(0);
		a_trk->setZ0(0);
		a_trk->setOmega(0.00105/MCP.Perp());
		a_trk->setPhi(MCP.Phi());
		a_trk->setTanLambda( 1.0/tan(MCP.Theta()) );
		SortedTracks.push_back(a_trk);
		mcpmimictrack->addElement(a_trk);

		MCPTrack_Type[a_trk] = a_MCP->getPDG();
		Track_Energy[a_trk] = a_MCP->getEnergy();
		Track_EndPoint[a_trk] = ImpactPoint( MCP, a_MCP->getCharge(), 2245, 1845, BField);	// ECALHalfZ, ECALRadius, 
		Track_Type[a_trk] = 999;
		Track_P3[a_trk] = MCP;
	}
	evtPP->addCollection(mcpmimictrack, "MIMICTrack");
	cout<<SumChEn<<"/"<<SumEn<<" Charged Energy percentage "<<SumChEn/SumEn<<endl;
 
	/*
	if(abs(SumEn - 250) > 2)
	{
		cout<<"Warning, Energy Double Counted at MCTruth Level!!"<<endl; 
	}
	*/
}


void BushConnect::TrackSort(LCEvent* evtPP) //, &std::map<Track*, int>Track_Tpye, &std::map<Track*, float> Track_Energy)
{

	LCCollection * TrackColl = evtPP->getCollection("MarlinTrkTracks");

	int NTrk = TrackColl->getNumberOfElements();
	float D0 = 0;
	float Z0 = 0;
	int NTrkHit = 0;
	const float mass = 0.139;	//Pion Mass
	TVector3 EndPointPos, StartPointPos; 
	int TrackType = 0; 

	std::vector<Track*> tracks_HQ_Barrel; 
	std::vector<Track*> tracks_HQ_Endcap;
	std::vector<Track*> tracks_HQ_Shoulder;
	std::vector<Track*> tracks_HQ_Forward; 
	std::vector<Track*> tracks_MQ_Barrel;
	std::vector<Track*> tracks_MQ_Endcap;
	std::vector<Track*> tracks_MQ_Shoulder;
	std::vector<Track*> tracks_MQ_Forward;
	std::vector<Track*> tracks_Vtx; 
	std::vector<Track*> tracks_LQ; 
	std::vector<Track*> tracks_LE; 
	std::vector<Track*> curr_tracks;

	Track_EndPoint.clear();

	tracks_HQ_Barrel.clear();
	tracks_HQ_Endcap.clear();
	tracks_HQ_Shoulder.clear();
	tracks_HQ_Forward.clear();
	tracks_MQ_Barrel.clear();
	tracks_MQ_Endcap.clear();
	tracks_MQ_Shoulder.clear();
	tracks_MQ_Forward.clear();
	tracks_Vtx.clear();
	tracks_LQ.clear();
	tracks_LE.clear();

	std::vector<Track*> tracks_ILL;
	tracks_ILL.clear();
	std::vector<Track*> tracks_preInteraction;
	tracks_preInteraction.clear();	//Used to denote pion and electron interaction inside TPC/Tracker. Simply vetoed for avoid double counting... but muon may still be problematic. Better way of treating would be find the cascade photons & tracks - clusters, and veto all the daughters instead of mother. Similar can done for Kshort...
	// Condition, tracks_head to others tail. head position far from boundary. and, track energy >= sum of cascade

	std::vector<int> TrackOrder; 
	TrackOrder.clear();	
	std::map<Track*, int> Track_Index; 
	Track_Index.clear();
	Track_Energy.clear();
	Track_Type.clear();
	Track_P3.clear();
	Track_EndPoint.clear();
	TrackStartPoint.clear();

	for(int i0 = 0; i0 < NTrk; i0++)
	{
		Track* a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i0 ) );
		NTrkHit = a_Trk->getTrackerHits().size();		
		EndPointPos = (a_Trk->getTrackerHits()[NTrkHit - 1])->getPosition();	
		StartPointPos = (a_Trk->getTrackerHits()[0])->getPosition();
		Track_EndPoint[a_Trk] = EndPointPos;
		TrackStartPoint[a_Trk] = StartPointPos;

		HelixClass * TrkInit_Helix = new HelixClass();
		TrkInit_Helix->Initialize_Canonical(a_Trk->getPhi(), a_Trk -> getD0(), a_Trk -> getZ0(), a_Trk -> getOmega(), a_Trk->getTanLambda(), BField);
		float TrackEn = mass*mass;

		for (int q3 = 0; q3 < 3; q3 ++)
		{
			TrackEn += (TrkInit_Helix->getMomentum()[q3])*(TrkInit_Helix->getMomentum()[q3]);
		}
		TVector3 TrkMom(TrkInit_Helix->getMomentum()[0],TrkInit_Helix->getMomentum()[1],TrkInit_Helix->getMomentum()[2]);

		TrackEn = sqrt(TrackEn);
		Track_Energy[a_Trk] = TrackEn;
		Track_Theta[a_Trk] = TrkMom.Theta();
		// Track_Phi[a_Trk] = TrkMom.Phi();
		Track_P3[a_Trk] = TrkMom;		

		delete TrkInit_Helix;
	}

	TVector3 currEp, currSp;
	float currMotherEn = 0;
	float sumDauEn = 0; 

	for(int i1 = 0; i1 < NTrk; i1++)
	{
		Track* a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i1 ) );		
		currEp = Track_EndPoint[a_Trk];

		if( currEp.Perp() < 1600 && currEp.Perp() > 400 && abs(currEp.Z()) < 2000 )	//Only check 
		{
			currMotherEn = Track_Energy[a_Trk];
			sumDauEn = 0;	
			for(int i2 = 0; i2 < NTrk; i2++)
			{
				Track* b_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i2 ) );
				if(i2 != i1)
				{
					currSp = TrackStartPoint[b_Trk];
					if( (currEp - currSp).Mag() < 40  )
						sumDauEn += Track_Energy[b_Trk];
				}
			}
			if(currMotherEn + 0.1 > 0.9*sumDauEn && currMotherEn > 3 && sumDauEn > 0 )	//Some protection is always needed...
			{
				tracks_preInteraction.push_back(a_Trk);
			}
		}
	}

	for(int t0 = 0; t0 < NTrk; t0++)
	{
		Track* a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( t0 ) );
		D0 = a_Trk->getD0();
		Z0 = a_Trk->getZ0();
		NTrkHit = a_Trk->getTrackerHits().size();
		TrackerHit * last_hit = a_Trk->getTrackerHits()[NTrkHit - 1];
		EndPointPos = last_hit->getPosition();
		Track_EndPoint[a_Trk] = EndPointPos;
		StartPointPos = (a_Trk->getTrackerHits()[0])->getPosition();
		Track_Index[a_Trk] = t0;

		if( NTrkHit > 9 || (fabs(EndPointPos.Z()) > LStar - 500 && EndPointPos.Perp() < TPCInnerRadius ) || fabs(EndPointPos.Z()) > ECALHalfZ - 200  )		// Min requirement for track quality
		{	// LStar - 500, suppose to be the last Disk Position

			if( find(tracks_preInteraction.begin(), tracks_preInteraction.end(), a_Trk ) != tracks_preInteraction.end() )
			{
				cout<<"So We Drop it! "<<Track_Energy[a_Trk]<<endl; 
				continue; 
			}

			TrackType = 0;
			if((Track_Energy[a_Trk] < 1.0 && fabs(Track_Theta[a_Trk]-1.57)< 0.4) || (fabs(Track_Theta[a_Trk]-1.57) >= 0.4 && log10(Track_Energy[a_Trk]) < -(fabs(Track_Theta[a_Trk]-1.57)-0.4)*0.2/0.3 ))
			{
				TrackType = 100;
			}
			else if( fabs(EndPointPos.Z()) > ECALHalfZ - 500 && EndPointPos.Perp() > TPCOuterRadius - 300  )	//Shoulder
			{
				TrackType = 30;
			}
			else if( fabs(EndPointPos.Z()) > LStar - 500 && EndPointPos.Perp() < TPCInnerRadius )		//Forward
			{
				TrackType = 40;
			}
			else if( EndPointPos.Perp() > TPCOuterRadius - 100 )		//Barrel
			{
				TrackType = 10;
			}
			else if( fabs(EndPointPos.Z()) > ECALHalfZ - 200 )		//Endcap
			{
				TrackType = 20; 
			}

			if( fabs(D0) < 1 && fabs(Z0) < 1 )
			{
				TrackType += 1;
			}

			Track_Type[a_Trk] = TrackType; 

			if(TrackType == 11)
				tracks_HQ_Barrel.push_back(a_Trk);
			else if(TrackType == 21)
				tracks_HQ_Endcap.push_back(a_Trk);
			else if(TrackType == 31)
				tracks_HQ_Shoulder.push_back(a_Trk);
			else if(TrackType == 41)
				tracks_HQ_Forward.push_back(a_Trk);
			else if(TrackType == 10)
				tracks_MQ_Barrel.push_back(a_Trk);
			else if(TrackType == 20)
				tracks_MQ_Endcap.push_back(a_Trk);
			else if(TrackType == 30)
				tracks_MQ_Shoulder.push_back(a_Trk);
			else if(TrackType == 40)
				tracks_MQ_Forward.push_back(a_Trk);
			else if(TrackType == 1)
				tracks_Vtx.push_back(a_Trk);
			else if(TrackType == 101)
				tracks_LE.push_back(a_Trk);
			else if( (StartPointPos.Mag() > 50 && EndPointPos.Mag() < 1000 && NTrkHit < 50) || TrackType == 100  )
				tracks_ILL.push_back(a_Trk);
			else
				tracks_LQ.push_back(a_Trk);
		}
	}

	std::vector<float > currTrkMomentum;
	std::vector<int> currTrkIndex;

	for(int t1 = 0; t1 < 11; t1++)
	{
		currTrkMomentum.clear();
		currTrkIndex.clear();
		curr_tracks.clear();
		if(t1 == 0)
			curr_tracks = tracks_HQ_Endcap;
		else if(t1 == 1)
			curr_tracks = tracks_HQ_Barrel;
		else if(t1 == 2)
			curr_tracks = tracks_MQ_Endcap;
		else if(t1 == 3)
			curr_tracks = tracks_MQ_Barrel;
		else if(t1 == 4)
			curr_tracks = tracks_HQ_Shoulder;
		else if(t1 == 5)
			curr_tracks = tracks_MQ_Shoulder;
		else if(t1 == 6)
			curr_tracks = tracks_HQ_Forward;
		else if(t1 == 7)
			curr_tracks = tracks_MQ_Forward;
		else if(t1 == 8)
			curr_tracks = tracks_Vtx;
		else if(t1 == 9)			
			curr_tracks = tracks_LQ; 
		else if(t1 == 10)			
			curr_tracks = tracks_LE; 


		int N_currTrack = curr_tracks.size();

		for(int t2 = 0; t2 < N_currTrack; t2++)
		{
			Track* tmpTrk = curr_tracks[t2];
			currTrkMomentum.push_back(Track_Energy[tmpTrk]);
		}

		currTrkIndex = SortMeasure(currTrkMomentum, 1);

		for(int t3 = 0; t3 < N_currTrack; t3++)
		{
			Track* b_tmpTrk = curr_tracks[currTrkIndex[t3]];
			if(t1 < 9 || Track_Energy[b_tmpTrk] < 10)
				TrackOrder.push_back(Track_Index[b_tmpTrk]);
		}
	}

	for(unsigned int t4 = 0; t4 < TrackOrder.size(); t4++)
	{
		Track* b_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( TrackOrder[t4] ) );
		SortedTracks.push_back(b_Trk);
	}
}

void BushConnect::BushSelfMerge(LCEvent * evtPP)
{
	LCCollection * CaloClu = evtPP->getCollection("EHBushes");      //A sort here should be helpful
	int NClu = CaloClu->getNumberOfElements();

	std::vector<Cluster* > Core_1st; 
	std::vector<Cluster* > Frag_1st;
	std::vector<Cluster* > UnId_1st; 
	Core_1st.clear();
	Frag_1st.clear();
	UnId_1st.clear();

	float CluDepth = 0; 
	float CluEn = 0;
	int CluSize = 0; 
	TVector3 PosCluSeed, PosSeedDiff, PosCoGDiff, PosSeedA, PosSeedB; 
	TVector3 AxisCluA,AxisCluB;

	int NJoints = 0; 	
	int SmallCluSize = 0; float DeeperDepth = 0; float LaterT0 = 0; 
	float Depth_A = 0; float Depth_B = 0;
	int Size_A = 0; int Size_B = 0; 

	TMatrixF FlagMerge(NClu, NClu);

	for(int i0 = 0; i0 < NClu; i0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CaloClu->getElementAt(i0));
		CluFD[a_clu] = FDV3(a_clu, ECALCellIDDecoder);
		CluT0[a_clu] = ClusterT0(a_clu);
		CluCoG[a_clu] = ClusterCoG(a_clu);
		PosCluSeed = a_clu->getPosition();
		if(_FLAG_MCPMIMIC)
		{
			Clu_Depth[a_clu] =  DisSeedSurfaceSimple(PosCluSeed);
		}
		else
		{
			Clu_Depth[a_clu] = DisSeedSurface(PosCluSeed);
		}
	}
	// CluCoG_Top20Percent Might be used to improve Photon Split Remerge performance

	float cell = 1;
	for(int s0 = 0; s0 < NClu; s0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CaloClu->getElementAt(s0));
		PosSeedA = a_clu->getPosition();
		Depth_A = Clu_Depth[a_clu];
		Size_A = a_clu->getCalorimeterHits().size();
		AxisCluA = PosSeedA-CluCoG[a_clu];

		for(int s1 = s0 + 1; s1 < NClu; s1++)
		{
			Cluster *b_clu = dynamic_cast<Cluster*>(CaloClu->getElementAt(s1));
			NJoints = JointsBetweenBush(a_clu, b_clu, 1.6*cell); //<10mm
			//NJoints = JointsBetweenBush(a_clu, b_clu, 1.2*cell); //>=10mm
			PosSeedB = b_clu->getPosition();
			Depth_B = Clu_Depth[b_clu];
			Size_B = b_clu->getCalorimeterHits().size();
			AxisCluB = PosSeedB-CluCoG[b_clu];

			PosSeedDiff = PosSeedA - PosSeedB;
			DeeperDepth = std::max(Depth_A, Depth_B);
			LaterT0 = std::max(CluT0[a_clu], CluT0[b_clu]);
			PosCoGDiff = CluCoG[a_clu] - CluCoG[b_clu];
			SmallCluSize = std::min( Size_A, Size_B );

			if(NJoints && PosSeedDiff.Mag() < 40 + 0.05*DeeperDepth )	//And depth...
			{

				//if( ( ( NJoints > 4 || (NJoints > 1 && SmallCluSize < 10) ) && DeeperDepth > 30 ) || NJoints > 4 )
				//if( SmallCluSize < 10 || LaterT0 > 10 || DeeperDepth > 40 || (NJoints > 8 && PosCoGDiff.Mag()*PosCoGDiff.Angle(PosSeedB) < 15) )
				//if( SmallCluSize < 10 || LaterT0 > 10 || DeeperDepth > 30 || (NJoints > 8 && PosSeedDiff.Mag()*PosSeedA.Angle(PosSeedB) < 15) )
				//if( (NJoints > 8 && PosCoGDiff.Mag()*PosCoGDiff.Angle(PosSeedB) < 15) || (NJoints>8 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 10*pow(cell,2))) && (pow(PosSeedDiff.Mag(),2)-pow(PosSeedDiff.Z(),2)<2.1*cell))
				if( ( ( NJoints>4 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 9.1*pow(cell,2) ) ) && (pow(PosSeedDiff.Mag(),2)-pow(PosSeedDiff.Z(),2)<2.1*pow(cell,2) ) ) || (NJoints>8 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 4.1*pow(cell,2))))//1mm,2mm
				//if( ( ( NJoints>8 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 25.1*pow(cell,2) ) ) && (pow(PosSeedDiff.Mag(),2)-pow(PosSeedDiff.Z(),2)<9.1*pow(cell,2) ) ) || (NJoints>4 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 16.1*pow(cell,2))))//Higher Energy
				//if( (NJoints>8 && ( (Depth_A<Depth_B && Abs( PosSeedDiff.Mag()*AxisCluA.Angle(PosSeedDiff)) < 25.1*pow(cell,2) ) || (Depth_A>Depth_B && Abs(PosSeedDiff.Mag()*AxisCluB.Angle(PosSeedDiff)) < 25.1*pow(cell,2) ) ) && AxisCluA.Angle(AxisCluB) < 1.5) || (NJoints>1 && CluCoG[a_clu].Angle(CluCoG[b_clu])<1.2)  )//Angle
				//if( ( ( NJoints>12 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 4.1*pow(cell,2) ) ) && (pow(PosSeedDiff.Mag(),2)-pow(PosSeedDiff.Z(),2)<1.1*pow(cell,2) ) ) || (NJoints>8 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 2.1*pow(cell,2))))//4mm,5mm,8mm
				//if( ( ( NJoints>12 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 2.1*pow(cell,2) ) ) && (pow(PosSeedDiff.Mag(),2)-pow(PosSeedDiff.Z(),2)<2.1*pow(cell,2) ) ) || (NJoints>8 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 1.5*pow(cell,2))))//8mm,10mm
				//if( ( ( NJoints>8 && (pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 1.1*pow(cell,2) ) ) && (pow(PosSeedDiff.Mag(),2)-pow(PosSeedDiff.Z(),2)<1.1*pow(cell,2) ) && (DeeperDepth<30 || SmallCluSize>10) ) || pow(PosCoGDiff.Mag(),2)-pow(PosCoGDiff.Z(),2) < 0.9*pow(cell,2))//20mm
				{	
					FlagMerge[s0][s1] = 1.0;
					FlagMerge[s1][s0] = 1.0;
				}
			}

			//Head Tail Connection. Could be more sophsticate && should be very strict.
			if( PosSeedA.Angle(PosSeedB) < 0.1 && PosSeedDiff.Mag() < 1000 && PosSeedDiff.Mag()*PosSeedA.Angle(PosSeedB) < 60 + 0.02*DeeperDepth && ((CluFD[a_clu] < 0.2 && Size_A > 6) || (CluFD[b_clu] < 0.2 && Size_B > 6)) )
			{
				if( (PosSeedA.Mag() > PosSeedB.Mag() && PosSeedA.Angle(PosSeedB - PosSeedA) < 0.2) || (PosSeedB.Mag() > PosSeedA.Mag() && PosSeedA.Angle(PosSeedA - PosSeedB) < 0.2) )
				{
					FlagMerge[s0][s1] = 2.0;
					FlagMerge[s1][s0] = 2.0;
				}
			}
		}
	}

	std::vector<Cluster*> OriInputEHBushes = CollClusterVec(CaloClu);
	TMatrixF MergeSYM = MatrixSummarize(FlagMerge);
	LCCollection* CloseMergedCaloClu = ClusterVecMerge( OriInputEHBushes, MergeSYM);

	std::map<Cluster*,float> MinDisSeedToBush;
	MinDisSeedToBush.clear();
	for(int i0 = 0; i0 < CloseMergedCaloClu->getNumberOfElements(); i0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CloseMergedCaloClu->getElementAt(i0));
		PosCluSeed = a_clu->getPosition();
		float tmpmindis = 1e10;
		for(int i1 = 0; i1 < CloseMergedCaloClu->getNumberOfElements(); i1++)
		{
			Cluster * b_clu = dynamic_cast<Cluster*>(CloseMergedCaloClu->getElementAt(i1));
			if(i1 != i0)
			{
				if(DisPointToBush(PosCluSeed,b_clu) < tmpmindis) 
					tmpmindis = DisPointToBush(PosCluSeed,b_clu);  
			}
		}
		MinDisSeedToBush[a_clu] = tmpmindis;
	}

	for(int i0 = 0; i0 < CloseMergedCaloClu->getNumberOfElements(); i0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CloseMergedCaloClu->getElementAt(i0));
		PosCluSeed = a_clu->getPosition();
		TVector3 PosCoG = ClusterCoG(a_clu);
		if(_FLAG_MCPMIMIC == 0)
		{
			CluDepth = DisSeedSurface(PosCluSeed);
		}
		else
		{
			CluDepth = DisSeedSurfaceSimple(PosCluSeed);
		}
		CluEn = a_clu->getEnergy();
		CluSize = a_clu->getCalorimeterHits().size();

		//if( CluSize > 20 && ( (CluEn > 2.0 + 0.002*CluDepth && PosCoG.Z()-PosCluSeed.Z()>22 ) || CluEn > 5.0) )
		//if( CluSize > 10 && (  PosCoG.Z()-PosCluSeed.Z()>22 ) && CluEn > 5.0) 
		if( CluEn > 5.0) 
		//if( CluEn > 1.5 )
		//if( CluEn > 0.4 )
		{
			Core_1st.push_back(a_clu);
		}
		//else if( (CluSize > 10 || CluEn > 0.2) && MinDisSeedToBush[a_clu] > 20 && ClusterT0(a_clu) < 1.2 ) // && (PhotonTag(a_clu) == 1 || ClusterFlag1st(a_clu) == 11  ))
		
		/*else if( (CluSize > 15 || CluEn > 0.5) && MinDisSeedToBush[a_clu] > 2*cell && ClusterT0(a_clu) < 1.2 ) // && (PhotonTag(a_clu) == 1 || ClusterFlag1st(a_clu) == 11  ))
		{
			Core_1st.push_back(a_clu);
		}
*/

		else if( CluSize < 10  || ClusterT0(a_clu) >0  || CluEn<2.0 || CluDepth>70 || PosCoG.Z()-PosCluSeed.Z()<15 )
		{
			Frag_1st.push_back(a_clu);
		}
		else
		{
			UnId_1st.push_back(a_clu);
		}
	}

	std::vector<Cluster* > UndefFrag_1stAB = ClusterAbsorbtion(UnId_1st, Frag_1st, 300, 0.02);
	selfmergedcluster = ClusterAbsorbtion(Core_1st, UndefFrag_1stAB, 300, 0.02);
	//LCCollection *CluAB_1st=ClusterVecColl(selfmergedcluster);
	//evtPP->addCollection(CluAB_1st, "CluAB_1st");
	//IsoHits
	LCCollection * col_IsoHit = evtPP->getCollection("AllIsolatedHits");
	std::vector<CalorimeterHit*> IsoHits = CollHitVec(col_IsoHit, 0);
	std::vector<Cluster*> NBBAbs = ClusterHitAbsorbtion(selfmergedcluster, IsoHits, 300); //_HitAbsCut);	// Huge??
	LCCollection *CluAB_1st=ClusterVecColl(NBBAbs);
	evtPP->addCollection(CluAB_1st, "CluAB_1st");
	//evtPP->addCollection(CloseMergedCaloClu, "CluAB_1st");

}

void BushConnect::TagCore(LCEvent * evtPP) 
{
	LCCollection *arborrecoparticle_ch = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollection *arborchargedcorecluster = new LCCollectionVec(LCIO::CLUSTER);
	arborchargedcorecluster->setFlag(Cluflag.getFlag());

	int NTrk = SortedTracks.size();
	int NClu = selfmergedcluster.size();
	int currTrackType = 0;
	float currTrkEn = 0;
	float DisMatrix_Track_Clu_E[NTrk][NClu];
	float TimeMatrix_Track_Clu_E[NTrk][NClu];
	float CluDepth = 0;
	float CoreMergeDistanceDepthCorrector = 0;

	TVector3 TrkEndPoint(0, 0, 0);	
	TVector3 CluPos(0, 0, 0);
	std::map<Cluster*, int> BushTouchFlag; 
	std::map<Track*, Cluster*> FCMap_Track_CHCore;
	std::map<int, int> Closest_Trk_Clu_Map;
	std::vector<Cluster*> TightLinkedCluster;
	TightLinkedCluster.clear();
	Closest_Trk_Clu_Map.clear();
	BushTouchFlag.clear();
	FCMap_Track_CHCore.clear();

	Clu_Depth.clear();

	for(int s0 = 0; s0 < NTrk; s0++)
	{
		for(int s1 = 0; s1 < NClu; s1++)
		{
			DisMatrix_Track_Clu_E[s0][s1] = 1.0E10;
			TimeMatrix_Track_Clu_E[s0][s1] = 1.0E10; 
		}
	}

	for(int t0 = 0; t0 < NClu; t0++)
	{
		Cluster * a_clu = selfmergedcluster[t0];
		TVector3 PosCluSeed = a_clu->getPosition();
		if(_FLAG_MCPMIMIC)
		{
			Clu_Depth[a_clu] =  DisSeedSurfaceSimple(PosCluSeed);
		}
		else
		{
			Clu_Depth[a_clu] = DisSeedSurface(PosCluSeed);
		}
	}

	//~~~~~~~ find the closest cluster first...

	for(int g0 = 0; g0 < NTrk; g0++)
	{
		Track* a_trk = SortedTracks[g0];
		float ClosestDis = 1.0E9;  
		int ClosestCluIndex = -1; 
		int ClosestNC = 1E9;
		float TrackEn = Track_Energy[a_trk];
		TrkEndPoint = Track_EndPoint[a_trk];

		currTrackType = Track_Type[a_trk];
		TVector3 TrkP3 = Track_P3[a_trk];

		for(int g1 = 0; g1 < NClu; g1++)
		{
			Cluster *fccand_bush = selfmergedcluster[g1]; //ecalbushes[g1];
			float* Dis = SimpleDisTrackClu(a_trk, fccand_bush);
			float Time = SimpleBushTimeTrackClu(a_trk, fccand_bush);
			int NC = SimpleBushNC(a_trk, fccand_bush);
			if(Dis[2] > -0.1)
			{
				DisMatrix_Track_Clu_E[g0][g1] = Dis[2];
				TimeMatrix_Track_Clu_E[g0][g1] = Time;
				if( Dis[2] < ClosestDis ) // && ThetaDiff < 0.05 + 0.1/Track_Energy[a_trk] )
				{
					ClosestDis = Dis[2]; 
					ClosestCluIndex = g1;
					ClosestNC = NC;
				}
			}
		}

		//Diag for mimic
		cout<<" Z R "<<TrkEndPoint.Z()<<" MM "<<TrkEndPoint.Perp()<< " ClosestCluIndex "<<ClosestCluIndex<<" ClosestDis "<<ClosestDis <<endl; 
		//End Diag

		if( ClosestDis < 15 + 15./TrackEn && ClosestCluIndex > -0.1 && (ClosestNC < 3 || abs(TrkP3.Theta() - 1.57) < 0.01 ) ) 
		{
			Cluster * candiclu = selfmergedcluster[ClosestCluIndex];//  ecalbushes[ClosestCluIndex];
			CluPos = candiclu->getPosition();
			float TrackEndPDis = (TrkEndPoint - CluPos).Mag();
			if( (TrackEndPDis < 400 + 400./TrackEn || (TrackEn < 1 && candiclu->getEnergy() < 2 )) && ( fabs(TrackEn - candiclu->getEnergy() ) < 3.0*sqrt(TrackEn) + 1.0 ||  candiclu->getEnergy() < 3) && (TrkEndPoint.Z()*CluPos.Z() > 0 || abs(TrkEndPoint.Z() - CluPos.Z()) < 100 ) ) // && AngDiff < 0.4 )
			{
				Closest_Trk_Clu_Map[g0] = ClosestCluIndex;
				BushTouchFlag[candiclu] = 0;
			}
		}
	}	//~~~~~~~ end of finding closest cluster

	for(int i0 = 0; i0 < NTrk; i0++)  //Dropped Size can exist
	{
		Track* a_trk = SortedTracks[i0];
		currTrackType = Track_Type[a_trk];
		currTrkEn = Track_Energy[a_trk];

		TrkEndPoint = Track_EndPoint[a_trk];
		TightLinkedCluster.clear();
		float fccanden = 0; 
		if( Closest_Trk_Clu_Map.find(i0) != Closest_Trk_Clu_Map.end() )
		{
			Cluster * closeClu = selfmergedcluster[Closest_Trk_Clu_Map[i0]];//ecalbushes[Closest_Trk_Clu_Map[i0]];
			TightLinkedCluster.push_back(closeClu);
		}

		for(int j0 = 0; j0 < NClu; j0++)
		{
			Cluster *fccand_bush = selfmergedcluster[j0];	//ecalbushes[j0];			
			float Dis = DisMatrix_Track_Clu_E[i0][j0]; //SimpleDisTrackClu(a_trk, fccand_bush);
			float BushTime = TimeMatrix_Track_Clu_E[i0][j0];
			CluDepth = Clu_Depth[fccand_bush];
			CluPos = fccand_bush->getPosition();
			int currCluType = ClusterFlag1st(fccand_bush);

			CoreMergeDistanceDepthCorrector = 0;
			if(CluDepth > 20)
				CoreMergeDistanceDepthCorrector = 20;
			else if(CluDepth > 10)
				CoreMergeDistanceDepthCorrector = 10;

			float ProjectiveDis_TrackEndP_CluPos = (CluPos - TrkEndPoint).Mag()*(CluPos - TrkEndPoint).Angle(CluPos);
			if(log10(BushTime) < 3.5 && BushTime > 0 && currTrackType != 101 && Dis < 7 + CoreMergeDistanceDepthCorrector && Dis > -0.1 && BushTouchFlag.find(fccand_bush) == BushTouchFlag.end() && (currTrkEn > 3 || fccand_bush->getEnergy() < 5 || currCluType == 13 ) && ProjectiveDis_TrackEndP_CluPos < 100 ) // && (TrackEndPDis < 400 || currTrackType != 11) && (CluPos.Z()*TrkEndPoint.Z() > 0 || fabs((CluPos-TrkEndPoint).Z()) < 20 ) && (fccanden - currTrkEn <  something)  )
			{
				TightLinkedCluster.push_back(fccand_bush);
				BushTouchFlag[fccand_bush] = currTrackType;
				fccanden += fccand_bush->getEnergy();		//Energy Constrain Potentially be added
			}
		}
		//Diag 4 MIMIC
		cout<<"Tagged Clusters "<<TightLinkedCluster.size()<<endl; 
		//end Diag
		if( TightLinkedCluster.size() > 0 ) // && EcalCoreEnergy + HcalCoreEnergy < 2.0*currTrkEn )...
		{
			ClusterImpl * chcorecluster_eh =  NaiveMergeClu(TightLinkedCluster);
			FCMap_Track_CHCore[a_trk] = chcorecluster_eh;
			chargedclustercore.push_back(chcorecluster_eh);
		}
	}

	// So far, all the rest is used for Diagnosis in Tagcore

	for(int t0 = 0; t0 < NClu; t0++)
	{
		Cluster *a_clu = selfmergedcluster[t0];
		if( BushTouchFlag.find(a_clu) == BushTouchFlag.end() ) //Might be a neutral core
		{
			non_chargedclustercore.push_back(a_clu);
		}
	}

	int Track_Core_ID = -99;

	for(int j5 = 0; j5 < NTrk; j5++)
	{
		Track* a_trk = SortedTracks[j5];
		Track_Core_ID = -99;
		currTrackType = Track_Type[a_trk];
		ReconstructedParticleImpl * chargeparticle = new ReconstructedParticleImpl();
		chargeparticle->setEnergy( Track_Energy[a_trk] );
		chargeparticle->setCharge(a_trk -> getOmega()/fabs(a_trk -> getOmega()));
		TVector3 Ptrack = Track_P3[a_trk];
		double currTrkP[3] = { Ptrack.X(), Ptrack.Y(), Ptrack.Z() };
		chargeparticle->setMomentum(currTrkP); 
		chargeparticle->addTrack( a_trk );
		Cluster* a_clu = FCMap_Track_CHCore[a_trk];
		if( FCMap_Track_CHCore[a_trk] )		// No really need to pertect, as quality will be controled in Part.Reco
		{
			ClusterImpl * chargedarborcluster =  NaiveCluImpl(FCMap_Track_CHCore[a_trk]);
			chargeparticle->addCluster(chargedarborcluster);
			arborchargedcorecluster->addElement(chargedarborcluster);
			if( _FLAG_MCPMIMIC==0 )
			{
				LCCollection * col_TPCTrk =  evtPP->getCollection("ClupatraTracks");	//really stupid!...
				Track_Core_ID = ClusterFlag(a_clu, a_trk, col_TPCTrk);
			}
			else
			{
				Track_Core_ID = MCPTrack_Type[a_trk];
			}
		}
		chargeparticle->setType(Track_Core_ID);
		arborrecoparticle_ch->addElement(chargeparticle);
		ChCoreID[chargeparticle] = Track_Core_ID;
	}

	evtPP->addCollection( arborchargedcorecluster, "ClusterChargedCore" );
	evtPP->addCollection( arborrecoparticle_ch, "ArborChargedCore" );

	if(_FLAG_DIAGNOSIS)	// Pre Identifications of the Neutral Core Type...
	{
		LCCollection *arborrecoparticle_ne = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		LCCollection *arborneutralcorecluster = new LCCollectionVec(LCIO::CLUSTER);
		arborneutralcorecluster->setFlag(Cluflag.getFlag());

		std::vector<Cluster*> EMCore; 
		std::vector<Cluster*> Undef_EcalCluster; 
		std::vector<Cluster*> NHCore; 
		std::vector<Cluster*> Frag;
		EMCore.clear();
		NHCore.clear();
		Frag.clear();
		Undef_EcalCluster.clear();

		std::map<Cluster*, float> DisToCloseCC; 
		std::map<Cluster*, float> DisToCloseTrk;
		DisToCloseCC.clear();
		DisToCloseTrk.clear();
		float DisToClosestCC = 1E10; 
		float CluEn = 0; int CluSize = 0; 

		CluFD.clear();
		CluT0.clear();

		for(unsigned int i1 = 0; i1 < non_chargedclustercore.size(); i1++)
		{
			Cluster *a_clu = non_chargedclustercore[i1]; //ecalbushes[i1];
			CluPos = a_clu->getPosition();
			CluEn = a_clu->getEnergy();
			CluSize = a_clu->getCalorimeterHits().size();
			CluFD[a_clu] = FDV3(a_clu, ECALCellIDDecoder);
			CluT0[a_clu] = ClusterT0(a_clu);
			if(_FLAG_MCPMIMIC == 0)
			{
				CluDepth = DisSeedSurface(CluPos);
			}
			else
			{
				CluDepth = DisSeedSurfaceSimple(CluPos);
			}

			DisToClosestCC = 1E10;
			for(unsigned int j2 = 0; j2 < chargedclustercore.size(); j2++)
			{
				Cluster* a_CC = chargedclustercore[j2];
				float tmpDis = BushDis(a_CC, a_clu);
				if(tmpDis < DisToClosestCC)
				{
					DisToClosestCC = tmpDis; 
				}
			}
			DisToCloseCC[a_clu] = DisToClosestCC;

			float EnCorrEM = 0; // DisToClosestCC; 
			float SizeCorrEM = 0; 
			float EnCorrHad = 0; 
			float SizeCorrHad = 0; 

			if(DisToClosestCC < 100 &&  DisToClosestCC > 50)
			{
				EnCorrEM = 0.5;
				SizeCorrEM = 5; 
				EnCorrHad = 1.5;	// + 0.001 * Depth?
				SizeCorrHad = 10; 
			}
			else if( DisToClosestCC < 50 )
			{
				EnCorrEM = 1.0; 
				SizeCorrEM = 10;
				EnCorrHad = 3.0;	// And Depth??
				SizeCorrHad = 15;
			}

			if( CluEn > 0.15 + EnCorrEM && CluSize > 10 + SizeCorrEM && CluDepth < 40 && CluT0[a_clu] < 2 && (DisToClosestCC > 80 || CluEn < 1.0 || CluFD[a_clu] > 0.1 + 0.15*log10(CluSize)) )
			{
				EMCore.push_back(a_clu);
			}
			else if( CluEn > 2 + EnCorrHad && CluSize > 10 + SizeCorrHad && CluDepth > 30 && CluT0[a_clu] < 10 )	// Relations...
			{
				NHCore.push_back(a_clu);
			}
			else if( CluEn > 2 && CluDepth < 50 && CluT0[a_clu] < 2 + 0.005*CluDepth)
			{
				Undef_EcalCluster.push_back(a_clu);
			}
			else
			{
				Frag.push_back(a_clu);
			}
		}

		cout<<"Classificiation Done"<<endl; 

		LCCollection *Non_CC=ClusterVecColl(non_chargedclustercore);
		evtPP->addCollection(Non_CC, "Non_CC");
		LCCollection *emcore=ClusterVecColl(EMCore);
		LCCollection *Undef_EC = ClusterVecColl(Undef_EcalCluster);
		LCCollection *nhcore = ClusterVecColl(NHCore);
		LCCollection *frag = ClusterVecColl(Frag);
		evtPP->addCollection(emcore, "EMCore");
		evtPP->addCollection(Undef_EC, "EM_Undef");
		evtPP->addCollection(nhcore, "NHCore");
		evtPP->addCollection(frag, "Frag");

		float NAMom[3] = {0, 0, 0};
		TVector3 BushSeedPos;
		int NeutralPFOID = 0;
		float CluEnergy = 0; 

		for(unsigned int j6 = 0; j6 < EMCore.size() + Undef_EcalCluster.size() +  NHCore.size(); j6++ )
		{
			Cluster * a_clu(0);
			NeutralPFOID = 0;
			if(j6 < EMCore.size())
			{
				a_clu = EMCore[j6];
				NeutralPFOID = 22; 	
			}
			else if(j6 < EMCore.size() + Undef_EcalCluster.size())
			{
				a_clu = Undef_EcalCluster[j6 - EMCore.size()];
				NeutralPFOID = 33;	//KaonLong
			}
			else
			{
				a_clu = NHCore[j6 - EMCore.size() - Undef_EcalCluster.size()];
				NeutralPFOID = 44; 
			}
			CluEnergy = ClusterEE(a_clu);// a_clu->getEnergy();
			BushSeedPos = a_clu->getPosition();
			ReconstructedParticleImpl * neutralparticle = new ReconstructedParticleImpl();
			neutralparticle->setEnergy( CluEnergy );

			neutralparticle->setMass( 0.0 );
			neutralparticle->setCharge( 0.0 );
			neutralparticle->setType(NeutralPFOID);
			NAMom[0] = CluEnergy*1.0/BushSeedPos.Mag()*BushSeedPos.X();
			NAMom[1] = CluEnergy*1.0/BushSeedPos.Mag()*BushSeedPos.Y();
			NAMom[2] = CluEnergy*1.0/BushSeedPos.Mag()*BushSeedPos.Z();
			neutralparticle->setMomentum( NAMom );

			ClusterImpl * neutralarborcluster =  NaiveCluImpl(a_clu);
			neutralparticle->addCluster(neutralarborcluster);
			arborneutralcorecluster->addElement(neutralarborcluster);
			arborrecoparticle_ne->addElement(neutralparticle);		
		}

		evtPP->addCollection( arborneutralcorecluster, "ClusterNeutralCore");		// Used Only For Diagnosis...
		evtPP->addCollection( arborrecoparticle_ne, "ArborNeutralCore");
	}
}

void BushConnect::ParticleReco( LCEvent * evtPP )
{
	LCCollection *arborrecoparticle = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollection *mergedclu_ch = new LCCollectionVec(LCIO::CLUSTER);
	LCCollection *mergedclu_ne = new LCCollectionVec(LCIO::CLUSTER);
	mergedclu_ch->setFlag(Cluflag.getFlag());
	mergedclu_ne->setFlag(Cluflag.getFlag());

	LCCollection * ChargedCore = evtPP->getCollection("ArborChargedCore");
	LCCollection * col_TPCTrk;
	if(_FLAG_MCPMIMIC==0)  col_TPCTrk = evtPP->getCollection("ClupatraTracks");
	LCCollection * col_IsoHit = evtPP->getCollection("AllIsolatedHits");
	std::vector<CalorimeterHit*> IsoHits = CollHitVec(col_IsoHit, 0);

	int NChargedObj = ChargedCore->getNumberOfElements();
	int NNeutralCluster = non_chargedclustercore.size();
	double DisMatrix_Core_Neutral[NChargedObj][NNeutralCluster][2];		//Define different types of distances; 

	Track * a_chargedTrk(0); 
	Track * a_neighbourTrk(0);
	Cluster * a_chargedClu(0), *a_NeCandiClu(0); 
	float CluDepth = 0;
	std::map<Cluster*, double> CluDepthMap; 
	CluDepthMap.clear();
	int currChargeCoreType = 0;  
	TVector3 CluPos; 

	// Per Track usage...
	std::vector<Cluster*> loosecandicluster; 
	std::vector<Cluster*> tightcandicluster;		//Muon potential candi?
	std::vector<Cluster*> mergedcluster; 			//tmp for each charged P
	std::vector<Cluster*> chargedclustercore_merged; 	//overall

	chargedclustercore_merged.clear();

	std::vector<double> reftightdis; 
	std::vector<double> refloosedis; 

	std::map<Cluster*, int> NNCTouchFlag; 
	std::vector<Track*> SecondIterTracks;
	SecondIterTracks.clear();

	TVector3 currTrkEnd, neighbourTrkEnd, LeadP; 

	for(int i = 0; i < NChargedObj; i++)
	{
		ReconstructedParticle * a_recoP_ch = dynamic_cast<ReconstructedParticle*>(ChargedCore->getElementAt(i));

		loosecandicluster.clear();
		tightcandicluster.clear();
		mergedcluster.clear();
		reftightdis.clear();
		refloosedis.clear();
		a_chargedTrk = a_recoP_ch->getTracks()[0];
		currTrkEnd = Track_EndPoint[a_chargedTrk];
		currChargeCoreType = ChCoreID[a_recoP_ch];
		int currTrkType = Track_Type[a_chargedTrk];

		float CurrClusterEnergy = 0;
		float CurrTrackEnergy = Track_Energy[a_chargedTrk];
		if(a_recoP_ch->getClusters().size() != 0)
		{
			a_chargedClu = a_recoP_ch->getClusters()[0];
			CurrClusterEnergy = a_chargedClu->getEnergy();
			mergedcluster.push_back(a_chargedClu);		//Actually can use this chance to question if previous energy are balance...
		}

		float MinDisToNoClusterTrk = 1.0E10; 
		float MinDisToOtherTrack = 1.0E10;

		for( int is = 0; is < NChargedObj; is++ )
		{
			if(is != i)
			{
				ReconstructedParticle * b_recoP_ch = dynamic_cast<ReconstructedParticle*>(ChargedCore->getElementAt(is));
				a_neighbourTrk = b_recoP_ch->getTracks()[0];
				neighbourTrkEnd = Track_EndPoint[a_neighbourTrk];
				float currDD = (neighbourTrkEnd - currTrkEnd).Mag();
				if( currDD < MinDisToOtherTrack )
				{
					MinDisToOtherTrack = currDD;
				}
			}
		}

		for(int j = 0; j < NNeutralCluster; j++)
		{
			a_NeCandiClu = non_chargedclustercore[j];
			float NeCandEn = a_NeCandiClu->getEnergy(); 
			CluPos = a_NeCandiClu->getPosition();
			if(_FLAG_MCPMIMIC == 0)
			{
				CluDepth = DisSeedSurface(CluPos);
			}
			else
			{
				CluDepth = DisSeedSurfaceSimple(CluPos);
			}
			CluDepthMap[a_NeCandiClu] = CluDepth; 	

			if( ClusterType_1stID[a_NeCandiClu] == 1 )   continue; 

			for(int k = 0; k < 2; k++)
			{
				DisMatrix_Core_Neutral[i][j][k] = 1.0E9;
			}

			if(CurrClusterEnergy > 1E-6)	//put by hand...
			{
				DisMatrix_Core_Neutral[i][j][0] = BushDis(a_chargedClu, a_NeCandiClu);
			}
			float* Dis = SimpleDisTrackClu(a_chargedTrk, a_NeCandiClu);
			DisMatrix_Core_Neutral[i][j][1] = Dis[2];

			if( NNCTouchFlag.find(a_NeCandiClu) == NNCTouchFlag.end() && ( currChargeCoreType == 0 || DisMatrix_Core_Neutral[i][j][0] < 1000 ) && currTrkType != 101)
			{			
				if( currChargeCoreType == 130 )			//Matched Muon, should ignore
				{
					if( DisMatrix_Core_Neutral[i][j][1] < 0.2*CluDepth && CluDepth > 200  )	//&& FD?
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );	//dependence on Cluster Flag & Clu Depth. use some more fancy sort algorithm...
					}
				}
				else if( currChargeCoreType == 131 )
				{
					if( DisMatrix_Core_Neutral[i][j][1] < 0.3*CluDepth && CluDepth > 150 )
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );	
					}
					else if( DisMatrix_Core_Neutral[i][j][1] < 0.5*CluDepth && CluDepth > 100 )
					{
						loosecandicluster.push_back(a_NeCandiClu);
						refloosedis.push_back( DisMatrix_Core_Neutral[i][j][1] );
					}
				}	
				else if( currChargeCoreType == 110  )		// Electron
				{
					if( DisMatrix_Core_Neutral[i][j][0] < 0.15*CluDepth + 15 )
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}			
				else if( currChargeCoreType == 111 )		// look behind... might be pion...
				{
					if( DisMatrix_Core_Neutral[i][j][0] < 0.1*CluDepth + 15 && DisMatrix_Core_Neutral[i][j][1] < 0.1*CluDepth + 10 )	//Define Brems Photon region for correct
					{
						tightcandicluster.push_back(a_NeCandiClu);
						if(DisMatrix_Core_Neutral[i][j][0] < DisMatrix_Core_Neutral[i][j][1])	// not fully adequate.
						{
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
						}
						else
						{
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
					else if( DisMatrix_Core_Neutral[i][j][0] < 0.2*CluDepth + 15 || DisMatrix_Core_Neutral[i][j][1] < 0.2*CluDepth + 15  )
					{	
						loosecandicluster.push_back(a_NeCandiClu);

						if(DisMatrix_Core_Neutral[i][j][0] < DisMatrix_Core_Neutral[i][j][1])   // not fully adequate.
						{
							refloosedis.push_back( DisMatrix_Core_Neutral[i][j][0] );
						}
						else
						{
							refloosedis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
				}
				else if( currChargeCoreType == 211 )	//Main Cluster distance oriented
				{
					if(DisMatrix_Core_Neutral[i][j][0] < 0.2*CluDepth)
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}
				else if( currChargeCoreType == 212 )	//Non_Matched
				{
					if(DisMatrix_Core_Neutral[i][j][0] < 10 + 0.5*CluDepth )	//Energy Dependence...
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
					else if(DisMatrix_Core_Neutral[i][j][0] < 10 + 0.4*CluDepth || DisMatrix_Core_Neutral[i][j][1] < 20 + 0.5*CluDepth )
					{
						loosecandicluster.push_back(a_NeCandiClu);
						refloosedis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}
				else if( currChargeCoreType == 0 ) // && a_recoP_ch->getEnergy() < 3 ) // && !FlagLowPtPz )	//
				{
					if(CluDepth < 20)
					{
						if(DisMatrix_Core_Neutral[i][j][1] < MinDisToNoClusterTrk)	//Tag minimal distance cluster... and see if it can be potentially linked.
						{
							MinDisToNoClusterTrk = DisMatrix_Core_Neutral[i][j][1];
						}
						if( MinDisToNoClusterTrk < 300 && abs(a_recoP_ch->getEnergy() - NeCandEn) < 1.5*a_recoP_ch->getEnergy() )	//some hard cut
						{
							tightcandicluster.push_back(a_NeCandiClu);
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
				}
				else
				{
					cout<<"Over balanced/Un matched/defined case: "<<a_recoP_ch->getEnergy()<<" ??? "<<currChargeCoreType<<endl; 
				}
			}
		}

		float totaltightcandiEn = 0; 
		float totalloosecandiEn = 0; 
		for(unsigned int s = 0; s < tightcandicluster.size(); s++)
		{
			totaltightcandiEn += tightcandicluster[s]->getEnergy();
		}

		for(unsigned int s = 0; s < loosecandicluster.size(); s++)
		{
			totalloosecandiEn += loosecandicluster[s]->getEnergy();
		}

		if( currChargeCoreType == 130 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 ) //  && CurrClusterEnergy < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy) )	//Frags...
				{
					mergedcluster.push_back( a_clu );		
					CurrClusterEnergy += a_clu->getEnergy();
				}
				else if( ClusterType_1stID[a_clu] < 10 && (CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy) ))
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 131 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && (CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy)))  ) 
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}

			//Maybe Some ID over here?...	//layers & numbers...	//BS Case ID

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				Cluster* a_clu = loosecandicluster[i2];
				if( ClusterType_1stID[a_clu] >= 10 || CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy))       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 110 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 111 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				Cluster* a_clu = loosecandicluster[i2];
				if( ClusterType_1stID[a_clu] >= 10 || fabs(CurrClusterEnergy + a_clu->getEnergy() - CurrTrackEnergy) < fabs(CurrClusterEnergy - CurrTrackEnergy) )       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}	
		}
		else if( currChargeCoreType == 211 )	// Matched
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 212)
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				Cluster* a_clu = loosecandicluster[i2];
				if( ClusterType_1stID[a_clu] >= 10 || fabs(CurrClusterEnergy + a_clu->getEnergy() - CurrTrackEnergy) < fabs(CurrClusterEnergy - CurrTrackEnergy) )       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 0 && reftightdis.size() > 0)
		{
			float mindis = 1.0E10;
			int minindex = 0; 

			for(unsigned int i1 = 0; i1 < reftightdis.size(); i1 ++)
			{
				if(reftightdis[i1] < mindis)
				{
					mindis = reftightdis[i1];
					minindex = i1; 
				}
			}

			Cluster* a_clu = tightcandicluster[minindex];	// Only 1? ...

			mergedcluster.push_back( a_clu );
		}
		else
		{
			cout<<"No_match, currChargeCoreType "<<currChargeCoreType<<endl; 
		}

		float CHCluEnergy = 0;

		for(int is = 0; is < int(mergedcluster.size()); is++)
		{       
			Cluster* a_TBM_clu = mergedcluster[is]; 
			CHCluEnergy += EnUltraHotCorr(a_TBM_clu->getEnergy(), a_TBM_clu);
			// CHCluEnergy += ClusterEE(a_TBM_clu);
		}

		if( !( CHCluEnergy < 1 && CurrTrackEnergy > 5 ) || (MinDisToOtherTrack < 100) )	// Need to check if exist nearby larger charged cluster: maybe absorbed together and only left tiny MIP tail in the ECAL //* bool closeTonearByEnergeticChargedShower = 0  // MIP like; should also protect against energies
		{
			for(int i2 = 0; i2 < int(mergedcluster.size()); i2++)
			{
				Cluster* a_TBM_clu = mergedcluster[i2];
				NNCTouchFlag[a_TBM_clu]	= 2; 		// can make use of this intereting flag...
			}

			float charge = a_chargedTrk -> getOmega()/fabs(a_chargedTrk -> getOmega());

			ReconstructedParticleImpl * chargeparticle = new ReconstructedParticleImpl();
			chargeparticle->setCharge(charge);
			chargeparticle->addTrack( a_chargedTrk );
			TVector3 Ptrack = Track_P3[a_chargedTrk];
			double currTrkP[3] = { Ptrack.X(), Ptrack.Y(), Ptrack.Z() };

			int flagEnergyFlow = 0;
			int currChargeCoreType2 = -99;

			ClusterImpl * chclustermerged =  NaiveMergeClu(mergedcluster);
			mergedclu_ch->addElement(chclustermerged);
			chargeparticle->addCluster(chclustermerged);
			chargedclustercore_merged.push_back(chclustermerged);

			if(_FLAG_MCPMIMIC) 
				currChargeCoreType2 = MCPTrack_Type[a_chargedTrk];				
			else
				currChargeCoreType2 = ClusterFlag(chclustermerged, a_chargedTrk, col_TPCTrk);

			if( (_FLAG_MCPMIMIC == 0 && ( currChargeCoreType2 == 130 || currChargeCoreType2 == 131 ) ) || ( _FLAG_MCPMIMIC == 1 && fabs(currChargeCoreType2) == 13  ) )
			{
				chargeparticle->setType( int(-13*charge) );
			}
			else if( (_FLAG_MCPMIMIC == 0 && (currChargeCoreType2 == 110 || currChargeCoreType2 == 111)) || ( _FLAG_MCPMIMIC == 1 && fabs(currChargeCoreType2) == 11 ) )
			{
				chargeparticle->setType( int(-11*charge) );
				if(CHCluEnergy > CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy) + 1)
				{
					flagEnergyFlow = 1; 
				}
			}
			else
			{
				if( _FLAG_MCPMIMIC == 0)
				{
					chargeparticle->setType( int(211*charge) );
				}
				else
				{
					chargeparticle->setType( MCPTrack_Type[a_chargedTrk] );
				}
				if(CHCluEnergy > CurrTrackEnergy + 1.2*sqrt(CurrTrackEnergy) + 1)
				{
					flagEnergyFlow = 2;
				}
			}

			if( (_FLAG_MCPMIMIC==0 && ( (currChargeCoreType2 != 130 && currChargeCoreType2 != 131 &&  CurrTrackEnergy > 15 && CHCluEnergy > 0.5 && CurrTrackEnergy > CHCluEnergy + 2*sqrt(CHCluEnergy) )  || (( currChargeCoreType2 == 130 || currChargeCoreType2 == 131 ) && CurrTrackEnergy > 100 && CHCluEnergy < 3 && chclustermerged->getCalorimeterHits().size() < 20 ) ) )  || (_FLAG_MCPMIMIC == 1 && abs(currChargeCoreType2) != 13 && CurrTrackEnergy > 15 && CHCluEnergy > 0.5 && CurrTrackEnergy > CHCluEnergy + 2*sqrt(CHCluEnergy)  ) )  
			{
				chargeparticle->setEnergy( CHCluEnergy );
				double CorrMom[3] = {CHCluEnergy/CurrTrackEnergy*currTrkP[0], CHCluEnergy/CurrTrackEnergy*currTrkP[1], CHCluEnergy/CurrTrackEnergy*currTrkP[2] }; 
				chargeparticle->setMomentum( CorrMom );
			}
			else
			{
				chargeparticle->setEnergy( CurrTrackEnergy );
				chargeparticle->setMomentum( currTrkP );
			}

			if( flagEnergyFlow )
			{
				ReconstructedParticleImpl * a_Ef_Ne_particle = new ReconstructedParticleImpl();
				a_Ef_Ne_particle->setEnergy( CHCluEnergy - CurrTrackEnergy );
				TVector3 corePos = chclustermerged->getPosition();
				float WFactor = (CHCluEnergy - CurrTrackEnergy)/corePos.Mag(); 
				float PFNEMom[3] = {WFactor*float(corePos.X()), WFactor*float(corePos.Y()), WFactor*float(corePos.Z())};
				a_Ef_Ne_particle->setMomentum(PFNEMom);
				a_Ef_Ne_particle->setMass( 0.0 );
				a_Ef_Ne_particle->setCharge( 0.0 );
				a_Ef_Ne_particle->setType(501);
				arborrecoparticle->addElement(a_Ef_Ne_particle);

				cout<<"Energy Flow Neutral Tagged "<<CHCluEnergy - CurrTrackEnergy<<endl; 
			}

			arborrecoparticle->addElement(chargeparticle);
		}
		else	// push non valid tracks, etc to second iteration, as those for PreInteracting ones
		{
			SecondIterTracks.push_back(a_chargedTrk);
			cout<<"Second Iter Track Found"<<endl; 
		}	
	}
	evtPP->addCollection(mergedclu_ch, "ArborCharged");

	std::vector<Cluster*> Ab_or_veto_clu;
	std::vector<Cluster*> BBCore; 
	Ab_or_veto_clu.clear();
	BBCore.clear();

	for(int p6 = 0; p6 < NNeutralCluster; p6 ++)
	{
		Cluster * c_clu = non_chargedclustercore[p6];
		if( NNCTouchFlag.find(c_clu) == NNCTouchFlag.end() )
		{
			if( ClusterType_1stID[c_clu] < 10 || c_clu->getEnergy() > 0.02 + 0.001*CluDepthMap[c_clu] )	//Cores
			{
				BBCore.push_back(c_clu);
			}
		} 
	}

	float NAMom[3] = {0, 0, 0};

	//Final Re-absorption
	std::vector<Cluster*> NBBNeutral; 
	NBBNeutral.clear();

	for(int s = 0; s < int (BBCore.size()); s++)
	{
		Cluster * a_clu = BBCore[s];
		TVector3 PosClu = a_clu->getPosition();
		float Depth = 0; 
		if(_FLAG_MCPMIMIC == 0)
		{
			Depth = DisSeedSurface(PosClu);
		}
		else
		{
			Depth = DisSeedSurfaceSimple(PosClu);
		}
		float CoreEnCorr = ClusterEE(a_clu);

		if(ClusterFlag1st(a_clu) == 11)
		{
			TVector3 BushSeedPos = a_clu->getPosition();
			ReconstructedParticleImpl * neutralparticle = new ReconstructedParticleImpl();
			neutralparticle->setType(22);
			TVector3 PP = ClusterCoG(a_clu);
			NAMom[0] = CoreEnCorr*1.0/PP.Mag()*PP.X();
			NAMom[1] = CoreEnCorr*1.0/PP.Mag()*PP.Y();
			NAMom[2] = CoreEnCorr*1.0/PP.Mag()*PP.Z();
			neutralparticle->setEnergy( CoreEnCorr );
			neutralparticle->setMass( 0.0 );
			neutralparticle->setCharge( 0.0 );
			neutralparticle->setMomentum( NAMom );
			ClusterImpl * a_neclu = NaiveCluImpl(a_clu);
			a_neclu->setEnergy( CoreEnCorr );	//Reset...
			neutralparticle->addCluster(a_neclu);
			mergedclu_ne->addElement(a_neclu);
			arborrecoparticle->addElement(neutralparticle);
		}
		else	// Distance to Charged Core > sth;
		{
			float MinDisToChCore = 1.0E9;
			float currDis = 0; 
			int NChCore = mergedclu_ch->getNumberOfElements();
			float closestChCluEn = 0; 			
			for(int t = 0; t < NChCore; t++)
			{
				Cluster * a_chclu = dynamic_cast<Cluster*>(mergedclu_ch->getElementAt(t));
				currDis = BushDis(a_chclu, a_clu);
				if(currDis < MinDisToChCore)
				{
					MinDisToChCore = currDis;
					closestChCluEn = a_chclu->getEnergy();	// Or the Trk En??
				}
			}
			if( MinDisToChCore > 0.4*(15 + closestChCluEn + Depth*0.01) || a_clu->getEnergy() > 2.0 )	//Joint Depth??
			{
				NBBNeutral.push_back(a_clu);
			}
		}
	}

	// Add: Neural Core Remerge & Energy Scale Recalculate, IsoHit Abso
	std::vector<Cluster*> NBBAbs = ClusterHitAbsorbtion(NBBNeutral, IsoHits, 100); //_HitAbsCut);	// Huge??
	std::vector<float> BBAbsEn; 
	BBAbsEn.clear();

	for(unsigned s1 = 0; s1 < NBBAbs.size(); s1++)
	{
		BBAbsEn.push_back(NBBAbs[s1]->getEnergy());
	}

	std::vector<int> BBAbsIndex = SortMeasure(BBAbsEn, 1);

	std::vector<Cluster *> NeutronCore;
	std::vector<Cluster *> NeutronFlag;
	NeutronCore.clear();
	NeutronFlag.clear();	

	for(unsigned int s2 = 0; s2 < NBBAbs.size(); s2++)	//Sort it; the first one must be a neutral core?
	{
		Cluster * p_clu = NBBAbs[BBAbsIndex[s2]];
		float currCluEn = p_clu->getEnergy();
		std::vector<float> CluTime = ClusterTime(p_clu);
		if( (currCluEn > 1.0 || (currCluEn > 0.5 && s2 < 2) )&& CluTime[0] < 10)
		{
			NeutronCore.push_back(p_clu);
		}
		else
		{
			NeutronFlag.push_back(p_clu);
		}
	}

	std::vector<Cluster *> Neutrons = ClusterAbsorbtion(NeutronCore, NeutronFlag, 200, 0.01);

	for(unsigned int s3 = 0; s3 < Neutrons.size(); s3++)
	{
		Cluster * a_clu = Neutrons[s3];
		float CoreEnCorr = ClusterEE(a_clu);
		TVector3 PosClu = a_clu->getPosition();
		float Depth = 0;
		if(_FLAG_MCPMIMIC == 0)
		{
			Depth = DisSeedSurface(PosClu);
		}
		else
		{
			Depth = DisSeedSurfaceSimple(PosClu);
		}

		if(CoreEnCorr > 0.1 + 0.003*Depth || a_clu->getCalorimeterHits().size() > 20)	//CluTimes[0] < _TimeCut &&
		{
			if( ClusterFlag1st(a_clu) == 11 )	// Photon ID
				cout<<"WARNING... Photons after neutron merge merged"<<endl; 
			ReconstructedParticleImpl * neutralparticle = new ReconstructedParticleImpl();
			neutralparticle->setType(2112);
			TVector3 PP = ClusterCoG(a_clu);
			NAMom[0] = CoreEnCorr*1.0/PP.Mag()*PP.X();
			NAMom[1] = CoreEnCorr*1.0/PP.Mag()*PP.Y();
			NAMom[2] = CoreEnCorr*1.0/PP.Mag()*PP.Z();
			neutralparticle->setEnergy( CoreEnCorr );
			neutralparticle->setMass( 0.0 );
			neutralparticle->setCharge( 0.0 );
			neutralparticle->setMomentum( NAMom );
			ClusterImpl * a_neclu = NaiveCluImpl(a_clu);
			a_neclu->setEnergy( CoreEnCorr );       //Reset...
			neutralparticle->addCluster(a_neclu);
			mergedclu_ne->addElement(a_neclu);
			arborrecoparticle->addElement(neutralparticle);
		}
	}

	evtPP->addCollection( arborrecoparticle, "ArborPFOs");
	evtPP->addCollection( mergedclu_ne, "ArborNeutral" );
}

void BushConnect::processEvent( LCEvent * evtP )
{
	if (evtP)
	{
		_eventNr = evtP->getEventNumber();
		if(_eventNr%1 == 0)
			cout<<endl<<"   "<<_eventNr<<"  events Processed: "<<endl; //<<" bfield: "<<BField<<endl;

		BushConnect::Clean();	
		if(_FLAG_MCPMIMIC)
			BushConnect::MCPTrackSort( evtP );
		else
			BushConnect::TrackSort( evtP );
		BushConnect::BushSelfMerge( evtP ); 	
		BushConnect::TagCore( evtP );		
		BushConnect::ParticleReco( evtP );
	}
}

void BushConnect::end()
{
	std::cout<<"Bush Connection Finished, ArborObject Formed"<<std::endl;	
}

