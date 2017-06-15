#ifndef ARBORTOOLLCIO_H_
#define ARBORTOOLLCIO_H_

#include "EVENT/LCCollection.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>
#include "EVENT/Cluster.h"
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include "HelixClass.hh" 
#include "TVector3.h"
#include "TMatrixF.h"
#include "lcio.h"
#include <vector>
#include <string>
#include "math.h"

using namespace lcio;

int NHScale( LCCollection *inputHit, Cluster * clu0, int RatioX, int RatioY, int RatioZ );


float ClusterT0(Cluster * a_Clu);

float FD( Cluster * clu, LCCollection *HitCollection );

int NHScaleV2( const std::string& encoder_str, std::vector<CalorimeterHit*> clu0, int RatioX, int RatioY, int RatioZ );

float FDV2( std::vector<CalorimeterHit*> clu, const std::string& encoder_str );

int NHScaleV3( const std::string& encoder_str, Cluster * clu0, int RatioX, int RatioY, int RatioZ );

float FDV3( Cluster * clu, const std::string& encoder_str );

float FD_I( std::vector<CalorimeterHit*> clu, const std::string& encoder_str, int InitSize );

void CheckConsistence(LCEvent *evtQQ, std::vector<std::string> SubClass);

int ActiveLayers(  std::vector<CalorimeterHit*> clu, const std::string& encoder_str );

float BushDis( Cluster * clu1, Cluster * clu2);

float* DisTrackClu(Track * a_trk, Cluster * a_clu);

float* SimpleDisTrackClu(Track * a_trk, Cluster * a_clu);

float SimpleBushTimeTrackClu(Track * a_trk, Cluster * a_clu);

int SimpleBushNC(Track * a_trk, Cluster * a_clu);

float DisPointToBush( TVector3 Pos1, Cluster * clu1);

TVector3 ImpactPoint( TVector3 MCP_Mom, int MCP_Charge, float Track_Half_Z, float Track_Radius, float B_Field ); //B_Field always along Z+ direction

TVector3 ClusterCoG(Cluster * inputCluser);

// HelixCluster_TightLink(float TrackEnergy, float currTrkTheta, float Dis);
// bool MIPFragFlag(Cluster * a_clu, float Dis, float TrkEn);
// std::pair<TVector3, TVector3> ClosestPointPair(Cluster *a_clu, Cluster *b_clu);
// TVector3 TrackOuterHit( Track* inputTrack ); // Find the outer most hit of a track
// LCCollection* ClusterMerge( LCCollection * inputClusterColl, TMatrixF ConnMatrix );

LCCollection* ClusterVecMerge( std::vector<Cluster*> inputClusters, TMatrixF ConnectorMatrix  );

LCCollection* ClusterVecColl( std::vector<Cluster*> inputClusters );

std::vector<Cluster*> CollClusterVec(LCCollection * input_coll );

std::vector<CalorimeterHit*> CollHitVec(LCCollection * input_coll, float EnergyThreshold);

std::vector<Cluster*> ClusterHitAbsorbtion( std::vector<Cluster*> MainClusters, std::vector<CalorimeterHit*> IsoHits, float DisThreshold );

std::vector<Cluster*> ClusterAbsorbtion( std::vector<Cluster*> MainClusters, std::vector<Cluster*> FragClusters, float thresholds, float DepthSlope );

ClusterImpl* NaiveMergeClu(std::vector<Cluster*> inputCluVec);

float DisHelixCluster(Cluster* a_clu, HelixClass* a_helix, float &Time);

int JointsBetweenBush(Cluster* a_Clu, Cluster* b_Clu, float CellSize);

ClusterImpl* NaiveCluImpl(Cluster* a0_clu);

TVector3 CluEndP(Cluster* a_clu);

std::vector<TVector3> BranchRefDir( Cluster* inputBranch );

std::vector<TVector3> CluRefDir(Cluster* a_clu);

void ClusterCollBuild( LCEvent * evtPP, std::string Name, std::vector<Cluster*> inputClusters );

void ClusterBuilding( LCEvent * evtPP, std::string Name,std::vector<CalorimeterHit*> Hits, std::vector< std::vector<int> > BranchOrder, int DHCALFlag );

int ClusterFlag(Cluster* a_tree, Track* a_trk, LCCollection *col_TPCTrk);

int ClusterFlag1st(Cluster* a_tree);

int PhotonTag(Cluster * a_clu);

float EnUltraHotCorr(float En, Cluster *a_tree);

float ClusterEE(Cluster* inputCluster);

float EMClusterEE( Cluster *inputCluster );

float HADClusterEE(float En, Cluster *inputCluster);

std::vector<float> ClusterTime(Cluster *inputCluster);
#endif //
