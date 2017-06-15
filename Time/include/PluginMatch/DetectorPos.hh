#ifndef DETECTORPOS_H_
#define DETECTORPOS_H_

#include "TVector3.h"
#include "TMatrixF.h"
#include <iostream>
#include <vector>
#include <string>


//#include "gear/BField.h"
//#include "gear/CalorimeterParameters.h"
//#include <marlin/Global.h>
//#include <EVENT/Cluster.h>
//#include "lcio.h"

//const TVector3 thetpcCenter(0.,0.,0.);
//const TVector3 bField = marlin::Global::GEAR->getBField().at(thetpcCenter);
//const float BField=bField.z();
const float BField = 3.5;
//ILD_Default
const float DHCALBarrelRadius = 1845.0;         //Octo
const float DHCALEndCapInnerZ = 2245.0;
const float ECALEndCapInnerZ = 2245.0;
const float ECALHalfZ = 2245.0; // mm, Endcap Ente
const float ECALRadius = 1845.; // mm... minimal part for the octagnle.
const float TPCRadius = 1845.0 ;

const float TPCOuterRadius = 1845.0; 
const float TPCInnerRadius = 325.0; 
const float LStar = 1500.0; 

//ILD_v2
/*
const float ECALHalfZ = 1900.0; // mm, Endcap Ente
const float ECALRadius = 1400.0;
const float TPCRadius = 1365.0 ;
const float DHCALBarrelRadius = 2058.0; // ?    these are only called by DepthFlag, which is not really used in Arbor and Merging.. Anyway
const float ECALEndCapInnerZ = 2000.0;  // ?
const float DHCALEndCapInnerZ = 2200.0; // ?
*/

const float DHCALCalibrationConstant = 0.11;

int BarrelFlag( TVector3 inputPos );

int DepthFlag( TVector3 inputPos );

TVector3 CalVertex( TVector3 Pos1, TVector3 Dir1, TVector3 Pos2, TVector3 Dir2 );

int TPCPosition( TVector3 inputPos );		//Used to tag MCParticle position, if generated inside TPC & Dead outside

float DisSeedSurface( TVector3 SeedPos );	//for a given position, calculate the distance to Calo surface ( ECAL )

float DisSeedSurfaceSimple( TVector3 SeedPos );

//float DisSeedSurfaceClu( Cluster* a_clu );

float DisTPCBoundary( TVector3 Pos );

float DistanceChargedParticleToCluster(TVector3 CPRefPos, TVector3 CPRefMom, TVector3 CluPosition);

#endif //
