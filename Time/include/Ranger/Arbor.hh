#ifndef _Arbor_hh_
#define _Arbor_hh_

#include "ArborHit.hh"
#include <string>
#include <iostream>
#include <TVector3.h>
#include "ArborTool.hh"

//namespace arbor {

  void init(float CellSize, float LayerThickness);
  
  void HitsCleaning(std::vector<ArborHit> inputHits );
  
  void HitsClassification( linkcoll inputLinks );
  
  void BuildInitLink(float Threshold);
  
  void LinkIteration(float Threshold);
  
  void BranchBuilding(float NeighbourThreshold); 
  
  void BushMerging();
  
  void BushAbsorbing();
  
  branchcoll Arbor( std::vector<ArborHit>, const float CellSize, const float LayerThickness); 
//}

#endif


