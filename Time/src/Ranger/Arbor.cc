#include <Arbor.hh>
#include <ArborHit.hh>
#include <TTree.h>
#include <algorithm>
#include <TMath.h>
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>

using namespace std;

typedef std::unordered_multimap<unsigned,unsigned> HitLinkMap;

std::vector<ArborHit> cleanedHits;

std::vector<int> LeafHitsIndex; 
std::vector<int> JointHitsIndex; 
std::vector<int> StarJointHitsIndex; 
std::vector<int> IsoHitsIndex;
std::vector<int> SimpleSeedHitsIndex;
std::vector<int> StarSeedHitsIndex;

HitLinkMap BackLinksMap;
linkcoll Links;
linkcoll InitLinks;
HitLinkMap alliterBackLinksMap;
linkcoll alliterlinks;
linkcoll links_debug; 
linkcoll IterLinks; 
HitLinkMap IterBackLinks;
branchcoll LengthSortBranchCollection;
branchcoll Trees; 

int NHits = 0; 
float InitLinkThreshold = 0; 
float IterLinkThreshold = 0;

void init( float CellSize, float LayerThickness ) 
{
	cleanedHits.clear();
	LeafHitsIndex.clear();
	JointHitsIndex.clear();
	StarJointHitsIndex.clear();
	IsoHitsIndex.clear();
	SimpleSeedHitsIndex.clear();
	StarSeedHitsIndex.clear();

	BackLinksMap.clear();
	Links.clear();
	InitLinks.clear();
	IterLinks.clear();
	IterBackLinks.clear();
	LengthSortBranchCollection.clear();
	Trees.clear();
	alliterlinks.clear();
	alliterBackLinksMap.clear();
	links_debug.clear();

	InitLinkThreshold = 2*LayerThickness - 0.01; 
	IterLinkThreshold = InitLinkThreshold ;	//2.5 a bit too large...
}

void HitsCleaning(std::vector<ArborHit> inputHits )
{
	cleanedHits = std::move(inputHits);        //Cannot Really do much things here. Mapping before calling
	NHits = cleanedHits.size();
}

void HitsClassification( linkcoll inputLinks )
{
	int NLinks =  inputLinks.size();

	LeafHitsIndex.clear();
	JointHitsIndex.clear();
	StarJointHitsIndex.clear();
	IsoHitsIndex.clear();
	SimpleSeedHitsIndex.clear();
	StarSeedHitsIndex.clear();

	std::pair<int, int> a_link;

	int BeginIndex[NHits];
	int EndIndex[NHits];

	for(int i0 = 0; i0 < NHits; i0++)
	{
		BeginIndex[i0] = 0;
		EndIndex[i0] = 0;
	}

	for(int j0 = 0; j0 < NLinks; j0++)
	{
		++BeginIndex[(inputLinks[j0].first)];
		++EndIndex[(inputLinks[j0].second)];
	}

	for(int i1 = 0; i1 < NHits; i1++)
	{
		if(BeginIndex[i1] == 0 && EndIndex[i1] == 1)
		{
			LeafHitsIndex.push_back(i1);
		}
		else if(BeginIndex[i1] == 1 && EndIndex[i1] == 1)
		{
			JointHitsIndex.push_back(i1);
		}
		else if(BeginIndex[i1] > 1 && EndIndex[i1] == 1)
		{
			StarJointHitsIndex.push_back(i1);
		}
		else if(BeginIndex[i1] == 1 && EndIndex[i1] == 0)
		{
			SimpleSeedHitsIndex.push_back(i1);
		}
		else if(BeginIndex[i1] > 1 && EndIndex[i1] == 0)
		{
			StarSeedHitsIndex.push_back(i1);
		}
		else if(BeginIndex[i1] == 0 && EndIndex[i1] == 0)
		{
			IsoHitsIndex.push_back(i1);
		}
		else
		{
			//edm::LogWarning("ArborWarning") <<"WARNING: UNCLASSIFIED HITS, Begin Index: "<<BeginIndex[i1]<<",  End Index:  "<<EndIndex[i1]<<endl; 
		}
	}
}

linkcoll LinkClean(const std::vector<ArborHit>& allhits, linkcoll alllinks )
{
	linkcoll cleanedlinks; 

	// int NLinks = alllinks.size();
	int Ncurrhitlinks = 0;
	int MinAngleIndex = -1;
	float MinAngle = 1E6;
	float tmpOrder = 0;
	float DirAngle = 0;

	std::pair<int, int> SelectedPair;

	TVector3 PosA, PosB, PosDiffAB;

	std::vector< std::vector<int> > LinkHits;
	LinkHits.clear();
	for(int s1 = 0; s1 < NHits; s1++) {
		std::vector<int> hitlink;

		auto range = BackLinksMap.equal_range(s1);
		for( auto itr = range.first; itr != range.second; ++itr ){
			hitlink.emplace_back(itr->second);
		}
		LinkHits.push_back(std::move(hitlink));
	}

	for(int i1 = 0; i1 < NHits; i1++) {
		PosB = cleanedHits[i1].GetPosition();
		MinAngleIndex = -10;
		MinAngle = 1E6;

		const std::vector<int>& currhitlink = LinkHits[i1];

		Ncurrhitlinks = currhitlink.size();

		for(int k1 = 0; k1 < Ncurrhitlinks; k1++)
		{
			PosA = cleanedHits[ currhitlink[k1] ].GetPosition();
			DirAngle = (PosA + PosB).Angle(PosB - PosA);		//Replace PosA + PosB with other order parameter ~ reference direction
			const float paddedAngle = (DirAngle + 0.1);
			tmpOrder = (PosB - PosA).Mag2() * paddedAngle * paddedAngle;
			if( tmpOrder < MinAngle ) // && DirAngle < 2.5 )
			{
				MinAngleIndex = currhitlink[k1];
				MinAngle = tmpOrder;
			}
		}

		if(MinAngleIndex > -0.5)
		{
			SelectedPair.first = MinAngleIndex;
			SelectedPair.second = i1;
			cleanedlinks.push_back(SelectedPair);
		}
	}
	// cout<<"CleanedLinks "<<cleanedlinks.size()<<endl;   
	return cleanedlinks;
}

void BuildInitLink()
{
	KDTreeLinkerAlgo<unsigned,3> kdtree;
	typedef KDTreeNodeInfoT<unsigned,3> KDTreeNodeInfo;

	Links.clear();	//all tmp links
	TVector3 PosDiffAB;
	std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
	std::vector<KDTreeNodeInfo> nodes, found;

	for(int i0 = 0; i0 < NHits; ++i0 ) 
	{     
		const auto& hit = cleanedHits[i0].GetPosition();
		nodes.emplace_back(i0,(float)hit.X(),(float)hit.Y(),(float)hit.Z());
		if( i0 == 0 ) 
		{
			minpos[0] = hit.X(); minpos[1] = hit.Y(); minpos[2] = hit.Z();
			maxpos[0] = hit.X(); maxpos[1] = hit.Y(); maxpos[2] = hit.Z();
		} 
		else 
		{
			minpos[0] = std::min((float)hit.X(),minpos[0]);
			minpos[1] = std::min((float)hit.Y(),minpos[1]);
			minpos[2] = std::min((float)hit.Z(),minpos[2]);
			maxpos[0] = std::max((float)hit.X(),maxpos[0]);
			maxpos[1] = std::max((float)hit.Y(),maxpos[1]);
			maxpos[2] = std::max((float)hit.Z(),maxpos[2]);
		}
	}
	KDTreeCube kdXYZCube(minpos[0],maxpos[0],
			minpos[1],maxpos[1],
			minpos[2],maxpos[2]);
	kdtree.build(nodes,kdXYZCube);
	nodes.clear();

	int NLayer_A = 0;
	int NLayer_B = 0; 
	int NStave_A = 0; 
	int NStave_B = 0; 

	const float InitLinkThreshold2 = std::pow(InitLinkThreshold,2.0);
	for(int i0 = 0; i0 < NHits; ++i0)
	{
		found.clear();
		const auto& PosA = cleanedHits[i0].GetPosition();
		NLayer_A = cleanedHits[i0].GetLayer();
		NStave_A = cleanedHits[i0].GetStave();

		const float side = InitLinkThreshold;
		const float xplus(PosA.X() + side), xminus(PosA.X() - side);
		const float yplus(PosA.Y() + side), yminus(PosA.Y() - side);
		const float zplus(PosA.Z() + side), zminus(PosA.Z() - side);
		const float xmin(std::min(xplus,xminus)), xmax(std::max(xplus,xminus));
		const float ymin(std::min(yplus,yminus)), ymax(std::max(yplus,yminus));
		const float zmin(std::min(zplus,zminus)), zmax(std::max(zplus,zminus));      
		KDTreeCube searchcube( xmin, xmax,
				ymin, ymax,
				zmin, zmax );
		kdtree.search(searchcube,found);
		for(unsigned j0 = 0; j0 < found.size(); ++j0)
		{
			if( found[j0].data <= (unsigned)i0 ) continue;
			const auto& PosB = cleanedHits[found[j0].data].GetPosition();
			NLayer_B = cleanedHits[found[j0].data].GetLayer();
			NStave_B = cleanedHits[found[j0].data].GetStave();

			PosDiffAB = PosA - PosB;

			if( PosDiffAB.Mag2() < InitLinkThreshold2 )
			{
				std::pair<int, int> a_Link;

				if( NStave_A != NStave_B || ( NLayer_A == 0 && NLayer_B != 0 ) || ( NLayer_B == 0 && NLayer_A != 0 ) )
				{
					if( PosA.Mag() > PosB.Mag() )
					{
						a_Link.first = found[j0].data;
						a_Link.second = i0;
					}
					else
					{
						a_Link.first = i0;
						a_Link.second = found[j0].data;
					}
					Links.push_back(a_Link);
					BackLinksMap.emplace(a_Link.second,a_Link.first);
				}
				else if( NLayer_A != NLayer_B && NStave_A == NStave_B )
				{
					if( NLayer_A > NLayer_B )
					{
						a_Link.first = found[j0].data;
						a_Link.second = i0;
					}
					else
					{
						a_Link.first = i0;
						a_Link.second = found[j0].data;
					}
					Links.push_back(a_Link);
					BackLinksMap.emplace(a_Link.second,a_Link.first);
				}

			}
		}
	}

	links_debug = Links; 
}

void LinkIteration()	//Energy corrections, semi-local correction
{
	KDTreeLinkerAlgo<unsigned,3> kdtree;
	typedef KDTreeNodeInfoT<unsigned,3> KDTreeNodeInfo;

	IterLinks.clear();
	IterBackLinks.clear();
	std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
	std::vector<KDTreeNodeInfo> nodes, found;

	alliterlinks = InitLinks;
	alliterBackLinksMap = BackLinksMap;
	int NInitLinks = InitLinks.size();
	int NLayer_A = 0; 
	int NLayer_B = 0;
	int NStave_A = 0; 
	int NStave_B = 0; 
	int AngleAccIndex = 0;
	int time = 1; 
	float MagA = 0;
	float MagB = 0;
	int FlagNoJoint = 0;
	int FlagNoIso = 0;
	TVector3 PosA, PosB, PosDiffAB, PosDiffBA, linkDir; 
	std::pair<int, int> currlink; 

	TVector3 RefDir[NHits];
	int Nin_hit[NHits];
	int Nout_hit[NHits];

	for(int i = 0; i < NHits; i++)
	{
		const auto& hit = cleanedHits[i].GetPosition();
		RefDir[i] = 1.0/hit.Mag() * hit;
		Nin_hit[i] = 0;
		Nout_hit[i] = 0;

		nodes.emplace_back(i,(float)hit.X(),(float)hit.Y(),(float)hit.Z());
		if( i == 0 ) 
		{
			minpos[0] = hit.X(); minpos[1] = hit.Y(); minpos[2] = hit.Z();
			maxpos[0] = hit.X(); maxpos[1] = hit.Y(); maxpos[2] = hit.Z();
		} 
		else 
		{
			minpos[0] = std::min((float)hit.X(),minpos[0]);
			minpos[1] = std::min((float)hit.Y(),minpos[1]);
			minpos[2] = std::min((float)hit.Z(),minpos[2]);
			maxpos[0] = std::max((float)hit.X(),maxpos[0]);
			maxpos[1] = std::max((float)hit.Y(),maxpos[1]);
			maxpos[2] = std::max((float)hit.Z(),maxpos[2]);
		}
	}

	KDTreeCube kdXYZCube(minpos[0],maxpos[0],
			minpos[1],maxpos[1],
			minpos[2],maxpos[2]);
	kdtree.build(nodes,kdXYZCube);
	nodes.clear();

	for(unsigned j = 0; j < (unsigned)NInitLinks; j++) 
	{
		currlink = InitLinks[j];
		PosA = cleanedHits[ currlink.first ].GetPosition();
		PosB = cleanedHits[ currlink.second ].GetPosition();
		linkDir = (PosA - PosB);		//Links are always from first point to second - verify
		linkDir *= 1.0/linkDir.Mag(); 
		RefDir[currlink.first] += 0.2*linkDir; 	//Weights... might be optimized...
		RefDir[currlink.second] += 0.2*linkDir; 
		Nin_hit[currlink.first] ++;
		Nout_hit[currlink.second] ++;
	}

	for(unsigned i1 = 0; i1 < (unsigned)NHits; i1++) 
	{
		found.clear();
		RefDir[i1] *= 1.0/RefDir[i1].Mag();
		PosA = cleanedHits[i1].GetPosition();
		NLayer_A = cleanedHits[i1].GetLayer();
		NStave_A = cleanedHits[i1].GetStave();

		const float side = IterLinkThreshold;
		const float xplus(PosA.X() + side), xminus(PosA.X() - side);
		const float yplus(PosA.Y() + side), yminus(PosA.Y() - side);
		const float zplus(PosA.Z() + side), zminus(PosA.Z() - side);
		const float xmin(std::min(xplus,xminus)), xmax(std::max(xplus,xminus));
		const float ymin(std::min(yplus,yminus)), ymax(std::max(yplus,yminus));
		const float zmin(std::min(zplus,zminus)), zmax(std::max(zplus,zminus));      
		KDTreeCube searchcube( xmin, xmax,
				ymin, ymax,
				zmin, zmax );
		kdtree.search(searchcube,found);

		for(unsigned j1 = 0; j1 < found.size(); j1++) 
		{
			if( found[j1].data <= i1 ) continue;

			FlagNoJoint = Nout_hit[j1] * Nin_hit[i1] * Nout_hit[i1] * Nin_hit[j1];
			FlagNoIso = Nout_hit[j1] + Nin_hit[i1] + Nout_hit[i1] + Nin_hit[j1];

			if( FlagNoJoint == 0 && FlagNoIso != 0 )
			{
				PosB = cleanedHits[found[j1].data].GetPosition();
				NLayer_B = cleanedHits[found[j1].data].GetLayer();
				NStave_B = cleanedHits[found[j1].data].GetStave();
				PosDiffAB = PosB - PosA;
				PosDiffBA = PosA - PosB;

				AngleAccIndex = 0;

				if( PosDiffAB.Angle(RefDir[i1]) < 0.6/time )
					AngleAccIndex = 1;
				else if( PosDiffAB.Angle(RefDir[j1]) < 0.6/time )
					AngleAccIndex = 2;
				else if( PosDiffBA.Angle(RefDir[i1]) < 0.6/time )
					AngleAccIndex = 3;
				else if( PosDiffBA.Angle(RefDir[j1]) < 0.6/time )
					AngleAccIndex = 4;

				if( PosDiffAB.Mag() > InitLinkThreshold && PosDiffAB.Mag() < IterLinkThreshold )  
				{

					MagA = PosA.Mag();
					MagB = PosB.Mag();

					if(NLayer_A != NLayer_B)
					{
						std::pair<int, int> a_Link;
						if( NStave_A != NStave_B || ( NLayer_A == 0 && NLayer_B != 0 ) || ( NLayer_B == 0 && NLayer_A != 0 ) )
						{
							if( MagA > MagB && AngleAccIndex < 3 )
							{
								a_Link.first = found[j1].data;
								a_Link.second = i1;
							}
							else if( MagA < MagB && AngleAccIndex > 2  )
							{
								a_Link.first = i1;
								a_Link.second = found[j1].data;
							}
							alliterlinks.push_back(a_Link);
							alliterBackLinksMap.emplace(a_Link.second,a_Link.first);
						}
						else if( NLayer_A != NLayer_B && NStave_A == NStave_B)
						{
							if( MagA > MagB && AngleAccIndex < 3 )
							{
								a_Link.first = found[j1].data;
								a_Link.second = i1;
							}
							else if( MagA < MagB && AngleAccIndex > 2 )
							{
								a_Link.first = i1;
								a_Link.second = found[j1].data;
							}
							alliterlinks.push_back(a_Link);
							alliterBackLinksMap.emplace(a_Link.second,a_Link.first);
						}
					}
				}
			}
		}
	}

	int MinAngleIndex = -10;
	int Ncurrhitlinks = 0; 
	float MinAngle = 1E6; 
	float tmpOrder = 0;
	float DirAngle = 0; 
	std::pair<int, int> SelectedPair; 

	std::vector< std::vector<int> > LinkHits;
	LinkHits.clear();
	for(int s1 = 0; s1 < NHits; s1++)
	{
		std::vector<int> hitlink;
		auto range = alliterBackLinksMap.equal_range(s1);
		for( auto itr = range.first; itr != range.second; ++itr )
		{
			hitlink.emplace_back(itr->second);
		}
		LinkHits.push_back(std::move(hitlink));
	}

	for(int i2 = 0; i2 < NHits; i2++)
	{
		PosB = cleanedHits[i2].GetPosition();
		MinAngleIndex = -10;
		MinAngle = 1E6;

		const std::vector<int>& currhitlink = LinkHits[i2];

		Ncurrhitlinks = currhitlink.size();

		for(int j2 = 0; j2 < Ncurrhitlinks; j2++)
		{
			PosA = cleanedHits[ currhitlink[j2] ].GetPosition();
			DirAngle = (RefDir[i2]).Angle(PosA - PosB);
			const float paddedAngle = (DirAngle + 1.0);
			tmpOrder = (PosB - PosA).Mag2() * paddedAngle * paddedAngle;
			if(tmpOrder < MinAngle) //  && DirAngle < 1.0)
			{
				MinAngleIndex = currhitlink[j2];
				MinAngle = tmpOrder;
			}
		}

		if(MinAngleIndex > -0.5)
		{
			SelectedPair.first = MinAngleIndex;
			SelectedPair.second = i2;
			IterLinks.push_back(SelectedPair);
			IterBackLinks.emplace(SelectedPair.second,SelectedPair.first);
		}
	}	
}

void BranchBuilding(float NeighbourThreshold) 
{

	int NLinks = IterLinks.size();
	int NBranches = 0;
	std::map <int, int> HitBeginIndex;
	std::map <int, int> HitEndIndex;
	std::vector< std::vector<int> > InitBranchCollection;
	std::vector< std::vector<int> > PrunedInitBranchCollection;
	std::vector< std::vector<int> > TmpBranchCollection;
	TVector3 PosA, PosB;

	for(int i1 = 0; i1 < NHits; i1++)
	{
		HitBeginIndex[i1] = 0;
		HitEndIndex[i1] = 0;
	}

	for(int j1 = 0; j1 < NLinks; j1++)
	{
		HitBeginIndex[ (IterLinks[j1].first) ] ++;
		HitEndIndex[ (IterLinks[j1].second) ] ++;
	}

	int iterhitindex = 0;
	int FlagInternalLoop = 0;
	int LL = 0; 	//Uplimt to be set to twice the total layer thickness...

	for(int i2 = 0; i2 < NHits; i2++)
	{
		if(HitEndIndex[i2] > 1)
			cout<<"WARNING OF INTERNAL LOOP with more than 1 link stopped at the same Hit"<<endl;

		if(HitBeginIndex[i2] == 0 && HitEndIndex[i2] == 1)        //EndPoint
		{
			NBranches ++;
			std::vector<int> currBranchhits;      //array of indexes 

			iterhitindex = i2;
			FlagInternalLoop = 0;
			currBranchhits.push_back(i2);
			LL = 0;

			while(FlagInternalLoop == 0 && HitEndIndex[iterhitindex] != 0 && LL < 300)	// 100 put by hand
			{
				auto iterlink_range = IterBackLinks.equal_range(iterhitindex);
				assert( std::distance(iterlink_range.first,iterlink_range.second) == 1 );
				iterhitindex = iterlink_range.first->second;
				currBranchhits.push_back(iterhitindex);
				LL++;
			}
			InitBranchCollection.push_back(std::move(currBranchhits) );
		}
	}

	PrunedInitBranchCollection.resize(InitBranchCollection.size());

	std::vector<float> BranchSize;
	std::vector<float> cutBranchSize;
	std::vector<int> SortedBranchIndex;
	std::vector<int> SortedcutBranchIndex;
	std::vector<int> currBranch;
	std::vector<int> iterBranch;
	std::vector<int> touchedHits;
	std::vector<bool> touchedHitsMap(NHits,false);
	std::vector<bool> seedHitsMap(NHits,false);
	std::vector<unsigned> seedHits;
	std::vector<int> leadingbranch;

	KDTreeLinkerAlgo<unsigned,3> kdtree;
	typedef KDTreeNodeInfoT<unsigned,3> KDTreeNodeInfo;
	std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
	std::vector<KDTreeNodeInfo> nodes, found;
	bool needInitPosMinMax = true;

	HitLinkMap seedToBranchesMap;
	std::unordered_map<int,int> branchToSeed;
	std::map<branch, int> SortedBranchToOriginal;
	SortedBranchToOriginal.clear();
	branchToSeed.clear();

	int currBranchSize = 0;
	int currHit = 0;

	for(int i3 = 0; i3 < NBranches; i3++)
	{
		currBranch = InitBranchCollection[i3];
		BranchSize.push_back( float(currBranch.size()) );
	}

	SortedBranchIndex = std::move( SortMeasure(BranchSize, 1) );

	for(int i4 = 0; i4 < NBranches; i4++)
	{
		currBranch = InitBranchCollection[SortedBranchIndex[i4]];
		currBranchSize = currBranch.size();
		iterBranch.clear();

		for(int j4 = 0; j4 < currBranchSize; j4++)
		{
			currHit = currBranch[j4];

			if( !touchedHitsMap[currHit] )
			{
				iterBranch.push_back(currHit);
				touchedHitsMap[currHit] = true;
			}
		}
		const auto theseed = currBranch[currBranchSize - 1];
		SortedBranchToOriginal[iterBranch] = theseed;     //Map to seed...
		branchToSeed.emplace(SortedBranchIndex[i4],theseed);
		seedToBranchesMap.emplace(theseed,SortedBranchIndex[i4]); // map seed to branches
		if( !seedHitsMap[theseed] ) 
		{
			seedHitsMap[theseed] = true;
			const auto& hit = cleanedHits[theseed].GetPosition();
			nodes.emplace_back(theseed,(float)hit.X(),(float)hit.Y(),(float)hit.Z());
			if( needInitPosMinMax ) {
				needInitPosMinMax = false;
				minpos[0] = hit.X(); minpos[1] = hit.Y(); minpos[2] = hit.Z();
				maxpos[0] = hit.X(); maxpos[1] = hit.Y(); maxpos[2] = hit.Z();
			} else {
				minpos[0] = std::min((float)hit.X(),minpos[0]);
				minpos[1] = std::min((float)hit.Y(),minpos[1]);
				minpos[2] = std::min((float)hit.Z(),minpos[2]);
				maxpos[0] = std::max((float)hit.X(),maxpos[0]);
				maxpos[1] = std::max((float)hit.Y(),maxpos[1]);
				maxpos[2] = std::max((float)hit.Z(),maxpos[2]);
			}
		}

		TmpBranchCollection.push_back(iterBranch);
		cutBranchSize.push_back( float(iterBranch.size()) );
		PrunedInitBranchCollection[SortedBranchIndex[i4]] = std::move(iterBranch);
	}

	SortedcutBranchIndex = std::move( SortMeasure(cutBranchSize, 1) );

	for(int i6 = 0; i6 < NBranches; i6++)
	{
		currBranch.clear();
		currBranch = TmpBranchCollection[ SortedcutBranchIndex[i6]];
		LengthSortBranchCollection.push_back(currBranch);;
	}

	std::vector<bool> link_helper(NBranches*NBranches,false);
	TVector3 DisSeed;
	KDTreeCube kdXYZCube(minpos[0],maxpos[0],
			minpos[1],maxpos[1],
			minpos[2],maxpos[2]);
	kdtree.build(nodes,kdXYZCube);
	nodes.clear();

	const float NeighbourThreshold2 = std::pow(NeighbourThreshold,2.0);

	QuickUnion qu(NBranches);

	for(int i7 = 0; i7 < NBranches; i7++)
	{
		auto SeedIndex_A = branchToSeed[i7];
		auto shared_branches = seedToBranchesMap.equal_range(SeedIndex_A);

		for( auto itr = shared_branches.first; itr != shared_branches.second; ++itr ) {
			const auto foundSortedIdx = itr->second;
			if( foundSortedIdx <= (unsigned)i7 ) continue;
			if( link_helper[NBranches*i7 + foundSortedIdx] ||
					link_helper[NBranches*foundSortedIdx + i7]    ) continue;
			if( !qu.connected(i7,foundSortedIdx) ) {
				qu.unite(i7,foundSortedIdx);
			}
			link_helper[NBranches*i7 + foundSortedIdx] = true;
			link_helper[NBranches*foundSortedIdx + i7] = true;
		}

		// do a kd tree search for seeds to merge together based on distance
		const auto& seedpos = cleanedHits[SeedIndex_A].GetPosition();
		const float side = NeighbourThreshold;
		const float xplus(seedpos.X() + side), xminus(seedpos.X() - side);
		const float yplus(seedpos.Y() + side), yminus(seedpos.Y() - side);
		const float zplus(seedpos.Z() + side), zminus(seedpos.Z() - side);
		const float xmin(std::min(xplus,xminus)), xmax(std::max(xplus,xminus));
		const float ymin(std::min(yplus,yminus)), ymax(std::max(yplus,yminus));
		const float zmin(std::min(zplus,zminus)), zmax(std::max(zplus,zminus));
		KDTreeCube searchcube( xmin, xmax,
				ymin, ymax,
				zmin, zmax );
		found.clear();
		kdtree.search(searchcube,found);
		for(unsigned j7 = 0; j7 < found.size(); j7++) {
			DisSeed = seedpos - cleanedHits[ found[j7].data ].GetPosition();

			if( DisSeed.Mag2() < NeighbourThreshold2 )
			{ 
				auto seed_branches = seedToBranchesMap.equal_range(found[j7].data);
				for( auto itr = seed_branches.first; itr != seed_branches.second; ++itr ){
					const auto foundSortedIdx = itr->second;
					if( foundSortedIdx <= (unsigned)i7 ) continue;
					if( link_helper[NBranches*i7 + foundSortedIdx] ||
							link_helper[NBranches*foundSortedIdx + i7]    ) continue;
					if( !qu.connected(i7,foundSortedIdx) ) {
						qu.unite(i7,foundSortedIdx);
					}
					link_helper[NBranches*i7 + foundSortedIdx] = true;
					link_helper[NBranches*foundSortedIdx + i7] = true;
				}
			}
		}
	}

	Trees.clear();
	std::unordered_map<unsigned,branch> merged_branches(qu.count());
	Trees.reserve(qu.count());

	for( unsigned i = 0; i < (unsigned)NBranches; ++i ) {
		unsigned root = qu.find(i);
		const auto& branch = PrunedInitBranchCollection[i];
		auto& merged_branch = merged_branches[root];
		merged_branch.insert(merged_branch.end(),branch.begin(),branch.end());
	}

	unsigned total_hits = 0;
	for( auto& final_branch : merged_branches ) {
		total_hits += final_branch.second.size();
		Trees.push_back(std::move(final_branch.second));
	}
}

void BushMerging()
{
}

void BushAbsorbing()
{
}

std::vector< std::vector<int> > Arbor(std::vector<ArborHit> inputHits, 
		const float CellSize, 
		const float LayerThickness)
{
	init(CellSize, LayerThickness);

	HitsCleaning( std::move(inputHits) );
	BuildInitLink();

	InitLinks = std::move( LinkClean( cleanedHits, Links ) );
	cout<<"InitLinks: "<<InitLinks.size()<<endl; 

	LinkIteration();
	HitsClassification(IterLinks);	
	cout<<"IterLinks: "<<IterLinks.size()<<endl; 

	BranchBuilding( 2.0*CellSize );	

	return Trees;
}

