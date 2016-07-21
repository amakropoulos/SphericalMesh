
#include "mirtk/Common.h"
#include "mirtk/M2SConnectivity.h"
#include "mirtk/PointSetUtils.h"


#include <vtkMassProperties.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkFastMarchingGeodesicDistance.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>

#include <boost/range/irange.hpp>

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================


M2SConnectivity::M2SConnectivity()
{
  _Input = vtkSmartPointer<vtkPolyData>::New();
  _Normalize = false;
  _ScaleFactor = -1.0;
  _Threshold = -1.0;
  _K = -1;

  _ForwardDirection = true;
  _BackwardDirection = true;
  _IdentityEdges = false;
}


// -----------------------------------------------------------------------------


M2SConnectivity::M2SConnectivity(const M2SConnectivity &other)
{
  *this = other;
}


// -----------------------------------------------------------------------------

M2SConnectivity &M2SConnectivity::operator =(const M2SConnectivity &other)
{
  if (this != &other) {
    M2SConnectivity::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

M2SConnectivity::~M2SConnectivity(){}

// -----------------------------------------------------------------------------

void M2SConnectivity::CopyAttributes(const M2SConnectivity &other){
  cout << "TODO copy attributes" << endl;
}

void progress(float p){
  int barWidth = 70;

  cout << "[";
  int pos = barWidth * p;
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(p * 100.0) << " %\r";
  std::cout.flush();


}






bool idPairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) {

  if (firstElem.first < secondElem.first){
    return true;
  }

  if (firstElem.first > secondElem.first){
    return false;
  }

  return firstElem.second < secondElem.second;
}

bool compareIdsLessThan(const Edge &l, const Edge &r){
  return l.ptId1 < r.ptId1 || ( l.ptId1 == r.ptId1 &&  l.ptId2 < r.ptId2);
};

bool compareIdsGreaterThan (const Edge &l, const Edge &r){
  return l.ptId1 > r.ptId1 || ( l.ptId1 == r.ptId1 &&  l.ptId2 > r.ptId2);
};

bool compareIdsEquivalent (const Edge &l, const Edge &r){
  return l.ptId1 == r.ptId1 && l.ptId2 == r.ptId2;
};

bool compareIdsNotEquivalent (const Edge &l, const Edge &r){
  return l.ptId1 != r.ptId1 || l.ptId2 != r.ptId2;
};

void uniqueSortedEdges(Array<Edge> & edges){
  int numEdges = edges.size();

  sort(edges.begin(), edges.end(), compareIdsLessThan);

  Array<Edge> uniqueEdges;

  int count = 0;
  double meanLength = 0.0;
  int lastEdgeId = numEdges -1;

  for (int edgeId=0; edgeId<numEdges; edgeId++){
    Edge edge = edges[edgeId];
    count++;

    Edge nextEdge;
    if (edgeId < lastEdgeId){
      nextEdge = edges[edgeId+1];
    }

    if (compareIdsNotEquivalent(edge, nextEdge) || edgeId == lastEdgeId){
      meanLength += edge.length;
      meanLength /= (double)count;
      edge.length = meanLength;
      uniqueEdges.push_back(edge);

      count = 0;
      meanLength = 0.0;
    }
    else{
      meanLength += edge.length;
    }
  }

  edges = uniqueEdges;
}

void M2SConnectivity::ComputeGeodesics(){

  int numPts = _Input->GetNumberOfPoints();

  string arrayName = "FMMDist";
  vtkSmartPointer< vtkFastMarchingGeodesicDistance > fmm =
      vtkSmartPointer< vtkFastMarchingGeodesicDistance >::New();

  vtkSmartPointer< vtkIdList > seed =
      vtkSmartPointer< vtkIdList >::New();

  int numSeeds = _SeedPtIds.size();
  if (numSeeds == 0){
    _SeedPtIds = Array<int>(numPts);
    std::iota (_SeedPtIds.begin(), _SeedPtIds.end(), 0);
    numSeeds = numPts;
  }

  int numOtherPts = _OtherPtIds.size();
  if (numOtherPts == 0){
    _OtherPtIds = Array<int>(numPts);
    std::iota (_OtherPtIds.begin(), _OtherPtIds.end(), 0);
    numOtherPts = numPts;
  }

  fmm->SetInputData(_Input);
  fmm->SetFieldDataName(arrayName.c_str());
  seed->SetNumberOfIds(1);
  fmm->SetSeeds(seed);
  fmm->SetDistanceStopCriterion(_Threshold);

  if (_K > 0){
    int k = _K;
    if (!_IdentityEdges)
      k++;
    fmm->SetNumberOfVisitedPointsStopCriterion(k);
  }


  int onePercentOfSeeds = numSeeds/100;
  onePercentOfSeeds++;

  cout << "computing geodesics" << endl;

  bool symmetric = _ForwardDirection && _BackwardDirection;


  for (int i=0; i<numSeeds; i++){
    int seedPtId = _SeedPtIds[i];
    //cout << ptId << endl;
    if (i % onePercentOfSeeds == 0){
      double p = double(i) / double(numSeeds);
      progress(p);
    }

    seed->SetId(0, seedPtId);
    seed->Modified();
    fmm->Update();


    vtkDoubleArray * array = (vtkDoubleArray *)
        fmm->GetOutput()->GetPointData()->GetArray(arrayName.c_str());

    //note otherPtId < ptId, as int otherPtId=ptId+1 below


    int numEdge = 0;
    for (int j=0; j<numOtherPts; j++){

      int otherPtId = _OtherPtIds[j];

      if ((!_IdentityEdges) && seedPtId == otherPtId)
        continue;

      double d = array->GetTuple1(otherPtId);
      if (d >= 0.0 && (_Threshold <= 0.0 || d < _Threshold)){

        int ptId1 = seedPtId;
        int ptId2 = otherPtId;

        if (symmetric && otherPtId < seedPtId){
          std::swap(ptId1, ptId2);
        }

        Edge e;
        e.ptId1 = ptId1;
        e.ptId2 = ptId2;
        e.length = d;
        e.weight1 = _ForwardDirection;
        e.weight2 = _BackwardDirection;
        _Edges.push_back(e);
      }
    }

  }

  if (symmetric){ uniqueSortedEdges(_Edges); }

}


/*
void M2SConnectivity::InitializeEdges(){

  //get edges

  //map original id to subset id
  //map submesh id to original id -> ptidhierarchy

  _Edges.clear();

  cout << "Initialising edges" << endl;
  Array< Pair<int, int> > edgePairs;

  int numPts = _Input->GetNumberOfPoints();


  GenericSparseMatrix<double>::Entries entries;

  for (int ptId1 = 0; ptId1< numPts; ptId1++){

    if ( ptId1 >= numPts){
      cout << "InitializeEdges::point id ptId1 " << ptId1 << " is greater than number of points! " << numPts << endl;
      exit(1);
    }

    _Geodesics.GetCol(ptId1, entries);

    int numAdjs = entries.size();

    for (int k=0; k<numAdjs; k++){
      int ptId2 = entries[k].first;


      if ( ptId2 >= numPts){
        cout << "InitializeEdges::point id ptId2 is greater than number of points!" << endl;
        exit(1);
      }

      if (ptId1 < ptId2){
        double d = _Geodesics.Get(ptId2,ptId1);

        if (d <= 0){
          cout << "M2SConnectivity::InitialiseEdges : invalid geodesic value "
              << d << endl;
          exit(1);
        }
        edgePairs.push_back(MakePair(ptId1, ptId2));
      }

      else{
        cout << "M2SConnectivity::InitialiseEdges : invalid id pair "
            << "<" << ptId1 << "," << ptId2 << ">" << endl;
      }
    }
  }
  //add mesh links if required
  if (_EdgeLinksOn){

    _Input->GetPolys()->InitTraversal();
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();

    while(_Input->GetPolys()->GetNextCell(idList)) {
      int numIds = idList->GetNumberOfIds();

      for (int i=0; i<numIds ;i++){
        int ptId1 = idList->GetId(i);

        for (int j=0; j<numIds ;j++){
          int ptId2 = idList->GetId(j);

          if (ptId1 < ptId2){
            edgePairs.push_back(MakePair(ptId1,ptId2));
          }
        }
      }
    }
  }
  //get only unique pairs
  uniqueSortedPairs(edgePairs);

  _Edges.Initialize(edgePairs, numPts);
}
*/
void M2SConnectivity::ScaleEdgeDistances(){

  //radius of sphere with same surface area gives normalisation to unit sphere

  if (_ScaleFactor < 0){
    vtkSmartPointer<vtkMassProperties> mp =
        vtkSmartPointer<vtkMassProperties>::New();
    mp->SetInputData(_Input);
    mp->Update();
    double area = mp->GetSurfaceArea();
    _ScaleFactor = sqrt(area / (4 * pi));
  }


  for (auto &edge:_Edges){
    edge.length /= _ScaleFactor;
  }

}


/*
void M2SConnectivity::InitializeEdgeDistances(){
  //extract vector of geodesics
  EdgeIterator ei(_Edges);
  ei.InitTraversal();
  vtkIdType pid1, pid2;

  _Distances.clear();
 //double maxd = 0;
  while ( ei.GetNextEdge(pid1, pid2) != -1 ){
    double d = _Geodesics.Get(pid1,pid2);

    if (d == 0){
      d = _Geodesics.Get(pid2,pid1);
    }

    if (d == 0){//not in geodesics must be mesh edge compute length

      //cout << "WARNING - distance not in geodesics" << endl;

      double p1 [3], p2 [3];

      _Input->GetPoint(pid1,p1);
      _Input->GetPoint(pid2,p2);
      //cout << "edge link found computing distance" << endl;
      d = sqrt(pow(p1[0]-p2[0],2)
                  + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2));

      //if (d > maxd) maxd = d;
    }
    _Distances.push_back(d);
  }
  //cout << "maxd " << maxd << endl;
  ScaleEdgeDistances();
}
*/

/*
void M2SConnectivity::ComputeEdgeLinks(){

  int numTargetPts = _TargetMesh->GetNumberOfPoints();


  EdgeTable edgeTable(_TargetMesh);
  EdgeIterator ei(edgeTable);

  Array<double> weights;

  for (int ptId=0; ptId<numTargetPts; ptId++){
    int numAdjs = edgeTable.NumberOfAdjacentPoints(ptId);
    double w = 1.0 / (double)numAdjs;
    weights.push_back(w);
  }

  ei.InitTraversal();
  int ptId1, ptId2;
  while ( ei.GetNextEdge(ptId1, ptId2) != -1 ){

    Pair idPair = MakePair(ptId1, ptId2);
    bool seen = _VisitedEdges.find(idPair) != _VisitedEdges.end();
    if (seen)
      continue;

    _VisitedEdges.insert(idPair);

    double p1 [3], p2 [3];
    _TargetMesh->GetPoint(ptId1,p1);
    _TargetMesh->GetPoint(ptId2,p2);
    //cout << "edge link found computing distance" << endl;
    double d = sqrt(pow(p1[0]-p2[0],2)
                + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2));

    Edge e;
    e.ptId1 = ptId1;
    e.ptId2 = ptId2;
    e.length = d;
    e.weight1 = weights[ptId1];
    e.weight2 = weights[ptId2];

    _Edges.push_back(e);
  }
}

*/
/*
void M2SConnectivity::GetEdges(Array<int> seedPtIds, Array<Edge> otherPtIds, int max, double maxDist = -1.0){

  int maxId = 0;

  for (auto e:_Edges){
    if (e.ptId1 > maxId){ maxId = e.ptId1; }
    if (e.ptId2 > maxId){ maxId = e.ptId2; }
  }

  Array<bool> seedPtMask;
  if (seedPtIds.empty()){
    seedPtMask = Array<bool>(maxId+1, true);
  }
  else{
    seedPtMask = Array<bool>(maxId+1, false);
    for (auto id:seedPtIds){
      seedPtMask[id] = true;
    }
  }

  Array<bool> otherPtMask;
  if (otherPtIds.empty()){
    otherPtMask = Array<bool>(maxId+1, true);
  }
  else{
    otherPtMask = Array<bool>(maxId+1, false);
    for (auto id:otherPtIds){
      otherPtMask[id] = true;
    }
  }

}
*/

void M2SConnectivity::Update(){
  ComputeGeodesics();
  if (_Normalize)
    ScaleEdgeDistances();

}



} // mirtk namespace

