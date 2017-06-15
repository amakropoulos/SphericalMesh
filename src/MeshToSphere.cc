#include "mirtk/MeshToSphere.h"

//mirtk includes
#include "mirtk/Common.h"
// #include "mirtk/DeformableSurfaceModel.h"
// #include "mirtk/InflationForce.h"
// #include "mirtk/InflationStoppingCriterion.h"
#include "mirtk/LocalOptimizer.h"
// #include "mirtk/EulerMethodWithMomentum.h"
// #include "mirtk/MetricDistortion.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/M2SRemesher.h"
#include "mirtk/SparseMatrix.h"

//my mirtk includes
#include "mirtk/M2SConnectivity.h"
#include "mirtk/M2SOptimizer.h"
#include "mirtk/M2SDiffuser.h"
#include "mirtk/M2SPostProcessor.h"
// #include "mirtk/M2SEdge.h"
#include "mirtk/M2SMeshMerger.h"
#include "mirtk/M2SGeodesicInterpolator.h"


//vtk includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkMassProperties.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkFastMarchingGeodesicDistance.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================

MeshToSphere::MeshToSphere(){
  _FinalInterpolation = false;
  _FinalInterpolationThreshold = -1.0;

  _PostProcessingIterations = 10;
  _PostProcessingPercentile = 99;
  _PostProcessingThreshold = -1.0;
  _LocalWeighting = 1.0;

  Initialize();
}


// -----------------------------------------------------------------------------


MeshToSphere::MeshToSphere(const MeshToSphere &other)
{
  *this = other;
}

// -----------------------------------------------------------------------------

MeshToSphere &MeshToSphere::operator =(const MeshToSphere &other)
{
  if (this != &other) {
    MeshToSphere::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeshToSphere::~MeshToSphere(){
}

void MeshToSphere::CopyAttributes(const MeshToSphere &other){
  cout << "TODO copy attributes" << endl;
}


// =============================================================================
// Class methods
// =============================================================================

void MeshToSphere::Initialize(){
  _Input = vtkSmartPointer<vtkPolyData>::New();
  _Output = vtkSmartPointer<vtkPolyData>::New();
  _Inflated = NULL;
}

// -----------------------------------------------------------------------------
/// Compute centroid of point set
void GetCentroid(vtkPolyData *mesh, double centroid[3])
{
  double p[3];
  centroid[0] = centroid[1] = centroid[2] = .0;
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    mesh->GetPoint(ptId, p);
    centroid[0] += p[0];
    centroid[1] += p[1];
    centroid[2] += p[2];
  }
  centroid[0] /= npoints;
  centroid[1] /= npoints;
  centroid[2] /= npoints;
}

// -----------------------------------------------------------------------------
/// Get approximate scale of point set
void GetScale(vtkPolyData *mesh, const double centroid[3], double scale[3])
{
  double p[3];
  scale[0] = scale[1] = scale[2] = .0;
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    mesh->GetPoint(ptId, p);
    scale[0] += abs(p[0] - centroid[0]);
    scale[1] += abs(p[1] - centroid[1]);
    scale[2] += abs(p[2] - centroid[2]);
  }
  scale[0] /= npoints;
  scale[1] /= npoints;
  scale[2] /= npoints;
}

void MeshToSphere::EmbedPointsOnSphere(vtkPolyData * mesh){

  int numPoints = mesh->GetNumberOfPoints();
  double p [3];

  for (int ptId = 0; ptId < numPoints; ++ptId) {
    mesh->GetPoint(ptId, p);

    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

    if (norm > 0.0){
      for (int i = 0; i < 3; ++i){ p[i] /= norm; }
      mesh->GetPoints()->SetPoint(ptId, p);
    }
  }


}

void Smooth(vtkSmartPointer<vtkPolyData> mesh){

  vtkSmartPointer<vtkSmoothPolyDataFilter> smoother =
      vtkSmartPointer<vtkSmoothPolyDataFilter>::New();

  smoother->SetInputData(mesh);
  smoother->SetRelaxationFactor(1.0);
  smoother->SetFeatureEdgeSmoothing(0);
  smoother->SetNumberOfIterations(1);
  smoother->SetConvergence(0);
  mesh = smoother->GetOutput();

}

void MeshToSphere::LowPassFilter(vtkPolyData * mesh){
  //TODO investigate pass bands
  double lowPassBand = 0.5;
  int lowPassIterations = 10;

  double c1[3], s1[3];
  GetCentroid(mesh, c1);
  GetScale(mesh, c1, s1);

  // Perform low-pass filtering
     //
     // The translation and scale "fix" is due to a bug in vtkWindowedSincPolyDataFilter:
     // http://vtk.1045678.n5.nabble.com/Bug-in-vtkWindowedSincPolyDataFilter-td1234055.html

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> filter =
      vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  filter->SetPassBand(lowPassBand);
  filter->SetNumberOfIterations(lowPassIterations);
  filter->NormalizeCoordinatesOn();
  SetVTKInput(filter, mesh);
  filter->Update();

  double c2[3], s2[3];
  GetCentroid(filter->GetOutput(), c2);
  GetScale(filter->GetOutput(), c2, s2);

  // Adjust point displacements
  double p1[3], p2[3], sf[3];

  for (int i=0; i<0; i++) { sf[i] = s2[i] / s1[i]; } //approximate scale factor

  for (vtkIdType ptId = 0; ptId < mesh->GetNumberOfPoints(); ++ptId) {
    mesh->GetPoint(ptId, p1);
    filter->GetOutput()->GetPoint(ptId, p2);

    for (int i=0; i<0; i++){
      p2[i] = ((p2[i] - c2[i]) / sf[i]) + c1[i];
    }

    mesh->GetPoints()->SetPoint(ptId, p2);
  }

}


/*

void MeshToSphere::LowPassFilter(vtkPolyData * mesh){

  LowPassFilterSlow(mesh);
  return;



  double lowPassBand = 0.01;
  int lowPassIterations = 10;

  double p0[3], p1[3], c0[3], c1[3], scale[3];

  // original scale of spherical mesh is 1,1,1 as all points lie on the surface r=1
  // centriod is not 0,0,0 however due to uneven distribution of points


  int numActivePts = 0; // _PtIdHierarchy[_CurrentLevel].size();

  c0[0] = c0[1] = c0[2] = .0;
  int numPts = _Output->GetNumberOfPoints();
  for (int ptId = 0; ptId < numPts; ptId++) {
    if (_ActivePoints[ptId]){
      _Output->GetPoint(ptId, p0);
      c0[0] += p0[0];
      c0[1] += p0[1];
      c0[2] += p0[2];
      numActivePts++;
    }
  }
  c0[0] /= numActivePts;
  c0[1] /= numActivePts;
  c0[2] /= numActivePts;

  // Perform low-pass filtering
     //
     // The translation and scale "fix" is due to a bug in vtkWindowedSincPolyDataFilter:
     // http://vtk.1045678.n5.nabble.com/Bug-in-vtkWindowedSincPolyDataFilter-td1234055.html
  vtkSmartPointer<vtkPolyData> relaxedMesh =
      vtkSmartPointer<vtkPolyData>::New();

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> filter =
      vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  filter->SetPassBand(lowPassBand);
  filter->SetNumberOfIterations(lowPassIterations);
  //filter->NormalizeCoordinatesOn();
  SetVTKInput(filter, _Output);
  filter->Update();
  relaxedMesh = filter->GetOutput();


  c1[0] = c1[1] = c1[2] = .0;
  for (int ptId = 0; ptId < numPts; ++ptId) {
    if (_ActivePoints[ptId]){
      relaxedMesh->GetPoint(ptId, p1);
      c1[0] += p1[0];
      c1[1] += p1[1];
      c1[2] += p1[2];
    }
  }
  c1[0] /= numActivePts;
  c1[1] /= numActivePts;
  c1[2] /= numActivePts;


  scale[0] = scale[1] = scale[2] = .0;
  for (int ptId = 0; ptId < numPts; ++ptId) {
    if (_ActivePoints[ptId]){
      relaxedMesh->GetPoint(ptId, p1);
      scale[0] += abs(p1[0] - c1[0]);
      scale[1] += abs(p1[1] - c1[1]);
      scale[2] += abs(p1[2] - c1[2]);
    }
  }
  scale[0] /= numActivePts;
  scale[1] /= numActivePts;
  scale[2] /= numActivePts;


  // Adjust point displacements

  for (int ptId = 0; ptId < numPts; ptId++) {
    if (_ActivePoints[ptId]){
      _Output->GetPoint(ptId, p0);
      relaxedMesh->GetPoint(ptId, p1);

      for (int i=0; i<0; i++){
        p1[i] = ((p1[i] - c1[i]) / scale[i]) + c0[i];
      }
      _Output->GetPoints()->SetPoint(ptId, p1);
    }
  }


}
*/

void MeshToSphere::Remesh(double minEdgeLength, double maxEdgeLength, bool mapPoints){

  int maxIterations = 100;
  int numOriginalPoints = _Input->GetNumberOfPoints();

  Remesher r;
  if (_Inflated != NULL){
    cout << "Using inflated mesh for to aid decimation" << endl;
    r.Input(_Inflated);
  }
  else{
    r.Input(_Input);
  }

  r.MinEdgeLength(minEdgeLength);
  r.MaxEdgeLength(maxEdgeLength);
  r.MaxIterations(maxIterations);
  r.Verbose(true);
  r.Run();

  vtkSmartPointer<vtkPolyData> newMesh =
      vtkSmartPointer<vtkPolyData>::New();
  newMesh->DeepCopy(r.Output());

  int numNewPts = newMesh->GetNumberOfPoints();

  vtkSmartPointer<vtkKdTreePointLocator> tree
      = vtkSmartPointer<vtkKdTreePointLocator>::New();

  if (_Inflated != NULL){
    cout << "Using inflated mesh to lookup original points" << endl;
    tree->SetDataSet(_Inflated);
  }
  else{
    tree->SetDataSet(_Input);
  }

  tree->BuildLocator();

  double np [3], p [3];
  Array<int> mapping;




/*
// Original code by Rob
  for (int nid=0; nid<numNewPts; nid++) {
    //replace remeshed point with nearest mesh point
    newMesh->GetPoint(nid, np);
    int id = tree->FindClosestPoint(np);
    _Input->GetPoint(id, p);

    if (mapPoints){
      newMesh->GetPoints()->SetPoint(nid, p);
    }

    mapping.push_back(id); // also keep mesh id
  }

  _MeshHierarchy.push_back(newMesh);
  _PtIdHierarchy.push_back(mapping);

  OrderedSet<int> seen;
  for (auto id: mapping){
    if (seen.find(id) != seen.end()){
      cout << "ERROR: MeshToSphere::Remesh - duplicate id in remeshing" << endl;
      cout << "TODO: modify code so that each vertex is assigned a unique vertex in the parent mesh" << endl;
      exit(1);
    }
    seen.insert(id);
  }
  */




  int nid, id;
  double dist, maxdist=-1;
  OrderedMap< double , Pair<int, int> > mappingDist;
  for (nid=0; nid<numNewPts; nid++) {
    //replace remeshed point with nearest mesh point
    newMesh->GetPoint(nid, np);
    id = tree->FindClosestPoint(np);

    if (_Inflated != NULL){
      _Inflated->GetPoint(id, p);
    }
    else{
      _Input->GetPoint(id, p);
    }

    dist = sqrt(vtkMath::Distance2BetweenPoints(p, np));

    Pair<int,int> pair = MakePair(nid, id);
    mappingDist[dist] = pair;
    mapping.push_back(id);

    if(dist>maxdist) maxdist=dist;

  }

  OrderedSet<int> seen;
  for(auto iter: mappingDist){
    nid=iter.second.first;
    id=iter.second.second;

    // if we have a point in the old mesh already assigned 
    // find iteratively points connected to the original 
    // and use the first unassigned with the closest distance
    if (seen.find(id) != seen.end()){
      newMesh->GetPoint(nid, np);
      dist=iter.first;

      int pointid, nnpointid;
      double pointdist, nnpointdist;
      bool found = false;

      OrderedMap<double, int > nnpoints;
      OrderedSet<int> visited;
      nnpoints[dist] = id;

      while(nnpoints.size()>0){
        auto it=nnpoints.begin();
        pointdist = it->first;
        pointid = it->second;
        nnpoints.erase(it);    
        if(visited.find(pointid) != visited.end()) continue;
        visited.insert(pointid);    

        if(seen.find(pointid) != seen.end()){
          vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
          _Input->GetPointCells(pointid, cellIds);
          for (int i=0; i<cellIds->GetNumberOfIds(); i++){
            vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
            _Input->GetCellPoints(cellIds->GetId(i), pointIds);
            for (int j=0; j<pointIds->GetNumberOfIds(); j++){
              nnpointid=pointIds->GetId(j);
              if(visited.find(nnpointid) != visited.end()) continue;
              _Input->GetPoint(nnpointid, p);
              nnpointdist = sqrt(vtkMath::Distance2BetweenPoints(p, np));
              nnpoints[nnpointdist] = nnpointid;
            }
          }
        }else{
          found = true;
          break;
        }
      }


      // The distance must be within the range of the existing distances
      if( found && pointdist <= 1.5*maxdist ){
        cout<<"replaced id due to duplicate for point in new mesh "<<nid<<": point in old mesh "<<id<<"->"<<pointid<<",  distance "<<dist<<"->"<<pointdist<<endl;
        cout<<"maxdist: "<<maxdist<<endl;
        id = pointid;
        mapping[nid]=id;
      }
      // otherwise fail
      else{
        cout << "ERROR: MeshToSphere::Remesh - duplicate id in remeshing" << endl;
        cout << "TODO: modify code so that each vertex is assigned a unique vertex in the parent mesh" << endl;
        exit(1);
      }
    }

    seen.insert(id);
  }

  //change points of new mesh back to original as opposed to inflated
  if (_Inflated != NULL){
    for (int nid=0; nid<numNewPts; nid++) {
        id = mapping[nid];
        _Input->GetPoint(id, p);
        newMesh->GetPoints()->SetPoint(nid, p);
    }
  }



  _MeshHierarchy.push_back(newMesh);
  _PtIdHierarchy.push_back(mapping);


  if (mapPoints){
    for (int nid=0; nid<numNewPts; nid++) {
        id = mapping[nid];
        _Input->GetPoint(id, p);
        newMesh->GetPoints()->SetPoint(nid, p);
    }
  }

}






void MeshToSphere::ComputeMeshHierarchy(){

  _MeshHierarchy.clear();

  int l=0;
  for (auto p: _Parameters){
    l++;

    double minEdgeLength = p.remesh.minEdgeLength.value;
    double maxEdgeLength = p.remesh.maxEdgeLength.value;

    if (minEdgeLength != -1 && maxEdgeLength != -1){
      Remesh(minEdgeLength, maxEdgeLength, p.diffusion.active);
      cout <<  "level " << l << "," << _PtIdHierarchy[l].size() << "points" << endl;
    }

    else{
      cout << "pushing back full mesh" <<  endl;
      vtkSmartPointer<vtkPolyData> inputCopy =
          vtkSmartPointer<vtkPolyData>::New();
      inputCopy->DeepCopy(_Input);
      _MeshHierarchy.push_back(inputCopy);

      Array<int> ptIds;
      for (int ptId=0; ptId<_Input->GetNumberOfPoints(); ptId++){
        ptIds.push_back(ptId);
      }
      _PtIdHierarchy.push_back(ptIds);
    }

    if (_Debug) {
        string outName = "level" + std::to_string(l) + "-input-mesh.vtk";
        WritePolyData(outName.c_str(), _MeshHierarchy[_MeshHierarchy.size()-1]);
    }
  }
}











void MeshToSphere::InterpolateAngles(int previousLevel, int nextLevel = -1){

  int diffusionIterations = 1000;

  //compute edge lengths

  int numInputPts = _Input->GetNumberOfPoints();

  //initialise source mask
  Array<bool> sourceMask;
  for (int i=0; i<numInputPts; i++){
    sourceMask.push_back(false);
  }

  Array<int> sourcePtIds = _PtIdHierarchy[previousLevel];
  int numSourcePts = sourcePtIds.size();

  for (int i=0; i<numSourcePts; i++){
    int ptId = sourcePtIds[i];
    sourceMask[ptId] = true;
  }

  //initialise full mesh onto sphere
  vtkSmartPointer<vtkPolyData> fullOutput =
      vtkSmartPointer<vtkPolyData>::New();

  fullOutput->DeepCopy(_Input);
  for (int ptId=0; ptId<numInputPts; ptId++){
    double p[3] = {0.0,0.0,0.0};
    fullOutput->GetPoints()->SetPoint(ptId,p);
  }

  if (_Debug){
    string outName = "level" + std::to_string(previousLevel) + "_pre_diffusion.vtk";
    WritePolyData(outName.c_str(), _Output);
  }

  for (int i=0; i<numSourcePts; i++){
    double p[3];
    _Output->GetPoint(i,p);
    int ptId =  sourcePtIds[i];
    fullOutput->GetPoints()->SetPoint(ptId,p);
  }

  _Output = fullOutput;





  //diffuse coordinates to non source points over sphere
  M2SDiffuser diffuser;
  diffuser.Input(_Output);
  diffuser.SourceMask(sourceMask);
  diffuser.MaxIterations(diffusionIterations);
  diffuser.AdaptiveTermination(true);
  //diffuser.Debug(_Debug);
  diffuser.Run();
  _Output = diffuser.Output();
  EmbedPointsOnSphere(_Output);

  if (_Debug){
    string outName = "level" + std::to_string(previousLevel) + "_post_diffusion.vtk";
    WritePolyData(outName.c_str(), _Output);
  }

  //retriangulate discontinuities areas



  if (nextLevel != -1){
    vtkPolyData * nextMesh = _MeshHierarchy[nextLevel];

    //initialise full mesh onto sphere
    vtkSmartPointer<vtkPolyData> nextSphere =
        vtkSmartPointer<vtkPolyData>::New();

    nextSphere->DeepCopy(nextMesh);

    for (int i=0; i<nextMesh->GetNumberOfPoints(); i++){
      int ptId = _PtIdHierarchy[nextLevel][i];
      double p[3];
      _Output->GetPoint(ptId,p);
      nextSphere->GetPoints()->SetPoint(i,p);
    }

    _Output = nextSphere;
  }

  if (_Debug){
    string outName = "level" + std::to_string(previousLevel) + "_downsampled.vtk";
    WritePolyData(outName.c_str(), _Output);
  }

}






void MeshToSphere::InitialiseSphericalMesh(){


  cout << "out has " << _MeshHierarchy[0]->GetNumberOfPoints() << " potins" << endl;

  _Output->DeepCopy(_MeshHierarchy[0]);
  cout << "out has " << _Output->GetNumberOfPoints() << " potins" << endl;

  Center(_Output);
  EmbedPointsOnSphere(_Output);
  for (int i=0; i<10; i++){
    LowPassFilter(_Output);
    EmbedPointsOnSphere(_Output);
  }

}



Array<Edge> ComputeNeighbourEdges(vtkPolyData * mesh, double weightScaleFactor = 1.0){

  Array<Edge> out;

  EdgeTable edges(mesh);

  Array<double> edgeWeights;
  int numPts = mesh->GetNumberOfPoints();
  for (int ptId = 0; ptId < numPts; ptId++){
    double w = 1.0 / edges.NumberOfAdjacentPoints(ptId);
    w *= weightScaleFactor;
    edgeWeights.push_back(w);
  }

  EdgeIterator ei(edges);
  ei.InitTraversal();
  int ptId1, ptId2;
  while (ei.GetNextEdge(ptId1, ptId2) != -1){
    Edge e;
    e.ptId1 = ptId1;
    e.ptId2 = ptId2;

    double pt1 [3], pt2 [3];
    mesh->GetPoint(ptId1, pt1);
    mesh->GetPoint(ptId2, pt2);

    e.length = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
    e.weight1 = edgeWeights[ptId1];
    e.weight2 = edgeWeights[ptId2];
    out.push_back(e);
  }
  return out;
}




void NormalizeEdgeWeights(Array<Edge> &edges, int numPts, double weightScaleFactor = 1.0 ){
  Array<double> ptWeightTotal(numPts, 0.0);

  for (auto e:edges){
    if (e.ptId1 < numPts)
      ptWeightTotal[e.ptId1] += e.weight1;
    if (e.ptId2 < numPts)
      ptWeightTotal[e.ptId2] += e.weight2;
  }

  for (auto &e:edges){
    double total1 = ptWeightTotal[e.ptId1];
    if (total1 > 0){
      e.weight1 /= total1;
    }

    double total2 = ptWeightTotal[e.ptId2];
    if (total2 > 0){
      e.weight2 /= total2;
    }
  }

}


void MapEdgePtIds(Array<Edge> &edges, Array<int> pt1Ids,
    Array<int> pt2Ids, int offset1 = 0,  int offset2 = 0){

  OrderedMap<int,int> pt1Mapping;
  for (int i=0; i<pt1Ids.size(); i++){
    int ptId = pt1Ids[i];
    pt1Mapping[ptId] = i + offset1;
  }

  OrderedMap<int,int> pt2Mapping;
  for (int i=0; i<pt2Ids.size(); i++){
    int ptId = pt2Ids[i];
    pt2Mapping[ptId] = i + offset2;
  }

  for(auto &edge: edges){
   auto it = pt1Mapping.find(edge.ptId1);
   bool found = it != pt1Mapping.end();
   if (found) { edge.ptId1 = it->second; }

   auto it2 = pt2Mapping.find(edge.ptId2);
   bool found2 = it2 != pt2Mapping.end();
   if (found2) { edge.ptId2 = it2->second; }
  }
}




/*
bool firstLessThan(const Pair<int, double>& l, const Pair<int, double>& r) {
  return l.first < r.first;
}

bool firstEquivalent(const Pair<int, double>& l, const Pair<int, double>& r) {
  return l.first == r.first;
}
*/


vtkSmartPointer<vtkPoints> InitializePtsFromSeeds(vtkPolyData * seedPts, Array<Edge> edges, int numPts){

  vtkSmartPointer<vtkPoints> pts =
      vtkSmartPointer<vtkPoints>::New();


  //extract adjs and edge lengths
  Array<  Array<int>  > adjsArray(numPts, Array<int>(0));
  Array<  Array<double>  > edgeLengthsArray(numPts, Array<double>(0));

  for (auto e: edges){
    if (e.weight1 > 0){
      adjsArray[e.ptId1].push_back(e.ptId2);
      edgeLengthsArray[e.ptId1].push_back(e.length);
    }
    if (e.weight2 > 0){
      adjsArray[e.ptId2].push_back(e.ptId1);
      edgeLengthsArray[e.ptId2].push_back(e.length);
    }

  }

  //compute weighted average of seed pts
  for (int otherPtId=0; otherPtId<numPts; otherPtId++){
    Array<int> adjs = adjsArray[otherPtId];
    Array<double> edgeLengths = edgeLengthsArray[otherPtId];


    if (adjs.size() == 0){
      cout << "MeshToSphere::InitializePtsFromSeeds" << endl;
      cout << "point " << otherPtId << " has no adjacencies"  << endl;
      exit(1);
    }

    double meanEdgeLength = 0.0;
    for (auto l: edgeLengths){ meanEdgeLength += l; }
    meanEdgeLength /= adjs.size();

    Array<double> weights;
    for (auto &l: edgeLengths){
      double w = 1 / (l + meanEdgeLength*0.01);
      weights.push_back(w);
    }

    //normalize weights
    double weightSum = 0.0;
    for (auto w: weights){ weightSum += w; }
    for (auto &w: weights){ w /= weightSum; }

    //reconstruct point
    double pt[3] = {0.0, 0.0, 0.0};

    for (int i=0; i<adjs.size(); i++){
      int seedPtId =  adjs[i];
      double w = weights[i];
      double seedPt[3];
      seedPts->GetPoint(seedPtId, seedPt);


      pt[0] += seedPt[0] * w;
      pt[1] += seedPt[1] * w;
      pt[2] += seedPt[2] * w;
    }
    double norm = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);
    pt[0] /= norm;     pt[1] /= norm;     pt[2] /= norm;

    pts->InsertNextPoint(pt);
  }
  return pts;
}


void MeshToSphere::Run(){

  int numLevels = _Parameters.size();

  cout << "running" << endl;

  ComputeMeshHierarchy();


  for (_CurrentLevel=0; _CurrentLevel < numLevels; _CurrentLevel++){

    vtkSmartPointer<vtkPolyData> mesh = _MeshHierarchy[_CurrentLevel];
    int numPts = mesh->GetNumberOfPoints();

    LevelParameters p = _Parameters[_CurrentLevel];
    Array<int> ptIds = _PtIdHierarchy[_CurrentLevel];

    /*
    M2SConnectivity conn;
    conn.Input(mesh);
    conn.Threshold(p.smds.distThreshold.value);
    //conn.SeedPtIds(seedPtIds);
    conn.IdentityEdgesOff();
    conn.ForwardDirectionOn();
    conn.BackwardDirectionOn();
    conn.NormalizeOn();
    conn.Update();

    Array<Edge> edges = conn.Edges();
    */

    if (_CurrentLevel == 0){
      InitialiseSphericalMesh();
    }
    else{
      //Interpolate angles

      //add spherical coordinates to output as pointdata

      /*
      vtkSmartPointer<vtkPolyData> previousMesh =
          vtkSmartPointer<vtkPolyData>::New();
      previousMesh->DeepCopy(_MeshHierarchy[_CurrentLevel-1]);

      string arrayName = "sphericalData";
      vtkDoubleArray * sphericalData = (vtkDoubleArray *) _Output->GetPoints()->GetData();
      sphericalData->SetName(arrayName.c_str());
      previousMesh->GetPointData()->AddArray(sphericalData);


      M2SGeodesicInterpolator gi;
      gi.Target(mesh);
      gi.Source(previousMesh);
      gi.K(p.interp.k.value);
      gi.Lambda(p.interp.lambda.value);
      gi.Iterations(p.interp.iterations.value);
      gi.ArrayName(arrayName);
      gi.OutArrayName(arrayName);
      gi.Norm(2); //normalise spherical coordinates with l2 norm
      gi.Run();

      sphericalData = (vtkDoubleArray *) gi.Output()->GetPointData()->GetArray(arrayName.c_str());

      vtkSmartPointer<vtkPoints> pts =
          vtkSmartPointer<vtkPoints>::New();
      pts->SetData(sphericalData);
      _Output->DeepCopy(mesh);
      _Output->SetPoints(pts);
      */

      InterpolateAngles(_CurrentLevel-1, _CurrentLevel);

      if (p.postprocess.active){
        M2SPostProcessor pp;
        pp.MaxDiffusionIterations(p.postprocess.diffusionIterations.value);
        pp.MaxIterations(p.postprocess.iterations.value);
        pp.Percentile(p.postprocess.percentile.value);
        pp.Threshold(p.postprocess.threshold.value);
        pp.Input(_Output);
        pp.Run();
        _Output = pp.Output();
        EmbedPointsOnSphere(_Output);


        if (_Debug){
          string outName = "level" + std::to_string(_CurrentLevel-1) + "_post_processed.vtk";
          WritePolyData(outName.c_str(), _Output);
        }

      }




    }

    if (! p.smds.active || p.smds.maxIterations.value == 0){ continue; }

  if (_Debug){
    string outName = "level" + std::to_string(_CurrentLevel) + "_conn_mesh.vtk";
      WritePolyData(outName.c_str(), mesh);
  }

    M2SConnectivity conn;
    conn.Input(mesh);
    conn.Threshold(p.smds.distThreshold.value);
    //conn.SeedPtIds(seedPtIds);
    conn.IdentityEdgesOff();
    conn.ForwardDirectionOn();
    conn.BackwardDirectionOn();
    conn.NormalizeOn();
    conn.Update();

    Array<Edge> edges = conn.Edges();

  if (_Debug){
    string outName = "level" + std::to_string(_CurrentLevel) + "_edges.txt";
    ofstream out(outName);
    for(int ei=0;ei<edges.size();ei++){
    Edge e=edges[ei];
          out<<"edge "<<e.ptId1<<" "<<e.ptId2<<" "<<e.length<<" "<<e.weight1<<" "<<e.weight2<<endl;
    }
    out.close();
  }

    M2SOptimizer optimizer;
    //optimizer.Edges(conn.Edges());
    optimizer.Edges(edges);
    optimizer.Input(_Output);
    optimizer.MaxIterations(p.smds.maxIterations.value);
    //optimizer.Debug(false);
    optimizer.Run();
    _Output = optimizer.Output();
  }


  /*
  if (_FinalInterpolation){
      //add spherical coordinates to output as pointdata

      vtkSmartPointer<vtkPolyData> previousMesh =
          vtkSmartPointer<vtkPolyData>::New();
      previousMesh->DeepCopy(_MeshHierarchy.back());

      string arrayName = "sphericalData";
      vtkDoubleArray * sphericalData = (vtkDoubleArray *) _Output->GetPoints()->GetData();
      sphericalData->SetName(arrayName.c_str());
      previousMesh->GetPointData()->AddArray(sphericalData);


      M2SGeodesicInterpolator gi;
      gi.Target(_Input);
      gi.Source(previousMesh);
      gi.K(100);
      gi.Lambda(0.01);
      gi.ArrayName(arrayName);
      gi.OutArrayName(arrayName);
      gi.Norm(2); //normalise spherical coordinates with l2 norm
      gi.Run();



      sphericalData = (vtkDoubleArray *) gi.Output()->GetPointData()->GetArray(arrayName.c_str());

      vtkSmartPointer<vtkPoints> pts =
          vtkSmartPointer<vtkPoints>::New();
      pts->SetData(sphericalData);
      _Output->DeepCopy(_Input);
      _Output->SetPoints(pts);
  }
  */




}





} // mirtk namespace



