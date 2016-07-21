

#include "mirtk/Common.h"
#include "mirtk/M2SOptimizer.h"
//#include <mirtkPointSetUtils.h>
#include "mirtk/Vtk.h"

#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================


M2SOptimizer::M2SOptimizer()
{
  _Energy = 0.0;
  _MaxIterations = 100;
  _LowPassIterations = 0;
  _LowPassBandwidth = 0.5;
  _Debug = false;
}


// -----------------------------------------------------------------------------


M2SOptimizer::M2SOptimizer(const M2SOptimizer &other)
{
  *this = other;
}

// -----------------------------------------------------------------------------

M2SOptimizer &M2SOptimizer::operator =(const M2SOptimizer &other)
{
  if (this != &other) {
    M2SOptimizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

M2SOptimizer::~M2SOptimizer(){}

// -----------------------------------------------------------------------------

void M2SOptimizer::CopyAttributes(const M2SOptimizer &other){
  cout << "TODO copy attributes" << endl;
}



void M2SOptimizer::NormalizeWeights(){
  Array<double> ptWeightTotal(_NumPts, 0.0);

  for (auto e:_Edges){
    ptWeightTotal[e.ptId1] += e.weight1;
    ptWeightTotal[e.ptId2] += e.weight2;
  }

  for (auto &e:_Edges){
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

void M2SOptimizer::ComputeTotalWeight(){
  _TotalWeight = 0.0;
  for (auto e:_Edges){
    _TotalWeight += e.weight1;
    _TotalWeight += e.weight2;
  }
}

void M2SOptimizer::Initialize(){

  _NumPts = _Input->GetNumberOfPoints();
  _NumEdges = _Edges.size();

  /*
  if (_StationaryPts.size() != _NumPts){
    _StationaryPts.clear();
    for (int ptId=0; ptId<_NumPts; ptId++){
      _StationaryPts.push_back(true);
    }
  }
  */




  //clear all arrays and reinitialise

  _CurrentSquareDistances.clear();
  _Errors.clear();
  _OriginalSphericalCoordinates.clear();
  _OldSphericalCoordinates.clear();
  _SphericalCoordinates.clear();
  _BestSphericalCoordinates.clear();
  _Gradient.clear();

  SphericalVector v;
  for (int ptId=0; ptId<_NumPts; ptId++){
    _OldSphericalCoordinates.push_back(v);
    _SphericalCoordinates.push_back(v);
    _BestSphericalCoordinates.push_back(v);
    _Gradient.push_back(v);
    //_StepSizeWpush_back(v);
  }


  for (int edgeId=0; edgeId<_NumEdges; edgeId++){
    _CurrentSquareDistances.push_back(0.0);
    _Errors.push_back(0.0);
  }



  _Output = vtkSmartPointer<vtkPolyData>::New();
  _Output->DeepCopy(_Input);
  UpdateSphericalCoordinatesFromMesh();
  NormalizeWeights();
  ComputeTotalWeight();
}

void M2SOptimizer::UpdateSphericalCoordinatesFromMesh(){

  ofstream of;
  if (_Debug){
    of.open ("/home/rob/playground/converting-spherical-coords-c++.txt", std::ofstream::out | std::ofstream::app);
  }

  int numPoints = _Output->GetNumberOfPoints();
  double pt [3];

  for (int ptId=0; ptId<numPoints; ptId++){
    _Output->GetPoint(ptId,pt);
    double norm = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);
    SphericalVector spt;
    spt.phi = atan2(pt[1], pt[0]); //azimuth
    spt.theta = acos(pt[2] / norm); //inclination
    _SphericalCoordinates[ptId] = spt;

    if (spt.phi != spt.phi || spt.theta != spt.theta ){
      cout << pt[0] << " " << pt[1] << " " << pt[2] << endl;
      cout << ptId << " " << spt.phi <<"," << spt.theta << endl;
      exit(1);
    }

    if (_Debug){
      of << "pt " << ptId << " xyz<"  << pt[0] << "," << pt[1] << "," << pt[2] << "> " << " phi/theta<"  << spt.phi << "," << spt.theta << ">" << "\n";
    }
  }

  if (_Debug){
    of << "\n";
    of.close();
  }

}


void M2SOptimizer::UpdateMeshFromSphericalCoordinates(){

  int numPoints = _Output->GetNumberOfPoints();

  double p[3];

  double max_phi = -1000;
  double min_phi = 1000;
  double max_theta = -1000;
  double min_theta = 1000;


  vtkSmartPointer <vtkDoubleArray> array;

  if (_Output->GetPointData()->HasArray("coords")){
    array = (vtkDoubleArray*) _Output->GetPointData()->GetArray("coords");
  }
  else{
    array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetName("coords");
    array->SetNumberOfComponents(2);
    array->SetNumberOfTuples(numPoints);
    _Output->GetPointData()->AddArray(array);
  }

  for (int i=0; i<numPoints; i++){
    double phi =  _SphericalCoordinates[i].phi;
    double theta =  _SphericalCoordinates[i].theta;

    //cout << "<" << theta << "," << phi << endl;
    /*
    if (theta > max_theta) max_theta = theta;
    if (phi > max_phi) max_phi = phi;
    if (theta < min_theta) min_theta = theta;
    if (phi < min_phi) min_phi = phi;
    */

    p[0] = sin(theta) * cos(phi);
    p[1] = sin(theta) * sin(phi);
    p[2] = cos(theta);
    _Output->GetPoints()->SetPoint(i,p);

    array->SetTuple2(i, phi, theta);

  }

  /*
  cout << "max theta " << max_theta << endl;
  cout << "min theta " << min_theta << endl;
  cout << "max phi " << max_phi << endl;
  cout << "min phi " << min_phi << endl;
  */

}


void M2SOptimizer::EmbedPointsOnSphere(){

  int numPoints = _Output->GetNumberOfPoints();
  double p [3];

  for (int ptId = 0; ptId < numPoints; ++ptId) {
    _Output->GetPoint(ptId, p);

    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

    for (int i = 0; i < 3; ++i){
      p[i] /= norm;
    }
    _Output->GetPoints()->SetPoint(ptId, p);
  }
}
// -----------------------------------------------------------------------------
/// Compute centroid of point set
void M2SOptimizer::GetCentroid(vtkPolyData *mesh, double centroid[3])
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
void M2SOptimizer::GetScale(vtkPolyData *mesh, const double centroid[3], double scale[3])
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




void M2SOptimizer::SmoothDisoriented(){
  UpdateMeshFromSphericalCoordinates();


  vtkCellArray * cells = _Output->GetPolys();
  cells->InitTraversal();


  Array<bool> disoriented;

  for (int ptId=0; ptId<_NumPts; ptId++){
    disoriented.push_back(false);
  }

  vtkSmartPointer<vtkIdList> cellIds =
      vtkSmartPointer<vtkIdList>::New();

  double p0[3], p1[3], p2[3], a[3], b[3], cross[3], centriod[3];

  while (cells->GetNextCell(cellIds)){
    int ptId0 = cellIds->GetId(0);
    int ptId1 = cellIds->GetId(1);
    int ptId2 = cellIds->GetId(2);
    _Output->GetPoint(ptId0, p0);
    _Output->GetPoint(ptId1, p1);
    _Output->GetPoint(ptId2, p2);

    for (int i=0; i<3; i++){
      a[i] = p2[i]-p0[i];
      b[i] = p1[i]-p0[i];
      centriod[i] = (p0[i] + p1[i] + p2[i]);
    }

    cross[0] = (a[1] * b[2]) - (a[2] * b[1]);
    cross[1] = (a[2] * b[0]) - (a[0] * b[2]);
    cross[2] = (a[0] * b[1]) - (a[1] * b[0]);


    for (int i=0; i<3; i++){
      a[i] = p2[i]-p0[i];
      b[i] = p1[i]-p0[i];
      centriod[i] = (p0[i] + p1[i] + p2[i])/3;
    }


    //don't need to normalize vector lengths as we're only interested in the sign
    double dot = centriod[0] * cross[0];
    dot += centriod[1] * cross[1];
    dot += centriod[2] * cross[2];

    if (dot < 0){
      disoriented[ptId0];
      disoriented[ptId1];
      disoriented[ptId2];
    }
  }
  AdaptiveSmooth(disoriented);
  EmbedPointsOnSphere();
  UpdateSphericalCoordinatesFromMesh();


}

void M2SOptimizer::AdaptiveSmooth(Array<bool> disoriented){

  int numIterations = 10;
  Array< vtkSmartPointer<vtkIdList> > adjs;

  for (int ptId=0; ptId<_NumPts; ptId++){
    adjs.push_back( vtkSmartPointer<vtkIdList>::New() );
  }
  vtkCellArray * cells = _Output->GetPolys();
  cells->InitTraversal();

  vtkSmartPointer<vtkIdList> cellIds =
      vtkSmartPointer<vtkIdList>::New();

  while (cells->GetNextCell(cellIds)){
    int ptId0 = cellIds->GetId(0);
    int ptId1 = cellIds->GetId(1);
    int ptId2 = cellIds->GetId(2);

    adjs[ptId0]->InsertUniqueId(ptId1);
    adjs[ptId0]->InsertUniqueId(ptId2);

    adjs[ptId1]->InsertUniqueId(ptId0);
    adjs[ptId1]->InsertUniqueId(ptId2);

    adjs[ptId2]->InsertUniqueId(ptId0);
    adjs[ptId2]->InsertUniqueId(ptId1);
  }


  vtkPoints * pts = _Output->GetPoints();


  vtkSmartPointer<vtkPoints> backup =
      vtkSmartPointer<vtkPoints>::New();

  double centriod [3], aPt[3], pt[3];
  double alpha = 0.5, beta = 1-alpha;

  for (int iteration=0; iteration < numIterations; iteration++){


    backup->DeepCopy(pts);
    cout << "ad smoothing" << endl;
    for (int ptId=0; ptId<_NumPts; ptId++){
      if (disoriented[ptId]){
        cout << "ad smoothing" << endl;
        vtkIdList * aPtIds = adjs[ptId];
        int k = aPtIds->GetNumberOfIds();

        centriod[0] = 0; centriod[1] = 0; centriod[2] = 0;
        for (int i=0; i<k;i++){
          int aPtId = aPtIds->GetId(i);
          backup->GetPoint(aPtId, aPt);
          centriod[0] += aPt[0];
          centriod[1] += aPt[1];
          centriod[2] += aPt[2];
        }

        centriod[0] /= k; centriod[1] /= k; centriod[2] /= k;

        backup->GetPoint(ptId, pt);

        for (int i=0; i<3; i++){
          pt[i] = alpha * centriod[i] + beta * pt[i];
        }

        pts->SetPoint(ptId, pt);
      }
    }

  }

}


void M2SOptimizer::LowPassFilter(){

  double c1[3], s1[3];
  GetCentroid(_Output, c1);
  GetScale(_Output, c1, s1);

  // Perform low-pass filtering
     //
     // The translation and scale "fix" is due to a bug in vtkWindowedSincPolyDataFilter:
     // http://vtk.1045678.n5.nabble.com/Bug-in-vtkWindowedSincPolyDataFilter-td1234055.html

  vtkNew<vtkWindowedSincPolyDataFilter> filter;
  filter->SetPassBand(_LowPassBandwidth);
  filter->SetNumberOfIterations(_LowPassIterations);
  filter->NormalizeCoordinatesOn();
  SetVTKInput(filter, _Output);
  filter->Update();

  double c2[3], s2[3];
  GetCentroid(filter->GetOutput(), c2);
  GetScale(filter->GetOutput(), c2, s2);

  // Adjust point displacements
  double p1[3], p2[3], sf[3];

  for (int i=0; i<0; i++) { sf[i] = s2[i] / s1[i]; } //approximate scale factor

  for (vtkIdType ptId = 0; ptId < _NumPts; ++ptId) {
    _Output->GetPoint(ptId, p1);
    filter->GetOutput()->GetPoint(ptId, p2);

    for (int i=0; i<0; i++){
      p2[i] = ((p2[i] - c2[i]) / sf[i]) + c1[i];
    }

    _Output->GetPoints()->SetPoint(ptId, p2);
  }
}



void rescaleTargetDistances(Array<double> &currentSquareDistances, Array<double> &targetDistances){


  double totalCurrentDist = 0.0;
  double totalTargetDist = 0.0;

  for (int i=0; i<targetDistances.size(); i++){
    double g = acos((2 - currentSquareDistances[i]) / 2);
    double d = targetDistances[i];

    totalCurrentDist += g;
    totalTargetDist += d;

    //cout << "g " << totalCurrentDist << " d " << totalTargetDist << endl;



  }

  double scale = totalCurrentDist / totalTargetDist;

  cout << "scaling target distances by " << scale << endl;


  for (int i=0; i<targetDistances.size(); i++){
    targetDistances[i] *= scale;
  }

}


void M2SOptimizer::Optimize(){


  double maxStepSize = 1.0;
  int maxEvaluations = 50;
  int minEvaluations = 5;
  double stepSize = 1.0;
  double epsilon =   0.000005;
  double bestStepSize = -1.0;

  for (int iteration=0; iteration<_MaxIterations; iteration++){
    ComputeCurrentSquareDistances();
    ComputeErrors();
    ComputeEnergy();
    ComputeGradient();
    //SmoothGradient();
    //cout << "smooth gradients" << endl;

    double bestEnergy = _Energy; //TO_REMOVE
    double oldBestEnergy = _Energy;
    _OldSphericalCoordinates = _SphericalCoordinates;


    cout << "iteration:  " << iteration << " energy:  "
        << _Energy << " step_size:   " << bestStepSize  <<  endl;

    int e;
    for (e=0; e<maxEvaluations; e++){
      //compute new angles
      UpdateSphericalCoordinates(stepSize);
      ComputeCurrentSquareDistances();
      ComputeErrors();
      ComputeEnergy();

      if (_Energy < bestEnergy){
        bestStepSize = stepSize;
        _BestSphericalCoordinates = _SphericalCoordinates;
        bestEnergy = _Energy;
      }
      stepSize /= 2.0;

      if (e > minEvaluations && bestStepSize != -1.0){
        break;
      }
    }

    if (bestStepSize == -1.0){
      cout << "no best step size found" << endl;
      cout << e << endl;
      break;

      //rescaleTargetDistances(_CurrentSquareDistances, _TargetDistances);
    }

    _SphericalCoordinates = _BestSphericalCoordinates;



    /*
    double improvement = (oldBestEnergy - bestEnergy) / oldBestEnergy;
    */


    stepSize = bestStepSize * 2.0;

    if (stepSize > maxStepSize)
      stepSize = maxStepSize;


    oldBestEnergy = bestEnergy;

    /*
    if (improvement < epsilon){
        cout << "Finished - improvement not sufficient" << endl;
        cout << improvement << endl;
        break;
    }
    */


    if (_LowPassIterations > 0){
      cout << "lowpass" << endl;
      UpdateMeshFromSphericalCoordinates();
      LowPassFilter();
      //SmoothDisoriented();
      //Smooth(_Output);
      UpdateSphericalCoordinatesFromMesh();
    }

  }
  UpdateMeshFromSphericalCoordinates();
}




double squareSphericalDistance(double i1, double  i2,double  j1,double  j2){

    double squareDist = 2 * (1 -
      sin(i2) * sin(j2) * cos(i1-j1) -
      cos(i2) * cos(j2));

    return squareDist;
}


void M2SOptimizer::ComputeCurrentSquareDistances(){


  double totalDistance = 0;
  double totalTargetDistance = 0;

  ofstream of;
  if (_Debug){
    of.open ("/home/rob/playground/sqr-distances-c++.txt", std::ofstream::out | std::ofstream::app);
  }

  for(int edgeId=0; edgeId<_Edges.size(); edgeId++){
    Edge e = _Edges[edgeId];

    SphericalVector spt1 = _SphericalCoordinates[e.ptId1]; //spt = spherical point
    SphericalVector spt2 = _SphericalCoordinates[e.ptId2];

    double d = squareSphericalDistance(spt1.phi,spt1.theta, spt2.phi,spt2.theta);

    if (d<=0){ d=numeric_limits<double>::epsilon(); }
    _CurrentSquareDistances[edgeId] = d;

    if (d!=d){
      cout << e.ptId1 << " " << spt1.phi <<"," << spt1.theta << endl;
      cout << e.ptId2 << " " << spt2.phi <<"," << spt2.theta << endl;
      exit(1);
    }


    if (_Debug){
      of << "edge " << edgeId << " <"  << e.ptId1 << "," << e.ptId2 << "> " << _CurrentSquareDistances[edgeId] << "\n";
    }

  }

  if (_Debug){
    of << "\n";
    of.close();
  }


}


void M2SOptimizer::ComputeErrors(){

  ofstream of;
  if (_Debug){
    of.open ("/home/rob/playground/errors-c++.txt", std::ofstream::out | std::ofstream::app);
  }


  double totalError = 0.0;
  double maxError = 0;
  int maxId = 0;
  int numEdges = _Edges.size();
  for (int edgeId=0; edgeId<numEdges; edgeId++){
    Edge e = _Edges[edgeId];
    double d2 = _CurrentSquareDistances[edgeId];

    if (d2 <= 0){
      _Errors[edgeId] = 0;
    }
    else{
      _Errors[edgeId] = acos((2 - d2) / 2) - e.length;
    }


    if (_Errors[edgeId]  != _Errors[edgeId]){
      cout << "edge " << edgeId << " error " << _CurrentSquareDistances[edgeId] << endl;
      cout << e.ptId1 << "," << e.ptId2 << " " << e.length << " " << e.weight1 << " " << e.weight2 <<  endl;
      exit(1);
    }


    if (_Debug){
      of << "edge " << edgeId << " " << _Errors[edgeId] << "\n";
    }

    //totalError += abs(error);



    /*
    if (abs(error) > maxError){
      maxId = i;
      maxError = abs(error);
    }*/

  }

  if (_Debug){
    of << "\n";
    of.close();
  }


/*
  double averageError = totalError / numEdges ;

  cout << "totalError " << totalError << endl;
  cout << "averageError " << averageError << endl;
  cout << "maxError " << maxError << endl;
  cout << "maxId " << maxId << endl;
  cout << _CurrentSquareDistances[maxId] << "," << _TargetDistances[maxId] << endl;
*/
  /*
  EdgeIterator ei(_Edges);
  ei.InitTraversal();

  vtkIdType pid1, pid2;
  int i=0;

  double totalDistance = 0;
  int edgeid = 0;
  while ( ei.GetNextEdge(pid1, pid2) != -1 ){
    if (edgeid == maxId){
      cout << pid1 << "<" << pid2 << endl;

      cout << "<"<< _SphericalCoordinates[pid1].theta << "," <<
          _SphericalCoordinates[pid1].phi << ">" << endl;

      cout << "<"<< _SphericalCoordinates[pid2].theta << "," <<
          _SphericalCoordinates[pid2].phi << ">" << endl;

      double d = squareSphericalDistance(_SphericalCoordinates[pid1].theta,
                              _SphericalCoordinates[pid1].phi,
                              _SphericalCoordinates[pid2].theta,
                              _SphericalCoordinates[pid2].phi);

      cout <<"spohd " << d << endl;
      break;
    }
    edgeid++;
  }
  */

  //TODO implement as a optimization parameter
  /*
  bool centerErrors = true;
  if (centerErrors){
    totalError/=numEdges;
    for (int i=0; i<numEdges; i++){
      _Errors[i]-=totalError;
    }
  }
  */


}


void M2SOptimizer::ComputeEnergy(){

  int numEdges = _Edges.size();
  _Energy = 0.0;
  for (int edgeId=0; edgeId<numEdges; edgeId++){
    Edge e = _Edges[edgeId];
    _Energy += (e.weight1 + e.weight2) * _Errors[edgeId] * _Errors[edgeId];

    if (_Errors[edgeId]  != _Errors[edgeId]){
      cout << "edge " << edgeId << " error " << _Errors[edgeId] << "d2 " << _CurrentSquareDistances[edgeId]  << endl;
      cout << e.ptId1 << "," << e.ptId2 << " " << e.length << " " << e.weight1 << " " << e.weight2 <<  endl;
      exit(1);

    }

  }

  _Energy /= _TotalWeight;

  //cout << "energy " << _Energy << endl;
}





void M2SOptimizer::ComputeGradient(){

  for (int ptId=0; ptId < _NumPts; ptId++){
    _Gradient[ptId].theta = .0;
    _Gradient[ptId].phi = .0;
  }



  ofstream of;
  if (_Debug){
    of.open ("/home/rob/playground/gradients-c++.txt", std::ofstream::out | std::ofstream::app);
  }


  for (int edgeId=0; edgeId<_Edges.size(); edgeId++){
    Edge e = _Edges[edgeId];

    double i1 = _SphericalCoordinates[e.ptId1].phi;
    double i2 = _SphericalCoordinates[e.ptId1].theta;

    double j1 = _SphericalCoordinates[e.ptId2].phi;
    double j2 = _SphericalCoordinates[e.ptId2].theta;

    double d =  _CurrentSquareDistances[edgeId];
    double a = 2 * _Errors[edgeId];

    double b = -2 / (sqrt(1 - pow(((2 - d) / 2),2)) + numeric_limits<double>::epsilon());


    if (a!=a){ cout << edgeId << "a" << a << endl; exit(1); }
    if (b!=b){ cout << edgeId << "b" << b << endl; cout << edgeId << "d" << d << endl; exit(1); }
    if (d!=d){ cout << edgeId << "d" << d << endl; exit(1); }


    if (e.weight1 > 0){
      double c1 = -sin(i2) * sin(j2) * sin(i1 - j1);
      double c2 = -cos(i2)*sin(j2)*cos(i1-j1) + sin(i2)*cos(j2);
      _Gradient[e.ptId1].phi += a*b*c1 * e.weight1;
      _Gradient[e.ptId1].theta += a*b*-c2 * e.weight1;
    }


    if (e.weight2 > 0){
      double c1 = -sin(j2) * sin(i2) * sin(j1 - i1);
      double c2 = -cos(j2)*sin(i2)*cos(j1-i1) + sin(j2)*cos(i2);

      _Gradient[e.ptId2].phi += a*b*c1 * e.weight2;
      _Gradient[e.ptId2].theta += a*b*-c2 * e.weight2;
    }
  }


  /*

  */

  //normalize gradients
  if (_Debug){
    for (int ptId=0; ptId < _NumPts; ptId++){
      of << "pt " << ptId << " <"  << _Gradient[ptId].phi << "," << _Gradient[ptId].theta << "> " << "\n";
    }
    of << "\n";
    of.close();
  }

}









void M2SOptimizer::UpdateSphericalCoordinates(double stepSize){
  int numCoordinates = _SphericalCoordinates.size();


  for (int ptId = 0; ptId<numCoordinates; ptId++){
    //if (_ActivePoints[ptId]){

      double phi = _OldSphericalCoordinates[ptId].phi;
      double theta = _OldSphericalCoordinates[ptId].theta;

      phi -= stepSize * _Gradient[ptId].phi;
      theta -= stepSize * _Gradient[ptId].theta;

      if (phi!=phi || theta!=theta){
        cout << ptId << " " << phi <<"," << theta << endl;
        cout << _Gradient[ptId].phi <<"," << _Gradient[ptId].theta << endl;
        cout << _OldSphericalCoordinates[ptId].phi <<"," << _OldSphericalCoordinates[ptId].theta << endl;
        exit(1);
      }


      //_OldSphericalCoordinates[ptId].phi - stepSize * _StepSizeW[ptId].phi * _Gradient[ptId].phi;


      //don't clamp inclination! messes things up

      /*
      if (theta > pi){
        theta = pi;
      }
      if (theta < 0){
        theta = 0;
      }
      */


      /*

      //wrap phi
      if (phi > pi)
        phi -= two_pi;
      if (phi < -pi)
        phi += two_pi;

      */

      /*
      //clamp inclination
      if (theta > pi || theta < 0){
        cout << "theta out of range for " << ptId << " " << theta << endl;
        cout <<"stepSize"  << stepSize << endl;
        cout <<"_Gradient[ptId].phi "  << _Gradient[ptId].phi << endl;
        exit(1);
      }
      */

      /*
      //no need to wrap phi!
      /*
      if (phi > pi || phi < -pi){
        cout << "phi out of range "<< phi << endl;
        exit(1);
      }
       */





      _SphericalCoordinates[ptId].phi = phi;
      _SphericalCoordinates[ptId].theta = theta;


  }
}


void M2SOptimizer::Run(){
  Initialize();


  if (_Edges.size() < 1){
    cout << "no edges!" << endl;
    exit(1);
  }

  cout << "Optimizing " << _Input->GetNumberOfPoints() << " points and "
      << _Edges.size()<< " edges." << endl;


  //Smooth(_Output);
  //EmbedPointsOnSphere();
  //UpdateSphericalCoordinatesFromMesh();

  Optimize();

}

} // mirtk namespace

