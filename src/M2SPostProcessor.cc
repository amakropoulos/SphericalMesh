

//mirtk includes
#include "mirtk/Common.h"
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/InflationForce.h"
#include "mirtk/InflationStoppingCriterion.h"
#include "mirtk/LocalOptimizer.h"
#include "mirtk/EulerMethodWithMomentum.h"
#include "mirtk/MetricDistortion.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/SparseMatrix.h"
#include "mirtk/MeshSmoothing.h"

//my mirtk includes
#include "mirtk/M2SPostProcessor.h"
#include "mirtk/M2SConnectivity.h"
#include "mirtk/M2SOptimizer.h"
#include "mirtk/M2SDiffuser.h"

//vtk includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkMassProperties.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkFastMarchingGeodesicDistance.h>
#include "vtkNew.h"
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCurvatures.h>


namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================

M2SPostProcessor::M2SPostProcessor(){
  _Debug = false;
  _Threshold = -1.0;
  _Percentile = 99;
  _MaxIterations = 1;
  _MaxDiffusionIterations = 1000;
  _Input = vtkSmartPointer<vtkPolyData>::New();
  _Output = vtkSmartPointer<vtkPolyData>::New();
}

// -----------------------------------------------------------------------------


M2SPostProcessor::M2SPostProcessor(const M2SPostProcessor &other)
{
  *this = other;
}


// -----------------------------------------------------------------------------

M2SPostProcessor &M2SPostProcessor::operator =(const M2SPostProcessor &other)
{
  if (this != &other) {
    M2SPostProcessor::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
M2SPostProcessor::~M2SPostProcessor(){
}

void M2SPostProcessor::CopyAttributes(const M2SPostProcessor &other){
  cout << "TODO copy attributes" << endl;
}


// =============================================================================
// Class methods
// =============================================================================

void projectToSphere(vtkPolyData * poly){
  int numPoints = poly->GetNumberOfPoints();

  double p[3];

  for (int i=0; i<numPoints; i++){
    poly->GetPoint(i,p);
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    p[0] /= norm; p[1] /= norm; p[2] /= norm;
    poly->GetPoints()->SetPoint(i,p);
  }
}

vtkSmartPointer<vtkPolyData> smoothScalars(vtkPolyData * mesh, string arrayName){
  cout << "smoothing scalars" << arrayName.c_str() << endl;
  if (! mesh->GetPointData()->HasArray(arrayName.c_str())){
    cout << "array " << arrayName.c_str() << "not found" << endl;
  }


 MeshSmoothing smoother;
 smoother.SetInputData(mesh);
 smoother.SmoothPointsOff();
 smoother.SetNumberOfIterations(2);
 //smoother.SignedSmoothing(true);
 //smoother.Weighting(PolyDataSmoothing::WeightFunction::Combinatorial);
 smoother.SmoothArray(arrayName.c_str());
 smoother.Update();
 smoother.Run();
 return smoother.GetOutput();
}

void M2SPostProcessor::DetectDiscontinuities(){


  cout << "Detecting Discontinuities " << endl;

  int numPts = _Output->GetNumberOfPoints();
  _DiscontinuityMask.clear();

  for (int ptId=0; ptId<numPts; ptId++){
    _DiscontinuityMask.push_back(false);
  }

  string k1Name = "Maximum_Curvature";
  string k2Name = "Minimum_Curvature";


  vtkNew<vtkCurvatures> k1Filter;
  k1Filter->SetInputData(_Output);
  k1Filter->SetCurvatureTypeToMaximum();
  k1Filter->Update();

  vtkDoubleArray * k1Array = (vtkDoubleArray *)
      k1Filter->GetOutput()->GetPointData()->GetArray(k1Name.c_str());

  _Output->GetPointData()->AddArray(k1Array);

  vtkNew<vtkCurvatures> k2Filter;
  k2Filter->SetInputData(_Output);
  k2Filter->SetCurvatureTypeToMinimum();
  k2Filter->Update();

  vtkDoubleArray * k2Array = (vtkDoubleArray *)
      k2Filter->GetOutput()->GetPointData()->GetArray(k2Name.c_str());

  _Output->GetPointData()->AddArray(k2Array);

  for (int ptId=0; ptId<numPts; ptId++){
    double k1 = k1Array->GetTuple1(ptId);
    double k2 = k2Array->GetTuple1(ptId);
    if (k1 < 0)
      _DiscontinuityMask[ptId] = true;
    if (k2 < 0)
      _DiscontinuityMask[ptId] = true;
  }

  vtkSmartPointer<vtkDoubleArray> maxAbsCurvature =
      vtkSmartPointer<vtkDoubleArray>::New();

  maxAbsCurvature->SetNumberOfComponents(1);
  string maxAbsCurvatureArrayName = "max_abs_curvature";

  maxAbsCurvature->SetName(maxAbsCurvatureArrayName.c_str());

  for (int ptId=0; ptId<numPts; ptId++){
    double k1 = k1Array->GetTuple1(ptId);
    double k2 = k2Array->GetTuple1(ptId);
    double maximum = max(abs(k1),abs(k2));
    maxAbsCurvature->InsertNextTuple1(maximum);
  }
  maxAbsCurvature->Modified();


  if (_Threshold == -1.0){

    Array<double> curvatures;
    for (int ptId=0; ptId<numPts; ptId++){
      double k = maxAbsCurvature->GetTuple1(ptId);
      curvatures.push_back(abs(k));
    }
    sort(curvatures.begin(),curvatures.end());

    int thresholdId = floor( _Percentile * numPts / 100);
    _Threshold = curvatures[thresholdId];
  }
  cout << "threshold " << _Threshold << endl;


  _Output->GetPointData()->AddArray(maxAbsCurvature);


  _Output = smoothScalars(_Output, maxAbsCurvatureArrayName);

  maxAbsCurvature = (vtkDoubleArray *)
    _Output->GetPointData()->GetArray(maxAbsCurvatureArrayName.c_str());




  for (int ptId=0; ptId<numPts; ptId++){
    double k = maxAbsCurvature->GetTuple1(ptId);
    if (k > _Threshold)
      _DiscontinuityMask[ptId] = true;
  }

  //add vtk array to mesh
  vtkSmartPointer <vtkDoubleArray> array =
      vtkSmartPointer <vtkDoubleArray>::New();
  for (int ptId=0; ptId<numPts; ptId++){
    array->InsertTuple1(ptId, _DiscontinuityMask[ptId]);
  }
  array->SetName("curvature_mask");
  _Output->GetPointData()->AddArray(array);
  _Output->GetPointData()->AddArray(maxAbsCurvature);

}


void M2SPostProcessor::SmoothDiscontinuities(){
  //compute edge lengths

  cout << "Smoothing discontinuities" << endl;

  int numPts = _Output->GetNumberOfPoints();

  Array<bool> sourceMask; //invert smooth mask to make source pt mask
  for (int ptId=0; ptId< numPts; ptId++){
    sourceMask.push_back(!_DiscontinuityMask[ptId]);
  }

  M2SDiffuser diffuser;
  diffuser.Input(_Output);
  diffuser.SourceMask(sourceMask);
  diffuser.MaxIterations(_MaxDiffusionIterations);
  diffuser.AdaptiveTermination(true);
  if(_Debug){
    diffuser.Debug(true);
  }
  diffuser.Run();

  _Output = diffuser.Output();

  projectToSphere(_Output);
}

int M2SPostProcessor::GetNumberOfDiscontinuities(){
  int numPts = _DiscontinuityMask.size();

  int numDiscontinuities =0;
  for (int ptId=0; ptId<numPts; ptId++){
    if(_DiscontinuityMask[ptId])
      numDiscontinuities++;
  }
  return numDiscontinuities;
}



void M2SPostProcessor::Run(){

  _Output->DeepCopy(_Input);

  cout << "num points " << _Output->GetNumberOfPoints() << endl;

  for (int iteration=0; iteration<_MaxIterations; iteration++){
    DetectDiscontinuities();

    int numDiscontinuities = GetNumberOfDiscontinuities();
    cout << numDiscontinuities << " discontinuities detected";

    if (numDiscontinuities == 0){
      break;
    }

    SmoothDiscontinuities();
  }
}



} // mirtk namespace



