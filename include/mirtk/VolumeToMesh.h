/*
 * volumeToMesh.h
 *
 *  Created on: 15 Oct 2013
 *      Author: rob
 */

#ifndef VOLUMETOMESH_H_
#define VOLUMETOMESH_H_


#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/EuclideanDistanceTransform.h"
#include "mirtk/Remesher.h"
 //MIRTK Remeshing
// #include "mirtk/SurfaceRemeshing.h"
#include "mirtk/LaplacianSmoothFilter.h"

#include "mirtk/ImageFunction.h"

#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkModifiedBSPTree.h>
#include <vtkOBBTree.h>



namespace mirtk {

class VolumeToMesh{

  enum InitialMesh {input, hull, sphere};

private:
  bool _multiResolutionOn,_screenShotsOn,_selfIntersectionOn,_intersectionOn,
  _meshInitialized, _lowPassFilterOn,_finalRefinementOn,_smoothingOn,
  _flipNormals,_closeBinaryVolume,_gaussianInterpolationOn;
  char * _screenShotDir, *_smoothArrayName;

  double _smoothing, _averageEdgeLengthMM, _boundaryThreshold, _resolution,
  _maxStepSize, _maxStepSizeMM, _maxStepSizeVoxels, _selfIntersectionVoxels,
  _selfIntersectionMM, _boundingRadius, _improvementThreshold,
  _lowPassFilterBand, _refinementSmoothing, _stepSizeSmoothingSig,_stepSizeSmoothingSigMM;

  double _centreOfMass [3];
  int _maxIterations,_minIterations,_levels,_selfIntersectionTreeLeavesPerNode,_minNumPoints,_subdivisions,
  _lowPassFilterIterations,_refinementIterations, _stepSizeSmoothingIterations;
  int _boundingBoxP0 [3], _boundingBoxP1 [3];
  // vtkPolyData mesh;
  GreyImage _binaryVolume;
  RealImage _distanceField;
  vtkSmartPointer<vtkPolyData> _mesh;
  EuclideanDistanceTransform<RealPixel> *_edt;
  vtkSmartPointer<vtkPolyDataNormals> _normalFilter;
  Remesher* _remeshFilter;
 //MIRTK Remeshing
  // SurfaceRemeshing* _remeshFilter;
  LaplacianSmoothFilter* _smoothingFilter;
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> _lowPassFilter;
  ImageFunction *_interpolator;

  vtkSmartPointer<vtkModifiedBSPTree> _selfIntersectBSPTree;

  InitialMesh _initialMesh;

  void Initialize();
  void InitializeParameters();
  void InitializeMesh();
  void InitializeMeshAsBoundingSphere();
  void InitializeMeshAsConvexHull();
  void InitializeDistanceField();
  void InitializeFilters();

  bool Intersect(double [], double []);
  double GetVector(double [],double [],double []);
  void GetAdjacency(vtkSmartPointer<vtkCellArray>, vtkSmartPointer<vtkIdList> *);
  void ComputeMeshNormals();
  void SmoothMesh();
  void TransformPointsToWorldCoordinates();
  void TransformPointsToImageCoordinates();
  double UpdatePointPositions();
  void Subdivide();
  void Screenshot(string);
  void GetAdjs(vtkSmartPointer<vtkIdList> adj []);
  void InterpolateVelocities(double**, bool  [],  vtkSmartPointer<vtkIdList>  []);
  void AdaptivelySmoothMesh();
  void Remesh();
  void LowPassFilter();
  vtkSmartPointer<vtkPolyData> GetHull();
  void DeformMeshToTarget();
  bool SelfIntersect(double [3], double [3], double *);
  double Repulsion(double [3], double [3]);
  void SetInputAsTarget();
  void BuildSelfIntersectionTree();
  void CalculateCentreOfMass();
  bool InsideVolume(double [3]);
  void FinalRefinement();
  void PrintSettings();

  vtkSmartPointer<vtkPolyData> BinaryVolume2VTKSurfacePoints();
  vtkSmartPointer<vtkImageData> Mesh2VtkMask(vtkSmartPointer<vtkPolyData>);
public:

// Constructor
  VolumeToMesh();

/// Deconstuctor
  ~VolumeToMesh();




/// Set input image for filter
  void SetInput (RealImage *);
/// Set output image for filter
 void SetMaxIterations(int);
 void SetMinIterations(int);
 void SetEdgeLength(double);
 void SetNumberOfLevels(int);


 void SetScreenShotDirectory(char *);
 void SetInitialMeshToSphere();
 void SetInitialMeshToHull();
 void SetInitialMesh(vtkSmartPointer<vtkPolyData>);

 void SetEdgeSmoothing(double);
 void SetImprovementThreshold(double);

 void ScreenShotsOn();
 void ScreenShotsOff();
 void SelfIntersectionOn();
 void SelfIntersectionOff();
 void LowPassFilterOn();
 void LowPassFilterOff();
 void FinalRefinementOn();
 void FinalRefinementOff();
 void SmoothingOn();
 void SmoothingOff();
 void GaussianInterpolationOn();
 void GaussianInterpolationOff();
 void SetMaxStepSize(double);

 void SetStepSizeSmoothing(double, int);


 void SetLowPassFilterBand(double);
 void SetLowPassFilterIterations(int);


 void SetSmoothingValue(double);
 void SetRefinementIterations(int);
 void SetBoundaryThreshold(double);
 void CloseBinaryVolumeOn();
 void CloseBinaryVolumeOff();

 vtkSmartPointer<vtkPolyData>  GetOuput();


};

} //namespace mirtk

#endif /* VOLUMETOMESH_H_ */
