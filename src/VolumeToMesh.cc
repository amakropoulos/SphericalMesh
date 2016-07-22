/*
 * VolumeToMesh.cc
 *
 *  Created on: 15 Oct 2013
 *      Author: rob
 */
#include "mirtk/VolumeToMesh.h"

#include "mirtk/ParticleSampleSphere.h"
#include "mirtk/Resampling.h"

#include "mirtk/Erosion.h"
#include "mirtk/Dilation.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/GaussianInterpolateImageFunction.h"

#include <math.h>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkArray.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPointData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkPointSet.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkButterflySubdivisionFilter.h>
#include <math.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkDoubleArray.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>
#include <vtkGraphicsFactory.h>
// #include <vtkImagingFactory.h>
#include <vtkOpenGLRenderer.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkHull.h>

#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkStructuredPoints.h>
#include <vtkMarchingCubes.h>
#include <vtkImageResample.h>
#include <vtkReverseSense.h>

//Constructor and Destructor
namespace mirtk {

VolumeToMesh::VolumeToMesh(){
  _selfIntersectBSPTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
  _selfIntersectBSPTree->SetMaxLevel(100);
  _selfIntersectBSPTree->SetNumberOfCellsPerNode(1);
  _selfIntersectBSPTree->AutomaticOn();
  _intersectionOn=true;

  _improvementThreshold=0.0001;
  _boundaryThreshold=0.5;
  _maxIterations=200;
  _subdivisions=1;
  _averageEdgeLengthMM=0.5;
  _maxStepSizeMM=1;
  _selfIntersectionMM=3;

  _stepSizeSmoothingSigMM = 1;
  _stepSizeSmoothingIterations = 0;

  _smoothArrayName = "_smoothingWeight";

  _levels=1;
  _smoothing=0.5;
  _minNumPoints=1000;
  _multiResolutionOn=false;
  _normalFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
  _remeshFilter = new Remesher();
  //MIRTK Remeshing
  // _remeshFilter = new SurfaceRemeshing();
  _smoothingFilter = new LaplacianSmoothFilter();
  _screenShotsOn=false;
  _screenShotDir="./";
  _selfIntersectionOn=true;
  _lowPassFilterOn=false;
  _finalRefinementOn=false;
  _lowPassFilterBand=0.5;
  _lowPassFilterIterations=100;
  _refinementSmoothing=0.5;
  _smoothingOn=true;
  _lowPassFilter = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  _flipNormals=false;
  _refinementIterations=0;
  _closeBinaryVolume=false;
  _gaussianInterpolationOn=false;
}





VolumeToMesh::~VolumeToMesh(){
  delete _remeshFilter;
  delete _smoothingFilter;
  delete _interpolator;
}



//Get and Set Methods
void VolumeToMesh::SetMaxIterations(int maxIterations){
  _maxIterations=maxIterations;
}

void VolumeToMesh::SetMinIterations(int minIterations){
  _minIterations=minIterations;
}


void VolumeToMesh::SetEdgeLength(double averageEdgeLength){
  _averageEdgeLengthMM=averageEdgeLength;
}

void VolumeToMesh::SetImprovementThreshold(double improvementThreshold){
  _improvementThreshold=improvementThreshold;
}

void VolumeToMesh::SetNumberOfLevels(int levels){
  _levels=levels;
}

void VolumeToMesh::ScreenShotsOn(){
  _screenShotsOn=true;
}


void VolumeToMesh::GaussianInterpolationOn(){
  _gaussianInterpolationOn = true;
}

void VolumeToMesh::GaussianInterpolationOff(){
  _gaussianInterpolationOn = false;
}

void VolumeToMesh::ScreenShotsOff(){
  _screenShotsOn=false;
}

void VolumeToMesh::SetScreenShotDirectory(char * screenShotDir){
  _screenShotDir=screenShotDir;
}

void VolumeToMesh::SelfIntersectionOn(){
  _selfIntersectionOn=true;
}
void VolumeToMesh::SelfIntersectionOff(){
  _selfIntersectionOn=false;
}

void VolumeToMesh::SetMaxStepSize(double maxStepSize){
  _maxStepSizeMM=maxStepSize;
}

void VolumeToMesh::SetStepSizeSmoothing(double stepSizeSmoothingSig,
    int stepSizeSmoothingIterations){
  _stepSizeSmoothingSigMM = stepSizeSmoothingSig;
  _stepSizeSmoothingIterations = stepSizeSmoothingIterations;
}


void VolumeToMesh::LowPassFilterOn(){
  _lowPassFilterOn=true;
}
void VolumeToMesh::LowPassFilterOff(){
  _lowPassFilterOn=false;
}

void VolumeToMesh::SetLowPassFilterBand(double lowPassFilterBand){
  _lowPassFilterBand=lowPassFilterBand;
}

void VolumeToMesh::SetLowPassFilterIterations(int lowPassFilterIterations){
  _lowPassFilterIterations=lowPassFilterIterations;
}

void VolumeToMesh::FinalRefinementOn(){
  _finalRefinementOn=true;
}

void VolumeToMesh::FinalRefinementOff(){
  _finalRefinementOn=false;
}

void VolumeToMesh::SetSmoothingValue(double smoothing){
  _smoothing=smoothing;
}


void VolumeToMesh::SmoothingOn(){
  _smoothingOn=true;
}

void VolumeToMesh::SmoothingOff(){
  _smoothingOn=false;
}

void VolumeToMesh::SetRefinementIterations(int iterations){
  _refinementIterations=iterations;
}

void VolumeToMesh::SetBoundaryThreshold(double threshold){
  _boundaryThreshold=threshold;
}

void VolumeToMesh::CloseBinaryVolumeOn(){
  _closeBinaryVolume=true;
}

void VolumeToMesh::CloseBinaryVolumeOff(){
  _closeBinaryVolume=false;
}

//Initialization Methods

void VolumeToMesh::Initialize(){
  InitializeDistanceField();
  InitializeFilters();
  InitializeMesh();
  _selfIntersectionVoxels=_selfIntersectionMM/_resolution;
  _maxStepSizeVoxels=_maxStepSizeMM/_resolution;
  _stepSizeSmoothingSig=_stepSizeSmoothingSigMM/_resolution;
}


void VolumeToMesh::InitializeFilters(){
  _normalFilter->SplittingOff();
  _normalFilter->AutoOrientNormalsOff() ;
  _normalFilter->ConsistencyOn();
  _remeshFilter->MaxEdgeLength(2*_averageEdgeLengthMM/_resolution);
  _remeshFilter->MinEdgeLength(0.5*_averageEdgeLengthMM/_resolution);
  _smoothingFilter->SetLambda(_smoothing);
  _smoothingFilter->SetSigma(0);
  _smoothingFilter->SetPointWeighting(_smoothArrayName);
  if (_gaussianInterpolationOn){
    _interpolator = new GaussianInterpolateImageFunction(_resolution);
  }
  else{
    _interpolator = new LinearInterpolateImageFunction();
  }
  _interpolator->Input(&_distanceField);
  _interpolator->Initialize();
  _lowPassFilter->SetPassBand(_lowPassFilterBand);
  _lowPassFilter->SetNumberOfIterations(_lowPassFilterIterations);
  _lowPassFilter->FeatureEdgeSmoothingOff();


}


void  VolumeToMesh::CalculateCentreOfMass(){
  int x,y,z,i,j,k,numVoxelsInside;

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();

  //determine com
  _centreOfMass[0]=0;
  _centreOfMass[1]=0;
  _centreOfMass[2]=0;
  numVoxelsInside=0;
  for (i=0;i<x;i++){
    for (j=0;j<y;j++){
      for (k=0;k<z;k++){
        if (_binaryVolume(i,j,k)==1){
          numVoxelsInside++;
          _centreOfMass[0]+=i;
          _centreOfMass[1]+=j;
          _centreOfMass[2]+=k;
        }
      }
    }
  }
  _centreOfMass[0]/=numVoxelsInside;
  _centreOfMass[1]/=numVoxelsInside;
  _centreOfMass[2]/=numVoxelsInside;

}

void InitializeParameters(){

}

void VolumeToMesh::InitializeMesh(){
  if (_initialMesh==hull) { InitializeMeshAsConvexHull(); }
  if (_initialMesh==sphere) { InitializeMeshAsBoundingSphere(); }
  if (_initialMesh==input) {   TransformPointsToImageCoordinates(); }
}


void VolumeToMesh::InitializeMeshAsBoundingSphere(){
  _mesh=vtkSmartPointer<vtkPolyData>::New();
  ParticleSampleSphere pss;
  int x,y,z,i,j,k,numVoxelsInside;

  int minSamples,subdivisions,minSubdivisions;
  double boundingRadius,distance,furthestDistance,surfaceArea, phiAngle, thetaAngle;
  int numPoints,samples;
  vtkSmartPointer<vtkPoints> points;
  double point [3], furthestPoint [3];

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();


  furthestDistance=0;
  for (i=0;i<x;i++){
    for (j=0;j<y;j++){
      for (k=0;k<z;k++){
        if (_binaryVolume(i,j,k)==1){
          distance = sqrt(
          pow((double)i-_centreOfMass[0],2.0) +
          pow((double)j-_centreOfMass[1],2.0) +
          pow((double)k-_centreOfMass[2],2.0));
          if (distance > furthestDistance){
            furthestPoint[0]=i;
            furthestPoint[1]=j;
            furthestPoint[2]=k;
            furthestDistance=distance;
          }
        }
      }
    }
  }



  _binaryVolume.ImageToWorld(_centreOfMass[0],_centreOfMass[1],_centreOfMass[2]);
  _binaryVolume.ImageToWorld(furthestPoint[0],furthestPoint[1],furthestPoint[2]);
  boundingRadius = sqrt(
      pow(furthestPoint[0]-_centreOfMass[0],2) +
      pow(furthestPoint[1]-_centreOfMass[1],2) +
      pow(furthestPoint[2]-_centreOfMass[2],2));


  cout << "centre of mass: " << _centreOfMass[0] << "," << _centreOfMass[1] << "," << _centreOfMass[2] << endl;
  cout << "furthestPoint: " << furthestPoint[0] << "," << furthestPoint[1] << "," << furthestPoint[2] << endl;

  cout << "radius: " << boundingRadius << endl;

  surfaceArea=4*M_PI*boundingRadius*boundingRadius;
  cout << "surfaceArea: " << surfaceArea << endl;
  //F=2*SA/(EL^2) (number of faces based on equilateral triangles)
  //V−E+F=2(1-g) (Euler's formula)
  //g=1 (genus: number of handles)
  //V=2+E-F
  //2E=3F (Each face has 3 half edges)
  //V=2+F/2;
  //V=2+SA/(EL^2);

  numPoints=(int)(0.5+2+surfaceArea/(_averageEdgeLengthMM*_averageEdgeLengthMM));
  cout << numPoints << " points" << endl;
  //work out number of subdivisions and number of initial points


  minSamples=25; //max samples with be 4*minSamples

  //calc the max number of subdivisions s based on the minimum sample number v
  // F=f*4^s=f*2^2s
  // F=2*(V-2)
  // 2(V-2)=2(v-2)*2^2s
  // 2^2s=(V-2)/(v-2)
  // s=log2((V-2)/(v-2))/2
  // take floor
  subdivisions=(int)(log2((numPoints-2)/(minSamples-2))/2);
  //calc initial samples based on
  //2^2s=(V-2)/(v-2)
  //v=((V-2)/2^2s)+2


  samples=(int)(0.5+(((double)numPoints-2.0)/pow(4.0,subdivisions))+2);
  //recalc number of points
  //V=(2^2s)*(v-2)
  numPoints=pow(4.0,subdivisions)*(minSamples-2);
  pss.SetNumberOfPoints(samples);


  if (_levels>1){
    subdivisions-=_levels-1;
    _averageEdgeLengthMM*=pow(2.0,_levels-1);
  }

  pss.SetNumberOfSubdivisions(subdivisions);
  pss.SphericalCoordsOn();
  _mesh=pss.GetSphere();




  numPoints=_mesh->GetNumberOfPoints();

  points=_mesh->GetPoints();
  //Scale, Translate vertices and convert back to image coords
  for (i=0; i<numPoints; i++){
    points->GetPoint(i,point);
    //cout << point[0] << "," << point[1] << "," << point[2] << endl;
    for (j=0; j<3; j++)  point[j]=point[j]*boundingRadius+_centreOfMass[j];
    //_binaryVolume.WorldToImage(point[0],point[1],point[2]);
    points->SetPoint(i,point);
    //cout << point[0] << "," << point[1] << "," << point[2] << endl;
  }


  _binaryVolume.WorldToImage(_centreOfMass[0],_centreOfMass[1],_centreOfMass[2]);

  points->GetPoint(0,point);
  _mesh->SetPoints(points);
  // _mesh->Update();
  _mesh->Modified();
  cout << "scaled" << boundingRadius << endl;
  TransformPointsToImageCoordinates();
}


void VolumeToMesh::InitializeMeshAsConvexHull(){
  LaplacianSmoothFilter smoothingFilter;
  vtkSmartPointer<vtkMarchingCubes> mcubes
    = vtkSmartPointer<vtkMarchingCubes>::New();
  vtkSmartPointer<vtkImageData> vtkImage;

  vtkSmartPointer<vtkPolyData> hull = GetHull();

  vtkImage=Mesh2VtkMask(hull);
  mcubes->SetValue(0, 0.5);
  mcubes->ComputeNormalsOn();
  mcubes->ComputeGradientsOff();
  mcubes->SetInputData(vtkImage);
  mcubes->Update();
  _mesh = mcubes->GetOutput();
  //smooth hull

  smoothingFilter.SetLambda(1);
  smoothingFilter.SetSigma(0);
  for (int i=0;i<20;i++){
    smoothingFilter.SetInput(_mesh);
    smoothingFilter.Update();
    _mesh=smoothingFilter.GetOutput();
  }

  ComputeMeshNormals();
  _mesh->BuildCells();

/*
  // Write result

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(_mesh);
  writer->SetFileName("hull.vtk");
  writer->SetFileTypeToBinary();
  writer->Write();
  writer->Delete();

  exit(1);

*/


}




void VolumeToMesh::InitializeDistanceField(){
  RealImage input;
  _edt = new EuclideanDistanceTransform<RealPixel>(EuclideanDistanceTransform<RealPixel>::DT_3D);


  // Threshold image
  RealImage insideDistanceField;
  int x,y,z;

  input=(RealImage)_binaryVolume;
  _distanceField=input;

  cout << "Calculating distance map...";
  cout.flush();
  // Calculate EDT
  _edt->Input(&input);
  _edt->Output(&_distanceField);
  _edt->Run();

  //get inverse
  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
         input(x, y, z)=input(x, y, z)==0;
        }
      }
    }

  _edt->Input(&input);
  _edt->Output(&insideDistanceField);
  _edt->Run();

  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
        _distanceField(x, y, z)  = sqrt(_distanceField(x, y, z)) - sqrt(insideDistanceField(x, y, z));
      }
    }
  }

  //_distanceField.Write("dist_field.nii.gz");
  //_binaryVolume.Write("binary_volume.nii.gz");
  cout << "done" << endl;
}


//convert binary image to vtk point field.
//step 1 erode
//step 2 difference

vtkSmartPointer<vtkPolyData> VolumeToMesh::BinaryVolume2VTKSurfacePoints(){
  int i,j,k,x,y,z;
  GreyImage surfaceVoxels;
  Erosion<GreyPixel> erosion;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  double point [3];

  //Establish surface voxels
  erosion.Connectivity(CONNECTIVITY_26);

  surfaceVoxels = _binaryVolume;

  erosion.Input(&surfaceVoxels);
  erosion.Output(&surfaceVoxels);
  erosion.Run();
  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();
  for (i=0;i<x;i++){
    for (j=0;j<y;j++){
      for (k=0;k<z;k++){
        if (_binaryVolume(i,j,k) && ~surfaceVoxels(i,j,k)){
          points->InsertNextPoint(i,j,k);
        }
      }
    }
  }

  polydata->SetPoints(points);

  return polydata;
}

vtkSmartPointer<vtkImageData> VolumeToMesh::Mesh2VtkMask(vtkSmartPointer<vtkPolyData> polydata){
  vtkSmartPointer<vtkTriangleFilter> triFilter = vtkTriangleFilter::New();
  GreyImage binaryVolume;
  GreyPixel *ImagePtr;
  double resampleFactor,resampleFactorX,resampleFactorY,resampleFactorZ;

  int x,y,z,numVoxels;
  int vtkX,vtkY,vtkZ;
  int i, j, k;
  triFilter->SetInputData(polydata);
  triFilter->Update();
  polydata = triFilter->GetOutput();

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();


  vtkSmartPointer<vtkImageData> vtkImage
    = vtkSmartPointer<vtkImageData>::New();

  resampleFactor=_averageEdgeLengthMM/_resolution;
  vtkX=ceil((int)(x/resampleFactor));
  vtkY=ceil((int)(y/resampleFactor));
  vtkZ=ceil((int)(z/resampleFactor));
  numVoxels=vtkX*vtkY*vtkZ;
  vtkImage->SetDimensions(vtkX,vtkY,vtkZ);

  resampleFactorX=(float)x/(float)vtkX;
  resampleFactorY=(float)y/(float)vtkY;
  resampleFactorZ=(float)z/(float)vtkZ;

  vtkImage->SetSpacing(resampleFactorX,resampleFactorY,resampleFactorZ);
  vtkImage->SetOrigin(0,0,0);
  vtkImage->AllocateScalars(VTK_SHORT, 1);
  // vtkImage->SetScalarTypeToShort();
  // vtkImage->SetNumberOfScalarComponents(1);
  // vtkImage->AllocateScalars();

  short *vtkImagePtr= (short *)(vtkImage->GetScalarPointer());

  for(i = 0; i < numVoxels; i++) {
    *vtkImagePtr=0;
    vtkImagePtr++;
  }

  vtkSmartPointer<vtkPolyDataToImageStencil> dataToStencil
    = vtkPolyDataToImageStencil::New();
  dataToStencil->SetTolerance(0);
  dataToStencil->SetInputData(polydata);
  dataToStencil->SetInformationInput(vtkImage);

  vtkSmartPointer<vtkImageStencil> stencil = vtkImageStencil::New();
  stencil->SetInputData(vtkImage);
  stencil->SetStencilData(dataToStencil->GetOutput());
  stencil->SetBackgroundValue(1);
  stencil->Update();
  vtkImage = stencil->GetOutput();
  //vtkImage->Modified();
  //vtkImage->Update();

  return vtkImage;
}



vtkSmartPointer<vtkPolyData> VolumeToMesh::GetHull(){



  int levels = 3;
  vtkSmartPointer<vtkHull> hullFilter;
  vtkSmartPointer<vtkDelaunay3D> delaunay;
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter;
  vtkSmartPointer<vtkPolyData> surfacePoints;
  vtkSmartPointer<vtkPolyData> hull;
  hullFilter= vtkSmartPointer<vtkHull>::New();
  delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
  surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();


  surfacePoints=BinaryVolume2VTKSurfacePoints();




  hullFilter->SetInputData(surfacePoints);
  hullFilter->AddRecursiveSpherePlanes(levels);
  delaunay->SetInputConnection(hullFilter->GetOutputPort());
  surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
  surfaceFilter->Update();
  hull=surfaceFilter->GetOutput();

  /*
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(hull);
  writer->SetFileName("hull.vtk");
  writer->SetFileTypeToBinary();
  writer->Write();



  cout << "returning hull" << endl;
  */
  return hull;
}


void VolumeToMesh::SetInitialMesh (vtkSmartPointer<vtkPolyData> polydata){
  _mesh=polydata;
  _initialMesh=input;
}

void VolumeToMesh::SetInitialMeshToSphere(){
  _initialMesh=sphere;
}
void VolumeToMesh::SetInitialMeshToHull(){
  _initialMesh=hull;
}


void VolumeToMesh::SetInput (RealImage *in){
  //TODO resize non-isotropic volumes
  Matrix i2w;
  cout << _boundaryThreshold << endl;
  cout << "Setting Input" << endl;
  _binaryVolume=*in;

  int x,y,z;
  for (z = 0; z < in->GetZ(); z++) {
    for (y = 0; y < in->GetY(); y++) {
      for (x = 0; x < in->GetX(); x++) {
        _binaryVolume(x, y, z) = in->Get(x,y,z)>_boundaryThreshold;
      }
    }
  }
  if(_closeBinaryVolume){
    cout << "Closing mask...";
    cout.flush();
    ///Close mask
    Erosion<GreyPixel> erosion;
    erosion.Connectivity(CONNECTIVITY_26);
    Dilation<GreyPixel> dilation;
    dilation.Connectivity(CONNECTIVITY_26);

    dilation.Input(&_binaryVolume);
    dilation.Output(&_binaryVolume);
    dilation.Run();

    erosion.Input(&_binaryVolume);
    erosion.Output(&_binaryVolume);
    erosion.Run();
    cout << "done" << endl;
  }
  _resolution=_binaryVolume.GetXSize();
  i2w = _binaryVolume.GetImageToWorldMatrix();
  _flipNormals=i2w(0,0)*i2w(1,1)*i2w(2,2)>0;
  CalculateCentreOfMass();
}

bool VolumeToMesh::Intersect(double vector [], double point [] ){
  //check not outside of volume!

  double normVector [3];
  int i,x,y,z;
  double nf;
  bool searchBackwards,noIntersection, outOfBounds,outOfBounds_;

  int pointRounded [3];

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();

  for (i=0;i<3;i++)  pointRounded[i]=floor(point[i]+0.5);

  //if starting point is inside volume then return true
  //if ( pointRounded[0]>=0 && pointRounded[0]<x && pointRounded[1]>=0 && pointRounded[1]<y && pointRounded[2]>=0 && pointRounded[2]<z ){
  //  inBounds=false;
    if (_binaryVolume(pointRounded[0],pointRounded[1],pointRounded[2])) //{
      return true;





  int pointRounded_ [3];
  double vectorMag [3],searchPoint [3],searchPoint_ [3];




  for (i=0;i<3;i++)  vectorMag[i]=abs(vector[i]);

  if (vectorMag[0] > vectorMag[1] ){
    if (vectorMag[0] > vectorMag[2] )  nf=vectorMag[0];
    else nf=vectorMag[2];
  }
  else if(vectorMag[1] > vectorMag[2])
    nf=vectorMag[1];
  else
    nf=vectorMag[2];

  for (i=0;i<3;i++)  normVector[i]=vector[i]/nf;


  //step along normal direction
  for (i=0;i<3;i++){
    searchPoint[i]=point[i]+normVector[i];
    searchPoint_[i]=point[i]-normVector[i];
  }






  noIntersection=true;
  outOfBounds=false;
  outOfBounds_=false;

  //while no intersection and outside

  while (noIntersection){


    for (i=0;i<3;i++)  pointRounded[i]=floor(searchPoint[i]+0.5);
    if (pointRounded[0]>=x || pointRounded[0]<0 ||
        pointRounded[1]>=y || pointRounded[1]<0 ||
        pointRounded[2]>=z || pointRounded[2]<0){
      break; //out of bounds
    }
    else{
      if(_binaryVolume(pointRounded[0],pointRounded[1],pointRounded[2])){ //intersection detected!
        noIntersection=false;
        //cout << "binVolume" << _binaryVolume(0,0,0) << endl;
        //cout << "intersection" << endl;
      }
      else{
         //step forwards along normal  direction
        //cout << "step" << endl;
        for (i=0;i<3;i++) searchPoint[i]=searchPoint[i]+normVector[i];
      }
    }


    if (!outOfBounds_){// check backwards direction
      for (i=0;i<3;i++)  pointRounded_[i]=floor(searchPoint_[i]+0.5);
      if (pointRounded_[0]>=x || pointRounded_[0]<0 ||
          pointRounded_[1]>=y || pointRounded_[1]<0 ||
          pointRounded_[2]>=z || pointRounded_[2]<0){
        outOfBounds_=true;
      }
      else{
        if(_binaryVolume(pointRounded_[0],pointRounded_[1],pointRounded_[2])){ //intersection detected!
          break;
          //cout << "intersection" << endl;
        }
        else{
           //step backwards along normal  direction
          //cout << "step" << endl;
          for (i=0;i<3;i++)  searchPoint_[i]=searchPoint_[i]-normVector[i];
        }
      }
    }
  }


  return !noIntersection;

}




double VolumeToMesh::GetVector(double p1 [],double p2 [], double vector []){
  double nf=0;
  int j;
  for(j=0;j<3;j++){
    vector[j]=p2[j]-p1[j]; //head towards com
    nf+=pow(vector[j],2);
  }
  nf=sqrt(nf);
  for(j=0;j<3;j++) vector[j]=vector[j]/nf;
  return nf;
}



void printVector(double v []){
  for (int i = 0; i < 3; ++i){
    cout << v[i];
    if (i != 2) { cout << ","; }
  }
  cout << endl;
}


void copy3DVector(double*v, double*v2 ){
  for (int i = 0; i < 3; i++){
    v2[i]=v[i];
  }
}

void invert3DVector(double*v, double*v2 ){
  for (int i = 0; i < 3; ++i){
    v2[i]=-v[i];
  }
}

void normalize3DVector(double*v, double*v2 ){
  double nf;
  for (int i = 0; i < 3; ++i){
    nf+=pow(v[i],2);
  }
  nf=sqrt(nf);
  for (int i = 0; i < 3; ++i){
    v2[i]=v2[i]/nf;
  }

}



void VolumeToMesh::ComputeMeshNormals(){
  vtkSmartPointer<vtkPolyDataNormals> normalFilter =
      vtkSmartPointer<vtkPolyDataNormals>::New();
  normalFilter->SetInputData(_mesh);
  normalFilter->SplittingOff();
  normalFilter->AutoOrientNormalsOn();
  normalFilter->ConsistencyOn();
  normalFilter->Update();
  _mesh = normalFilter->GetOutput();
  //_mesh = NULL;
  //_mesh = vtkSmartPointer<vtkPolyData>::New();
  //_mesh->DeepCopy(normalFilter->GetOutput());
}


void VolumeToMesh::SmoothMesh(){
  LaplacianSmoothFilter smoothingFilter;
  //_mesh->SetSource(0);
  smoothingFilter.SetLambda(1);
  smoothingFilter.SetPointWeighting(_smoothArrayName);
  smoothingFilter.SetSigma(0);
  smoothingFilter.SetInput(_mesh);
  smoothingFilter.Update();
  _mesh = NULL;
  _mesh = vtkSmartPointer<vtkPolyData>::New();
  _mesh->DeepCopy(smoothingFilter.GetOutput());
   //cout << "smoothing mesh: " << _smoothing << endl;
}


void VolumeToMesh::TransformPointsToWorldCoordinates(){
  int ptIdx,numPoints;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkReverseSense> reverseCells;
  double normal [3],vertex [3];

  points=_mesh->GetPoints();
  numPoints=_mesh->GetNumberOfPoints();

  for(ptIdx=0;ptIdx<numPoints;ptIdx++){
    points->GetPoint(ptIdx,vertex);
    _binaryVolume.ImageToWorld(vertex[0],vertex[1],vertex[2]);
    points->SetPoint(ptIdx,vertex);
  }
  _mesh->SetPoints(points);


  if (_flipNormals) {
    reverseCells=  vtkSmartPointer<vtkReverseSense>::New();
    reverseCells->ReverseCellsOn();
    reverseCells->ReverseNormalsOff();
    //_mesh->SetSource(0);
    reverseCells->SetInputData(_mesh);
    reverseCells->Update();
    _mesh=reverseCells->GetOutput();
  }
  ComputeMeshNormals();

}

void VolumeToMesh::TransformPointsToImageCoordinates(){
  int ptIdx,numPoints;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkReverseSense> reverseCells;
  double vertex [3];

  points=_mesh->GetPoints();
  numPoints=_mesh->GetNumberOfPoints();

  for(ptIdx=0;ptIdx<numPoints;ptIdx++){
    points->GetPoint(ptIdx,vertex);
    _binaryVolume.WorldToImage(vertex[0],vertex[1],vertex[2]);
    points->SetPoint(ptIdx,vertex);
  }
  _mesh->SetPoints(points);

  if (_flipNormals) {
    reverseCells=  vtkSmartPointer<vtkReverseSense>::New();
    reverseCells->ReverseCellsOn();
    reverseCells->ReverseNormalsOff();
    //_mesh->SetSource(0);
    reverseCells->SetInputData(_mesh);
    reverseCells->Update();
    _mesh=reverseCells->GetOutput();
  }
  ComputeMeshNormals();
}

//freeze

bool VolumeToMesh::SelfIntersect(double vertex [3], double vector [3],
    double * stepsize){

  int i, subId=-1;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double squareDistance,t,tol;

  double projectionDist = _maxStepSizeVoxels*1.5 + *stepsize;

  vtkIdType cellId;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-projectionDist*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
     pcoords[i]=0;
     x[i]=0;
   }
   //outward intersection, possibly twistes do nothing
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward intersection
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+projectionDist*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){
     *stepsize=0;
     return true;
   }
   return false;
}



//constant repulsion
/*
bool VolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){
  int i, subId=-1;;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double squareDistance,t,tol;
  vtkIdType cellId;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-_selfIntersectionVoxels*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
   }
   //outward intersection, possibly twistes do nothing
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward intersection
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+_selfIntersectionVoxels*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){
     if (*stepsize>-0.25*_maxStepSizeVoxels) {
       *stepsize=-0.25*_maxStepSize;
     }
     return true;
   }
   else{
     return false;
   }
}
*/
//Progressoin repulsion
/*
bool VolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){
  int i, subId=-1;;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double distance,projectionLength,t,tol,maxstep;
  vtkIdType cellId;

  projectionLength=_selfIntersectionVoxels;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-projectionLength*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
   }
   //outward intersection, possibly gone to far already self-intersecting!
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward self intersection
   projectionLength=_selfIntersectionVoxels+_maxStepSizeVoxels;
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+projectionLength*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){//reset gap to _selfIntersectionVoxels by moving half of
     distance=projectionLength*t;

     if(distance<=_selfIntersectionVoxels){//too close :/
       *stepsize=-(_selfIntersectionVoxels-distance)/2; //step halfway to safety
     }
     else{ //in the danger zone, only allow half the step to the boundary
       maxstep=(distance-_selfIntersectionVoxels)/2;
       if (*stepsize>maxstep) {
         *stepsize=maxstep;
       }
     }
     return true;
   }

   return false;
}
*/

//Progressoin repulsion 2
/*
bool VolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){
  int i, subId=-1;;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double distance,projectionLength,t,tol,maxstep;
  vtkIdType cellId;

  projectionLength=_selfIntersectionVoxels;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-projectionLength*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
   }
   //outward intersection, possibly gone to far already self-intersecting!, do nothing
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward self intersection
   projectionLength=_selfIntersectionVoxels+_maxStepSizeVoxels;
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+projectionLength*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){//reset gap to _selfIntersectionVoxels by moving half of
     distance=projectionLength*t;

     if(distance<=_selfIntersectionVoxels){//too close :/
       *stepsize=-(_selfIntersectionVoxels-distance)/2; //step to safety
     }
//     else{ //in the danger zone, only allow half the step to the boundary
//       maxstep=(distance-_selfIntersectionVoxels)/2;
//       if ( *stepsize>maxstep) {
//         *stepsize=maxstep;
//       }
//     }

     return true;
   }


   return false;
}
*/


void VolumeToMesh::BuildSelfIntersectionTree(){
  //_mesh->SetSource(0);
  _selfIntersectBSPTree->SetDataSet(_mesh);
  _selfIntersectBSPTree->BuildLocator();
}

bool VolumeToMesh::InsideVolume(double vertex [3]){
  //check if outside volume
  bool inside;
  int x,y,z;
  x=_distanceField.GetX()-1;
  y=_distanceField.GetY()-1;
  z=_distanceField.GetZ()-1;
  inside=(vertex[0]>=0 && vertex[0]<=x &&
          vertex[1]>=0 && vertex[1]<=y &&
          vertex[2]>=0 && vertex[2]<=z);
  return inside;
}

//TODO refactor this into a seperate filter class
void smoothScalars(vtkPolyData * poly,char * name,  double kernel,
    int noOfIterations) {

  vtkIdType i;
  long j, k, n;
  bool ok;


  int noOfPoints, noOfCells;

  double v1[3], v2[3];

  vtkIdType *ptsInCell = NULL;
  vtkIdType noOfPointsInCell;
  int u, v;

  int noOfEdgesAtPt;
  double sumDist, dist, meanDist;
  int count;
  double sumVals, distSq;
  double w, sumW;

  noOfPoints = poly->GetNumberOfPoints();
  noOfCells  = poly->GetNumberOfCells();

  // Storage for edge information.  Really this is double the amount needed
  // as the edge information for each edge's endpoint are stored.
  vtkIdList **edges;

  edges = new vtkIdList *[noOfPoints];
  for (i = 0; i < noOfPoints; ++i){
    edges[i] = vtkIdList::New();
  }
  cerr << "Allocated edges " << endl;

  poly->BuildCells();
  poly->BuildLinks();
  // poly->Update();
  poly->Modified();

  // Find the adjacency information.
  for (i = 0; i < noOfCells; ++i){
    poly->GetCellPoints(i, noOfPointsInCell, ptsInCell);

    u = ptsInCell[0];
    v = ptsInCell[noOfPointsInCell-1];

    edges[u]->InsertUniqueId(v);
    edges[v]->InsertUniqueId(u);

    for (j = 0; j < noOfPointsInCell - 1 ; ++j){
      u = ptsInCell[j];
      v = ptsInCell[j + 1];

      edges[u]->InsertUniqueId(v);
      edges[v]->InsertUniqueId(u);
    }
  }

  // Calculate mean edge length. The following is a bit redundant as every edge will
  // appear twice in the calculation.
  sumDist = 0.0;
  count   = 0;

  for (i = 0; i < noOfPoints; ++i){
    noOfEdgesAtPt = edges[i]->GetNumberOfIds();

    poly->GetPoint(i, v1);

    for (j = 0; j < noOfEdgesAtPt; ++j){
      poly->GetPoint(edges[i]->GetId(j), v2);

      dist = sqrt(vtkMath::Distance2BetweenPoints(v1, v2));
      sumDist += dist;
      ++count;
    }
  }

  meanDist = sumDist / ((double) count);

  cerr << "Points      : " << noOfPoints << endl;
  cerr << "Cells       : " << noOfCells << endl;
  cerr << "Edges       : " << count / 2 << endl;
  cerr << "Mean length : " << meanDist << endl;


  double denom = 2 * kernel * kernel * meanDist * meanDist;

  cerr << "Iterating " << endl;


  // Get the current scalars.
  vtkFloatArray * scalarsIn
      = (vtkFloatArray *) poly->GetPointData()->GetScalars(name);
  vtkSmartPointer<vtkFloatArray> scalarsOut =
      vtkSmartPointer<vtkFloatArray>::New();
  scalarsOut->DeepCopy(scalarsIn);

  for (n = 0; n < noOfIterations; ++n){
    cerr << "."; cerr.flush();


    for (i = 0; i < noOfPoints; ++i){
      poly->GetPoint(i, v1);

      noOfEdgesAtPt = edges[i]->GetNumberOfIds();
      // Initialise by including the current point.
      sumVals = scalarsIn->GetTuple1(i);
      sumW    = 1.0;

      for (j = 0; j < noOfEdgesAtPt; ++j){
        k = edges[i]->GetId(j);
        poly->GetPoint(k, v2);

        distSq = vtkMath::Distance2BetweenPoints(v1, v2);
        w = exp(-1.0 * distSq / denom);

        // Add adjoining point's contribution.
        sumVals += w * scalarsIn->GetTuple1(k);
        sumW    += w;
      }
      // Updated value for current point.
      scalarsOut->SetTuple1(i, sumVals / sumW);
    }

    // Copy for next iteration.
    for (i = 0; i < noOfPoints; ++i){
      scalarsIn->SetTuple1(i, scalarsOut->GetTuple1(i));
    }

  }





  for (i = 0; i < noOfPoints; ++i){
    edges[i]->Delete();
  }
  delete edges;

}



double VolumeToMesh::UpdatePointPositions(){
  bool update;
  int i,ptId;
  int numPoints;
  char ch;
  double n [3],direction [3],p [3],d[3];
  double stepSize,cumulativeDistance;
  char * stepsizeArrayName = "_stepsize";

  double smooth;

  vtkSmartPointer<vtkDataArray> normals;
  vtkSmartPointer<vtkPoints> points, newPoints;



  numPoints=_mesh->GetNumberOfPoints();
  normals=_mesh->GetPointData()->GetNormals();


  if (_mesh->GetPointData()->HasArray(_smoothArrayName)){
    _mesh->GetPointData()->RemoveArray(_smoothArrayName);
  }

  vtkSmartPointer<vtkDoubleArray> smoothingWeight =
      vtkSmartPointer<vtkDoubleArray>::New();
  smoothingWeight->SetNumberOfComponents(1);
  smoothingWeight->SetNumberOfTuples(numPoints);
  smoothingWeight->SetName(_smoothArrayName);

  vtkSmartPointer<vtkDoubleArray> stepsizeArray =
      vtkSmartPointer<vtkDoubleArray>::New();
  stepsizeArray->SetNumberOfComponents(1);
  stepsizeArray->SetNumberOfTuples(numPoints);
  stepsizeArray->SetName(stepsizeArrayName);
  _mesh->GetPointData()->AddArray(stepsizeArray);
  cumulativeDistance = 0;

  if (_selfIntersectionOn) { BuildSelfIntersectionTree(); }
  for(ptId=0;ptId<numPoints;ptId++){

    _mesh->GetPoint(ptId,p); //Get p position
    normals->GetTuple(ptId,n); //Get normal vector
    for (i=0; i<3; i++){ d[i] = -n[i]; }

    smooth = 0; // _smoothing;

    if (!InsideVolume(p)){
      stepSize=_maxStepSizeVoxels;
    }
    else{
      stepSize=_interpolator->Evaluate(p[0], p[1], p[2]);
      cumulativeDistance+=abs(stepSize);
      if (stepSize > _maxStepSizeVoxels) stepSize = _maxStepSizeVoxels;
      if (stepSize < -_maxStepSizeVoxels) stepSize = -_maxStepSizeVoxels;
      if (_intersectionOn && ! Intersect(d,p)==1){
        smooth = 1; //_smoothing;
        stepSize=0;
      }

    }
    cumulativeDistance+=stepSize;

    if(_selfIntersectionOn && stepSize > 0 ){ //modify stepsize if in the proxmity itself
      smooth = SelfIntersect(p,d,&stepSize);
    }

    stepsizeArray->SetTuple(ptId,&stepSize);
    smoothingWeight->SetTuple1(ptId,smooth);
  }

  if (_stepSizeSmoothingIterations > 0){
    smoothScalars(_mesh,stepsizeArrayName, _stepSizeSmoothingSig,
        _stepSizeSmoothingIterations);
  }


  for(ptId=0;ptId<numPoints;ptId++){
    normals->GetTuple(ptId,n); //Get normal vector
    stepSize = stepsizeArray->GetTuple1(ptId);
    _mesh->GetPoint(ptId,p); //Get p position
    for (i=0; i<3; i++){ p[i] = p[i] - n[i]*stepSize; }
    _mesh->GetPoints()->SetPoint(ptId,p);
  }

  _mesh->GetPointData()->RemoveArray(stepsizeArrayName);
  _mesh->GetPointData()->AddArray(smoothingWeight);
  return cumulativeDistance/(double)numPoints;
}



void VolumeToMesh::LowPassFilter(){

  vtkSmartPointer<vtkWindowedSincPolyDataFilter> sincFilter =
      vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();

  double translation [3], scaleFactor [3], p [3];
  double newCentroid [3], centroid [3];
  double newScale [3], scale [3];
  int numPoints = _mesh->GetNumberOfPoints();

  //calculate "scale" and centroid
  for (int i=0; i<3;i++) centroid[i] = 0;
  for (int i=0; i<numPoints;i++){
    _mesh->GetPoint(i,p);
    for (int i=0; i<3;i++) centroid[i] += p[i];
  }
  for (int i=0; i<3;i++) centroid[i] /= numPoints;

  for (int i=0; i<3;i++) scale[i] = 0;
  for (int i=0; i<numPoints;i++){
    _mesh->GetPoint(i,p);
    for (int i=0; i<3;i++) scale[i] += abs(p[i]-centroid[i]);
  }
  for (int i=0; i<3;i++) scale[i] /= numPoints;


  sincFilter->SetPassBand(_lowPassFilterBand);
  sincFilter->SetNumberOfIterations(_lowPassFilterIterations);
  sincFilter->NormalizeCoordinatesOn();
  //_mesh->SetSource(0);
  sincFilter->SetInputData(_mesh);
  sincFilter->Update();
  _mesh=sincFilter->GetOutput();

  //calculate new "scale" and centroid
  for (int i=0; i<3;i++) newCentroid[i] = 0;
  for (int i=0; i<numPoints;i++){
    _mesh->GetPoint(i,p);
    for (int i=0; i<3;i++) newCentroid[i] += p[i];
  }
  for (int i=0; i<3;i++) newCentroid[i] /= numPoints;

  for (int i=0; i<3;i++) newScale[i] = 0;
  for (int i=0; i<numPoints;i++){
    _mesh->GetPoint(i,p);
    for (int i=0; i<3;i++) newScale[i] += abs(p[i]-newCentroid[i]);
  }
  for (int i=0; i<3;i++) newScale[i] /= numPoints;

  //fix translation

  for (int i=0; i<3;i++) translation[i] = newCentroid[i] - centroid[i];

  for (int i=0; i<numPoints;i++){
    _mesh->GetPoint(i,p);
    for (int i=0; i<3;i++) p[i] -= translation[i];
    _mesh->GetPoints()->SetPoint(i,p);
  }

  //fix scale

  for (int i=0; i<3;i++) scaleFactor[i] = newScale[i] / scale[i];

  for (int i=0; i<numPoints;i++){
    _mesh->GetPoint(i,p);
    for (int i=0; i<3;i++)
      p[i] = ((p[i] - centroid[i]) / scaleFactor[i]) + centroid[i];
    _mesh->GetPoints()->SetPoint(i,p);
  }

}



void VolumeToMesh::Subdivide(){
  vtkSmartPointer<vtkButterflySubdivisionFilter> subdivisionFilter;

  subdivisionFilter = vtkSmartPointer<vtkButterflySubdivisionFilter>::New();
  //_mesh->SetSource(0);
  subdivisionFilter->SetInputData(_mesh);
  subdivisionFilter->Update();
  _mesh=subdivisionFilter->GetOutput();
  _averageEdgeLengthMM=_averageEdgeLengthMM/2;
  _remeshFilter->MaxEdgeLength((_averageEdgeLengthMM/_resolution)*2);
  _remeshFilter->MinEdgeLength(_averageEdgeLengthMM/_resolution/2);
  _smoothing=_smoothing*2;
}

void VolumeToMesh::Screenshot(string fileName){
  // Visualize
  const char * fileNameChars=fileName.c_str();


    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(_mesh);

    vtkSmartPointer<vtkActor> actor =
      vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetAlphaBitPlanes(1); //enable usage of alpha channel
    renderWindow->SetOffScreenRendering(1);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    //renderer->SetBackground(1,1,1); // Background color white


    renderWindow->Render();

    // Screenshot
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(fileNameChars);
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();


}

/*
void VolumeToMesh::Screenshot(string fileName){
  const char * fileNameChars=fileName.c_str();



  vtkSmartPointer<vtkGraphicsFactory> graphics_factory =
    vtkSmartPointer<vtkGraphicsFactory>::New();
  graphics_factory->SetOffScreenOnlyMode( 1);
  //graphics_factory->SetUseMesaClasses( 1 );


    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(_mesh);

    vtkSmartPointer<vtkActor> actor =
      vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);


    // A renderer and render window
    vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetOffScreenRendering(1);
    renderWindow->AddRenderer(renderer);


    // Add the actors to the scene
    renderer->AddActor(actor);
    renderer->SetBackground(1,1,1); // Background color white
    renderWindow->Render();

    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(fileNameChars);
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();
}
*/

void VolumeToMesh::SetInputAsTarget(){
  //TODO set binary image to binary image!
}


void VolumeToMesh::DeformMeshToTarget(){
  // no self intersections
  //coarse tree for intersections with mesh
  int i,j,iteration, improvementIterations=10;
  double improvement,oldDistanceError,distanceError;
  double distanceErrors [improvementIterations];

  ostringstream os;
  string screenshotFileName;

  distanceError = 0;
  oldDistanceError=DBL_MAX;
  improvement=DBL_MAX;

  iteration = 0;
  while( iteration < _minIterations ||  (iteration < _maxIterations && improvement > _improvementThreshold) ){
    iteration++;
    cout << "iteration " << iteration << endl;
    ComputeMeshNormals();

    distanceError = UpdatePointPositions();


    //cout << "point positions updated" << endl;
    if (_smoothingOn) { SmoothMesh(); } //     cout << "smoothing done" << endl; }

    if (_lowPassFilterOn) { LowPassFilter(); } //    cout << "low pass filter done" << endl; }

    Remesh();

    if (_screenShotsOn){
      os << _screenShotDir << iteration << ".png";
      screenshotFileName = os.str();
      os.clear();
      os.str(std::string());
      Screenshot(screenshotFileName);
    }

    //oldDistanceError=distanceErrors[(iteration+1)%improvementIterations];
    distanceErrors[iteration%improvementIterations]=distanceError;
    cout << "distanceError: " << distanceError << endl;


    if (iteration>=improvementIterations){
      improvement=-distanceError;
      for (i=1;i<improvementIterations/2;i++){
        improvement-=distanceErrors[(iteration-i)%improvementIterations];
      }
      for (;i<improvementIterations;i++){
        improvement+=distanceErrors[(iteration-i)%improvementIterations];
      }
      improvement/=improvementIterations;
      //cout << "improvement: " << improvement << endl;
    }

  }

}


//original remeshing
void VolumeToMesh::Remesh(){
  Remesher remeshFilter;
  remeshFilter.Input(_mesh);
  remeshFilter.MaxEdgeLength(2*_averageEdgeLengthMM/_resolution);
  remeshFilter.MinEdgeLength(0.8*_averageEdgeLengthMM/_resolution);
  remeshFilter.Verbose(true);
  remeshFilter.Run();
  vtkSmartPointer<vtkPolyData> remeshed = remeshFilter.Output();

  _mesh = NULL;
  _mesh = vtkSmartPointer<vtkPolyData>::New();
  _mesh->DeepCopy(remeshed);
  _mesh = remeshed;
}




//MIRTK Remeshing
  /*
void VolumeToMesh::Remesh(){
  SurfaceRemeshing remesher;

  // remesher.SkipTriangulation(false);
  remesher.MinEdgeLength(2*_averageEdgeLengthMM/_resolution);
  remesher.MaxEdgeLength(0.8*_averageEdgeLengthMM/_resolution);
  remesher.Input(_mesh);
  remesher.Run();
  vtkSmartPointer<vtkPolyData> remeshed = remesher.Output();

  _mesh = NULL;
  _mesh = vtkSmartPointer<vtkPolyData>::New();
  _mesh->DeepCopy(remeshed);
  _mesh = remeshed;
}*/



void VolumeToMesh::FinalRefinement(){
  bool lowPassFilterOn=_lowPassFilterOn;
  bool lowPassFilterIterations=_lowPassFilterIterations;
  bool lowPassFilterBand=_lowPassFilterBand;
  int maxIterations=_maxIterations;
  _maxIterations=_refinementIterations;

  _lowPassFilterOn=true;
  _lowPassFilterIterations=25;
  _lowPassFilterBand=0.5;
  _smoothingOn=false;

  // set laplacian smoothing off

  cout << "-------------" << endl;
  cout << "Refining Mesh" << endl;
  cout << "-------------" << endl;
  DeformMeshToTarget();
  //restore filter settings
  //_smoothingFilter->SetLambda(_smoothing);


  _lowPassFilterOn=lowPassFilterOn;
  _lowPassFilterIterations=lowPassFilterIterations;
  _lowPassFilterBand=lowPassFilterBand;
  _smoothingOn=true;
  _maxIterations=maxIterations;

}


void VolumeToMesh::PrintSettings(){
  cout << "SETTINGS: " << endl;
  cout << "Edge length: " << _averageEdgeLengthMM << endl;
  cout << "Smoothing: " << _smoothing << endl;

  cout << "Maximum iterations: " << _maxIterations << endl;
  cout << "Epsilon: " << _improvementThreshold << endl;

  if (_selfIntersectionOn){
    cout << "Self intersection test: On" << endl;
  }
  else{
    cout << "Self intersection test: Off" << endl;
  }
  if (_lowPassFilterOn){
    cout << "Low pass filter On: " << _lowPassFilterIterations << " iterations, filtering band: " <<
        _lowPassFilterBand << endl;
  }

  if (_finalRefinementOn){
    cout << "Final refinement: On, smoothing: " << _refinementSmoothing << endl;
  }

}




vtkSmartPointer<vtkPolyData>  VolumeToMesh::GetOuput(){
  int iteration,level;
  double stepSize,maxStepSize,distanceError,nf;
  ostringstream os;
  string screenshotFileName;
  Initialize();
  PrintSettings();
  cout << "------------------------" << endl;
  cout << "Deforming mesh to target" << endl;
  cout << "------------------------" << endl;

  DeformMeshToTarget();
  if (_finalRefinementOn){ FinalRefinement(); }

  TransformPointsToWorldCoordinates();
  return _mesh;
}

} //namespace mirtk