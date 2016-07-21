/*
 * mirtkM2SOptimizer.h
 *
 *  Created on: 8 Feb 2016
 *      Author: raw11
 */

#ifndef MIRTKM2SOPTIMIZER_H_
#define MIRTKM2SOPTIMIZER_H_


#include "mirtk/Array.h"
#include "mirtk/Object.h"
#include "mirtk/EdgeTable.h"

//my includes
#include "mirtk/M2SEdge.h"


//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

namespace mirtk {

class M2SOptimizer {

  private:
    struct SphericalVector
    {
      double theta = 0.0; //inclination
      double phi = 0.0; //azimuthal angle
    };

  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);
  mirtkPublicAttributeMacro(Array<double>, CurrentSquareDistances);
  mirtkPublicAttributeMacro(Array<double>, Errors);
  mirtkPublicAttributeMacro(double, TotalWeight);
  mirtkPublicAttributeMacro(Array<SphericalVector>, OriginalSphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, OldSphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, SphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, BestSphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, Gradient);
  mirtkPublicAttributeMacro(Array<SphericalVector>, StepSizeW);

  mirtkPublicAttributeMacro(Array<bool>, StationaryPts);


  mirtkPublicAttributeMacro(int, LowPassIterations);
  mirtkPublicAttributeMacro(double, LowPassBandwidth);



  mirtkPublicAttributeMacro(Array<Edge>, Edges);
  mirtkPublicAttributeMacro(double, Energy);
  mirtkPublicAttributeMacro(int, MaxIterations);
  mirtkPublicAttributeMacro(int, NumEdges);
  mirtkPublicAttributeMacro(int, NumPts);


  mirtkPublicAttributeMacro(double, MinStepSize);
  mirtkPublicAttributeMacro(double, MaxStepSize);

  mirtkPublicAttributeMacro(bool, Debug);



  // Construction/Destruction
 public:

   /// Constructor
   M2SOptimizer();

   /// Copy constructor
   M2SOptimizer(const M2SOptimizer &);

   /// Assignment operator
   M2SOptimizer &operator =(const M2SOptimizer &);

   /// Copy attributes of this class from another instance
   void CopyAttributes(const M2SOptimizer &);

   /// Destructor
   virtual ~M2SOptimizer();



  void Initialize();
  void Run();
  void GetOutput();
  void UpdateSphericalCoordinatesFromMesh();
  void UpdateMeshFromSphericalCoordinates();
  void EmbedPointsOnSphere();
  void Optimize();
  void ComputeGradient();
  void ComputeEnergy();
  void ComputeErrors();
  void ComputeCurrentSquareDistances();
  void UpdateSphericalCoordinates(double);
  void LowPassFilter();
  void SmoothGradient();
  //set target distances
  void SmoothDisoriented();
  void AdaptiveSmooth(Array<bool>);
  void NormalizeWeights();
  void ComputeTotalWeight();

  void GetScale(vtkPolyData *, const double [3], double[3]);
  void GetCentroid(vtkPolyData *, double [3]);


};

} // namespace mirtk

#endif /* MIRTKM2SOPTIMIZER_H_ */
