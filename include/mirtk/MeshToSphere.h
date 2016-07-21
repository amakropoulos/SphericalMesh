
#ifndef MESHTOSPHERE_H
#define MESHTOSPHERE_H

// #include "mirtk/Object.h" // for attribute macros
// #include "mirtk/EdgeTable.h"
#include "mirtk/M2SParameters.h"

//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>



namespace mirtk {

class MeshToSphere {


  public:
    struct Params{
      //int targetNumPoints = -1; //inclination
      double minEdgeLength = -1; //inclination
      double maxEdgeLength = -1; //inclination
      double connectivityDistanceThreshold = -1;
      double epsilon = 0.0000001;
      int maxIterations = 100;
    };


  private:
    // ---------------------------------------------------------------------------
    // Types
    struct SphericalVector
    {
      double theta = 0; //inclination
      double phi = 0; //azimuthal angle
    };






	// ---------------------------------------------------------------------------
	// Attributes

  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);
  mirtkPublicAttributeMacro(EdgeTable, Edges);
  mirtkPublicAttributeMacro(Array<double>, TargetDistances);
  mirtkPublicAttributeMacro(Array<double>, CurrentSquareDistances);
  mirtkPublicAttributeMacro(Array<double>, Errors);

  mirtkPublicAttributeMacro(Array<LevelParameters>, Parameters);

  mirtkPublicAttributeMacro(Array< vtkSmartPointer<vtkPolyData> >, MeshHierarchy);
  mirtkPublicAttributeMacro(Array< Array<int> >, PtIdHierarchy);
  mirtkPublicAttributeMacro(Array< Array<int> >, fullMeshPtIdMapping);
  mirtkPublicAttributeMacro(Array<int>, Level);


  mirtkPublicAttributeMacro(Array<bool>, ActivePoints);
  mirtkPublicAttributeMacro(Array<bool>, StationaryPoints);

  mirtkPublicAttributeMacro(GenericSparseMatrix<double>, Geodesics);


  mirtkPublicAttributeMacro(Array<SphericalVector>, OriginalSphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, OldSphericalCoordinates);
	mirtkPublicAttributeMacro(Array<SphericalVector>, SphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, BestSphericalCoordinates);
  mirtkPublicAttributeMacro(Array<SphericalVector>, Gradient);
  mirtkPublicAttributeMacro(Array<SphericalVector>, StepSizeW);
  mirtkPublicAttributeMacro(double, Energy);
  mirtkPublicAttributeMacro(int, CurrentLevel);
  mirtkPublicAttributeMacro(bool, FinalInterpolation);

  //post processing settings
  mirtkPublicAttributeMacro(int, PostProcessingIterations);
  mirtkPublicAttributeMacro(double, PostProcessingPercentile);
  mirtkPublicAttributeMacro(double, PostProcessingThreshold);

  mirtkPublicAttributeMacro(double, LocalWeighting);
  mirtkPublicAttributeMacro(double, GlobalWeighting);
  mirtkPublicAttributeMacro(bool, Debug);
  mirtkPublicAttributeMacro(int, NumberOfInterpolationPoints);
  mirtkPublicAttributeMacro(double, FinalInterpolationThreshold);

  public:
    mirtkOnOffMacro(FinalInterpolation);



	 // Construction/Destruction
	public:

	  /// Constructor
    MeshToSphere();

	  /// Copy constructor
	  MeshToSphere(const MeshToSphere &);

	  /// Assignment operator
	  MeshToSphere &operator =(const MeshToSphere &);

	  /// Copy attributes of this class from another instance
	  void CopyAttributes(const MeshToSphere &);

	  /// Destructor
	  virtual ~MeshToSphere();



	  //Methods
	  void Initialize();

	  void Run();

	  void ComputeTargetDistances();

	  void ComputeCurrentSquareDistances();

	  void UpdateSphericalCoordinatesFromMesh();
    void UpdateSphericalCoordinates(double);
    void UpdateMeshFromSphericalCoordinates();

    void ComputeMeshHierarchy();
    void InterpolateAngles(int, int);
    void EmbedPointsOnSphere(vtkPolyData *);
    void LowPassFilter(vtkPolyData *);
    void InitialiseTargetDistances();
    void InitialiseSphericalMesh();
    void Remesh(double, double, bool);


}; // end class

} // namespace mirtk

#endif // MESHTOSPHERE_H
