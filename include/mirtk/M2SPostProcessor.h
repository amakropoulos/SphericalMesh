
#ifndef M2SPOSTPROCESSOR_H
#define M2SPOSTPROCESSOR_H

#include "mirtk/Object.h" // for attribute macros
#include "mirtk/EdgeTable.h"

//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


namespace mirtk {

class M2SPostProcessor {



	// ---------------------------------------------------------------------------
	// Attributes

	mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);

  //Which points are already smooth
  mirtkPublicAttributeMacro(Array<bool>, DiscontinuityMask);

	mirtkPublicAttributeMacro(double, Threshold);
  mirtkPublicAttributeMacro(double, Percentile);

  mirtkPublicAttributeMacro(int, MaxIterations);
  mirtkPublicAttributeMacro(int, MaxDiffusionIterations);

  mirtkPublicAttributeMacro(bool, Debug);

	 // Construction/Destruction
	public:

	  /// Constructor
	  M2SPostProcessor();

	  /// Copy constructor
	  M2SPostProcessor(const M2SPostProcessor &);

	  /// Assignment operator
	  M2SPostProcessor &operator =(const M2SPostProcessor &);

	  /// Copy attributes of this class from another instance
	  void CopyAttributes(const M2SPostProcessor &);

	  /// Destructor
	  virtual ~M2SPostProcessor();


	  //Methods
	  void DetectDiscontinuities();
	  void SmoothDiscontinuities();
	  int GetNumberOfDiscontinuities();
	  void Run();

}; // end class

} // namespace mirtk

#endif // M2SPOSTPROCESSOR_H
