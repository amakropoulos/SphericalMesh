
#ifndef M2SDIFFUSER_H
#define M2SDIFFUSER_H

#include "mirtk/Object.h" // for attribute macros
#include "mirtk/EdgeTable.h"

//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

namespace mirtk {

class M2SDiffuser {



	// ---------------------------------------------------------------------------
	// Attributes

	mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);

	//TODO Allow a "target" mesh to be used to compute edge lengths which can
	// be used to weight the diffusion.
  //mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Target);

	mirtkPublicAttributeMacro(Array<bool>, SourceMask);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);
  mirtkPublicAttributeMacro(int, MaxIterations);
  mirtkPublicAttributeMacro(bool, AdaptiveTermination);

  mirtkPublicAttributeMacro(bool, Debug);


  // TODO compute average movement of points and terminate when
  // not greater than some threshold
  //mirtkPublicAttributeMacro(double, ConvergenceThreshold);

	 // Construction/Destruction
	public:

	  /// Constructor
	  M2SDiffuser();

	  /// Copy constructor
	  M2SDiffuser(const M2SDiffuser &);

	  /// Assignment operator
	  M2SDiffuser &operator =(const M2SDiffuser &);

	  /// Copy attributes of this class from another instance
	  void CopyAttributes(const M2SDiffuser &);

	  /// Destructor
	  virtual ~M2SDiffuser();

	  //Methods

	  void Run();


}; // end class

} // namespace mirtk

#endif // M2SDIFFUSER_H
