
#ifndef M2SMESHMERGER_H
#define M2SMESHMERGER_H

#include "mirtk/Object.h"
#include "mirtk/EdgeTable.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

namespace mirtk {

class M2SMeshMerger{


	// ---------------------------------------------------------------------------
	// Attributes

	mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Target);
	mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, SourcePts);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);
	 // Construction/Destruction
	public:

	  /// Constructor
	  M2SMeshMerger();

	  /// Copy constructor
	  M2SMeshMerger(const M2SMeshMerger &);

	  /// Assignment operator
	  M2SMeshMerger &operator =(const M2SMeshMerger &);

	  /// Copy attributes of this class from another instance
	  void CopyAttributes(const M2SMeshMerger &);

	  /// Destructor
	  virtual ~M2SMeshMerger();

	  //Methods
	  void Run();
}; // end class

} // namespace mirtk

#endif // M2SMESHMERGER_H
