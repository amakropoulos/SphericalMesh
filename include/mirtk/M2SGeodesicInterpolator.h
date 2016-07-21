/*
 * mirtkM2SOptimizer.h
 *
 *  Created on: 8 Feb 2016
 *      Author: raw11
 */

#ifndef MIRTKM2SGEODESICINTERPOLATOR_H_
#define MIRTKM2SGEODESICINTERPOLATOR_H_

#include "mirtk/Array.h"
#include "mirtk/Object.h"
#include "mirtk/EdgeTable.h"

//my includes

//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

namespace mirtk {

class M2SGeodesicInterpolator {

  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Target);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Source);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);
  mirtkPublicAttributeMacro(string, ArrayName);
  mirtkPublicAttributeMacro(string, OutArrayName);
  mirtkPublicAttributeMacro(int, K);
  mirtkPublicAttributeMacro(double, Norm);
  mirtkPublicAttributeMacro(double, Lambda);
  mirtkPublicAttributeMacro(int, Iterations);

  // Construction/Destruction
 public:

   /// Constructor
   M2SGeodesicInterpolator();

   /// Copy constructor
   M2SGeodesicInterpolator(const M2SGeodesicInterpolator &);

   /// Assignment operator
   M2SGeodesicInterpolator &operator =(const M2SGeodesicInterpolator &);

   /// Copy attributes of this class from another instance
   void CopyAttributes(const M2SGeodesicInterpolator &);

   /// Destructor
   virtual ~M2SGeodesicInterpolator();

  void Initialize();
  void Run();
  void GetOutput();

};

} // namespace mirtk

#endif /* MIRTKM2SGEODESICINTERPOLATOR_H_ */
