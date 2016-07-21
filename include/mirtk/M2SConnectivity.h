/*
 * mirtkM2SOptimizer.h
 *
 *  Created on: 8 Feb 2016
 *      Author: raw11
 */

#ifndef MIRTKM2SCONNECTIVITY_H_
#define MIRTKM2SCONNECTIVITY_H_



#include "mirtk/Object.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/M2SEdge.h"

//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

namespace mirtk {

class M2SConnectivity {

public:



  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);
  mirtkPublicAttributeMacro(Array<int>, SeedPtIds);
  mirtkPublicAttributeMacro(Array<int>, OtherPtIds);
  mirtkPublicAttributeMacro(Array<Edge> , Edges);

  mirtkPublicAttributeMacro(double, ScaleFactor);
  mirtkPublicAttributeMacro(double, Threshold);
  mirtkPublicAttributeMacro(int, K);
  mirtkPublicAttributeMacro(int, MinK);

  mirtkPublicAttributeMacro(bool, Normalize);

  mirtkPublicAttributeMacro(bool, ForwardDirection);
  mirtkPublicAttributeMacro(bool, BackwardDirection);
  mirtkPublicAttributeMacro(bool, IdentityEdges);

  mirtkPublicAttributeMacro(Array<double>, DistanceComputed);


public:
  mirtkOnOffMacro(Normalize);
  mirtkOnOffMacro(ForwardDirection);
  mirtkOnOffMacro(BackwardDirection);
  mirtkOnOffMacro(IdentityEdges);

  // Construction/Destruction
public:

   /// Constructor
  M2SConnectivity();

   /// Copy constructor
  M2SConnectivity(const M2SConnectivity &);

   /// Assignment operator
  M2SConnectivity &operator =(const M2SConnectivity &);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const M2SConnectivity &);

  /// Destructor
  virtual ~M2SConnectivity();

  void Update();
  void NormalizeEdgeLengths();
  void Initialize();

  void ComputeGeodesics();
  void InitializeEdges();
  void InitializeEdgeDistances();
  void ScaleEdgeDistances();

};

} // namespace mirtk

#endif /* MIRTKM2SCONNECTIVITY_H_ */
