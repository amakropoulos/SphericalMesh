/*
 * SphereQuadrisect.cc
 *
 *  Created on: 15 Nov 2013
 *
 */


#ifndef SPHEREQUADRISECT_H_
#define SPHEREQUADRISECT_H_

#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>



class SphereQuadrisect{



private:

public:

  SphereQuadrisect();
  ~SphereQuadrisect();
  vtkSmartPointer<vtkPolyData> Quadrisect(vtkSmartPointer<vtkPolyData>,int);

};

#endif /* SPHEREQUADRISECT_H_ */
