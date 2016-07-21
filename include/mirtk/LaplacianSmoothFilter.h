/*
 * LaplacianSmoothFilter.h
 *
 *  Created on: 15 Nov 2013
 *
 */


#ifndef LAPLACIANSMOOTHFILTER_H_
#define LAPLACIANSMOOTHFILTER_H_

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



class LaplacianSmoothFilter{



private:
  vtkSmartPointer<vtkPolyData> _input, _output;
  int _iterations;
  double _lambda, _sigma;
  char * _pointWeightingName;
public:

  LaplacianSmoothFilter();
  ~LaplacianSmoothFilter();
  void SetIterations(int);
  void SetLambda(double);
  void SetSigma(double);
  void SetPointWeighting(char *);

  void SetInput(vtkSmartPointer<vtkPolyData>);
  void Update();
  vtkSmartPointer<vtkPolyData> GetOutput();


};

#endif /* LAPLACIANSMOOTHFILTER_H_ */
