/*
 * Remesher.h
 *
 *  Created on: 13 Feb 2015
 *      Author: rob
 */

#ifndef REMESHER_H_
#define REMESHER_H_

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPriorityQueue.h>

class Remesher{

  private:
    int _maxIterations;
    double _maxEdgeLength,_minEdgeLength;

    double GetEdgeLength(double [],double []);

    void Midpoint(double [],double [],double []);

    vtkSmartPointer<vtkPolyData> CleanPolyData(vtkPolyData *);

    vtkSmartPointer<vtkPolyData> ProcessInput(vtkPolyData *);

    vtkSmartPointer<vtkPolyData> BuildEdgeStructure(vtkPolyData *);

    vtkSmartPointer<vtkDoubleArray> ComputeEdgeLengths(vtkPolyData *);

    vtkSmartPointer<vtkPolyData> RemoveDegenerateCells(vtkPolyData * polydata);

    vtkSmartPointer<vtkPriorityQueue> QueueEdgesByLength(vtkPolyData *, double,
        double);

    vtkSmartPointer<vtkIdList> GetAdjacentCellIds(int,int,vtkPolyData *);

    void DeleteEdgesAfterCollapse(int, int, vtkPolyData *);

    void CreateNewCells(int, int, int, int, int, vtkCellArray *,vtkPolyData *);

    void CombineCells(vtkPolyData *,vtkCellArray *,vtkIntArray *);

    void AveragePointdata(int, int, int, vtkPolyData *);

    vtkSmartPointer<vtkPolyData> ProcessEdgeQueue(vtkPriorityQueue *,
        vtkPolyData *, vtkPolyData *);



  public:
    Remesher();
    ~Remesher();

    void SetMaxEdgeLength(double);
    void SetMinEdgeLength(double);
    void SetMaxIterations(int);
    vtkSmartPointer<vtkPolyData> Remesh(vtkPolyData * );

};

#endif /* REMESHER_H_ */
