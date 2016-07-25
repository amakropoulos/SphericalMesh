/*
 * Remesher.h
 *
 *  Created on: 13 Feb 2015
 *      Author: rob
 */

#ifndef MIRTKMESHTOSPHEREREMESHER_H_
#define MIRTKMESHTOSPHEREREMESHER_H_

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPriorityQueue.h>

#include "mirtk/Array.h"
#include "mirtk/Object.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PriorityQueue.h"


namespace mirtk {

class Remesher{


private:
  struct RemeshingOperation {
    double priority = 0.0;
    bool collapse = false;
    int ptId1 = -1;
    int ptId2 = -1;

    bool operator > (const RemeshingOperation& ro) const
    {
        return (priority > ro.priority);
    }

    bool operator < (const RemeshingOperation& ro) const
    {
        return (priority < ro.priority);
    }


  };


  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);
  mirtkPublicAttributeMacro(EdgeTable, Edges);

  mirtkPublicAttributeMacro(int, MaxIterations);
  mirtkPublicAttributeMacro(double, MaxEdgeLength);
  mirtkPublicAttributeMacro(double, MinEdgeLength);
  mirtkPublicAttributeMacro(Array<double>, EdgeLengths);
  mirtkPublicAttributeMacro(Array<RemeshingOperation>, RemeshingOperationQueue);

  mirtkPublicAttributeMacro(Array<bool>, MarkedPts);
  mirtkPublicAttributeMacro(Array<bool>, MarkedCells);

  mirtkPublicAttributeMacro(Array<bool>, DeletedPts);
  mirtkPublicAttributeMacro(Array<bool>, DeletedCells);

  mirtkPublicAttributeMacro(Array<int>, NumAdjacentPts);

  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPoints>, Points);
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkCellArray>, Cells);

  mirtkPublicAttributeMacro(bool, Debug);
  mirtkPublicAttributeMacro(bool, Verbose);





  private:

    void Midpoint(double [],double [],double []);
    void Initialize();
    void ComputeEdgeLengths();
    void QueueEdgesByLength();
    int ProcessEdgeQueue();

    void GetSharedNeighbours(int, int);
    vtkSmartPointer<vtkIdList> GetSharedPoints(int, int, int, int);


    bool CollapseEdge(int, int);
    bool CollapseTriad(int);
    bool SplitEdge(int, int);
    void MarkPtIds(vtkIdList *);

    void UpdateCells(Array<int> &);
    Array<int> UpdatePoints();

    void Finalize();
    vtkSmartPointer<vtkIdList> GetSharedNeighbours(
        const int *, int, const int *,  int);

    vtkSmartPointer<vtkIdList> GetSharedCells(int, int);

    bool CheckTopology();

  public:

    /// Constructor
    Remesher();

    /// Copy constructor
    Remesher(const Remesher &);

    /// Assignment operator
    Remesher &operator =(const Remesher &);

    /// Copy attributes of this class from another instance
    void CopyAttributes(const Remesher &);

    /// Destructor
    virtual ~Remesher();

    void Run();




};

} //end mirtk namespace

#endif /* MIRTKMESHTOSPHEREREMESHER_H_ */
