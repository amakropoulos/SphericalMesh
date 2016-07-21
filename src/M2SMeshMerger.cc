#include "mirtk/Common.h"
#include "mirtk/M2SMeshMerger.h"

#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>

#include <vtkMath.h>

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================


M2SMeshMerger::M2SMeshMerger()
{
  _Target = vtkSmartPointer<vtkPolyData>::New();
  _SourcePts = vtkSmartPointer<vtkPolyData>::New();
  _Output = vtkSmartPointer<vtkPolyData>::New();
}


// -----------------------------------------------------------------------------


M2SMeshMerger::M2SMeshMerger(const M2SMeshMerger &other)
{
  *this = other;
}


// -----------------------------------------------------------------------------

M2SMeshMerger &M2SMeshMerger::operator =(const M2SMeshMerger &other)
{
  if (this != &other) {
    M2SMeshMerger::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

M2SMeshMerger::~M2SMeshMerger(){}

// -----------------------------------------------------------------------------

void M2SMeshMerger::CopyAttributes(const M2SMeshMerger &other){
  cout << "TODO copy attributes" << endl;
}


void M2SMeshMerger::Run(){

  // Read input point set
  vtkSmartPointer<vtkCellLocator> tree =
      vtkSmartPointer<vtkCellLocator>::New();

  tree->SetDataSet(_Target);
  tree->BuildLocator();

  vtkSmartPointer<vtkPointLocator> pointLocator =
      vtkSmartPointer<vtkPointLocator>::New();

  pointLocator->SetDataSet(_Target);
  pointLocator->BuildLocator();


  Array<int> nearestCellIds;

  vtkSmartPointer<vtkPoints> newPts =
      vtkSmartPointer<vtkPoints>::New();

  _Target->BuildLinks();
  _Target->BuildCells();
  for (int ptId=0; ptId<_SourcePts->GetNumberOfPoints(); ptId++){
    double pt[3], closestPoint[3];
    double closestDist2;
    vtkIdType cellId;  int subId;

    _SourcePts->GetPoint(ptId, pt);
    tree->FindClosestPoint(pt, closestPoint, cellId, subId, closestDist2);

    double closestDist = sqrt(closestDist2);

    int closestPointId = pointLocator->FindClosestPoint(pt);
    _Target->GetPoint(closestPointId, closestPoint);

    double cellCentriod[3] = {0.0, 0.0, 0.0};
    vtkSmartPointer<vtkIdList> cell =
        vtkSmartPointer<vtkIdList>::New();
    _Target->GetCellPoints(cellId, cell);
    for (int i=0; i<3; i++){
      int cellPtId = cell->GetId(i);
      double cellPt [3];
      _Target->GetPoint(cellPtId, cellPt);
      for (int j=0; j<3; j++){ cellCentriod[j] += cellPt[j]; }
    }
    for (int j=0; j<3; j++){ cellCentriod[j] /= 3; }
    double centriodDist = vtkMath::Distance2BetweenPoints(cellCentriod, pt);

    double threshold = 0.01;
    //if points have the same position move source point towards centriod by a small amount.
    //stop fast marching complaining
    //if (closestDist = 0){
    if (closestDist < threshold){
      //double alpha = 0.00001;
      double alpha = threshold;
      double beta = 1 - alpha;
      for (int j=0; j<3; j++){
        pt[j] = cellCentriod[j]*alpha + pt[j]*beta;
      }
    }

    newPts->InsertNextPoint(pt);
    nearestCellIds.push_back(cellId);
  }

  {
    Array<int> sortedCellIds = nearestCellIds;
    sort(sortedCellIds.begin(), sortedCellIds.end());
    auto it = adjacent_find(sortedCellIds.begin(), sortedCellIds.end());
    bool duplicates = it != sortedCellIds.end();
    if (duplicates){
      cout << "cannot merge" << endl;
      exit(1);
    }
  }

  int numCells = _Target->GetNumberOfPolys();
  Array<bool> deleteMask(numCells, false);

  vtkSmartPointer<vtkCellArray> cells =
      vtkSmartPointer<vtkCellArray>::New();


  vtkSmartPointer<vtkPoints> pts =
      vtkSmartPointer<vtkPoints>::New();

  pts->DeepCopy(_Target->GetPoints());

  for (int ptId=0; ptId<_SourcePts->GetNumberOfPoints(); ptId++){

    double pt[3];
    //_SourcePts->GetPoint(ptId, pt);
    newPts->GetPoint(ptId, pt);
    int newPtId = pts->InsertNextPoint(pt);

    vtkSmartPointer<vtkIdList> cell =
        vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> newCell =
        vtkSmartPointer<vtkIdList>::New();

    int cellId = nearestCellIds[ptId];
    _Target->GetCellPoints(cellId, cell);

    for (int i=0; i<3; i++){
      newCell->DeepCopy(cell);
      newCell->SetId(i,newPtId);
      cells->InsertNextCell(newCell);
    }

    deleteMask[cellId] = true;
  }

  _Output = vtkSmartPointer<vtkPolyData>::New();
  _Target->GetPolys()->InitTraversal();
  vtkSmartPointer<vtkIdList> cell =
      vtkSmartPointer<vtkIdList>::New();

  int cellId = 0;
  while (_Target->GetPolys()->GetNextCell(cell)){
    if (!deleteMask[cellId]){
      cells->InsertNextCell(cell);
    }
    cellId++;
  }

  _Output->Allocate();
  _Output->SetPoints(pts);
  _Output->SetPolys(cells);
}

} // end mirtk namespace
