

#include "mirtk/Common.h"
#include "mirtk/M2SGeodesicInterpolator.h"
#include "mirtk/Vtk.h"

#include <vtkFastMarchingGeodesicDistance.h>
//#include <.h>

//#include <mirtkPointSetUtils.h>

#include <vtkCellArray.h>
#include <vtkCellLocator.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================


M2SGeodesicInterpolator::M2SGeodesicInterpolator()
{
  _Target = vtkSmartPointer<vtkPolyData>::New();
  _Source = vtkSmartPointer<vtkPolyData>::New();
  _Output = vtkSmartPointer<vtkPolyData>::New();
  _ArrayName = "";
  _OutArrayName = "";
  _K = 6;
  _Norm = -1.0;
  _Lambda = 0.1;
}


// -----------------------------------------------------------------------------


M2SGeodesicInterpolator::M2SGeodesicInterpolator(const M2SGeodesicInterpolator &other)
{
  *this = other;
}

// -----------------------------------------------------------------------------

M2SGeodesicInterpolator &M2SGeodesicInterpolator::operator =(const M2SGeodesicInterpolator &other)
{
  if (this != &other) {
    M2SGeodesicInterpolator::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

M2SGeodesicInterpolator::~M2SGeodesicInterpolator(){}

// -----------------------------------------------------------------------------

void M2SGeodesicInterpolator::CopyAttributes(const M2SGeodesicInterpolator &other){
  cout << "TODO copy attributes" << endl;
}



void M2SGeodesicInterpolator::Run(){

  _Output->DeepCopy(_Target);
  int numTargetPts = _Output->GetNumberOfPoints();
  int numSourcePts = _Source->GetNumberOfPoints();

  vtkDoubleArray * sourceDataArray = (vtkDoubleArray *)
      _Source->GetPointData()->GetArray(_ArrayName.c_str());
  int numComponents = sourceDataArray->GetNumberOfComponents();

  vtkSmartPointer<vtkDoubleArray> outDataArray =
      vtkSmartPointer<vtkDoubleArray>::New();
  outDataArray->SetNumberOfComponents(numComponents);
  if (_OutArrayName == "")
    outDataArray->SetName(_ArrayName.c_str());
  else
    outDataArray->SetName(_OutArrayName.c_str());

  //set up cell locator
  vtkSmartPointer<vtkCellLocator> tree =
      vtkSmartPointer<vtkCellLocator>::New();
  tree->SetDataSet(_Source);
  tree->BuildLocator();

  //set fast marching
  string arrayName = "FMMDist";
  vtkSmartPointer< vtkFastMarchingGeodesicDistance > fmm =
      vtkSmartPointer< vtkFastMarchingGeodesicDistance >::New();
  vtkSmartPointer<vtkIdList> seed =
      vtkSmartPointer<vtkIdList>::New();


  fmm->SetFieldDataName(arrayName.c_str());
  seed->SetNumberOfIds(1);
  fmm->SetSeeds(seed);
  fmm->SetDistanceStopCriterion(-1);






  for (int targetId = 0; targetId < numTargetPts; targetId++){

    double targetPt[3], closestPoint[3];
    vtkIdType nearestCellId;  int subId;
    double closestDist2;

    _Target->GetPoint(targetId, targetPt);
    tree->FindClosestPoint(targetPt, closestPoint, nearestCellId, subId, closestDist2);
    vtkSmartPointer<vtkIdList> nearestCell = vtkSmartPointer<vtkIdList>::New();
    _Source->GetCellPoints(nearestCellId, nearestCell);

    //merge target point with source mesh
    //for input to fmm
    vtkSmartPointer<vtkPolyData> fmmMesh =
        vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> fmmPts =
        vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> fmmCells =
        vtkSmartPointer<vtkCellArray>::New();
    fmmPts->DeepCopy(_Source->GetPoints());

    //check target point does not overlap with source points in cell
    int overlapId = -1;
    for (int i=0; i<3; i++){
      int sourceId = nearestCell->GetId(i);
      double sourcePt[3];
      _Source->GetPoint(sourceId, sourcePt);

      if ( sourcePt[0] == targetPt[0] &&
           sourcePt[1] == targetPt[1] &&
           sourcePt[1] == targetPt[2] ){
        overlapId = sourceId;
        break;
      }
    }


    if (overlapId != -1){
      fmmMesh->DeepCopy(_Source);
      seed->SetId(0, overlapId);
      seed->Modified();
    }
    else{
      int newPtId = fmmPts->InsertNextPoint(targetPt);
      seed->SetId(0, newPtId);
      seed->Modified();


      //extract source cells
      vtkSmartPointer<vtkIdList> cell
          = vtkSmartPointer<vtkIdList>::New();
      _Source->GetPolys()->InitTraversal();
      int cellId = 0;
      while (_Source->GetPolys()->GetNextCell(cell)){
        if (cellId != nearestCellId){
          fmmCells->InsertNextCell(cell);
        }
        cellId++;
      }


      vtkSmartPointer<vtkIdList> newCell
          = vtkSmartPointer<vtkIdList>::New();
      for (int i=0; i<3; i++){
        newCell->DeepCopy(nearestCell);
        newCell->SetId(i, newPtId);
        fmmCells->InsertNextCell(newCell);
      }

      fmmMesh->Allocate();
      fmmMesh->SetPoints(fmmPts);
      fmmMesh->SetPolys(fmmCells);
    }



    fmm->SetInputData(fmmMesh);
    fmm->SetNumberOfVisitedPointsStopCriterion(_K + 1); // k points other than seed point!
    fmm->SetDistanceStopCriterion(-1);
    fmm->Update();


    vtkDoubleArray * distArray = (vtkDoubleArray *)
        fmm->GetOutput()->GetPointData()->GetArray(arrayName.c_str());

    Array<double> distances;
    Array<int> sourceIds;

    for (int sourcePtId=0; sourcePtId<numSourcePts; sourcePtId++){
      double d = distArray->GetTuple1(sourcePtId);
      if (d > 0){
        distances.push_back(d);
        sourceIds.push_back(sourcePtId);
      }
    }

    double meanDist = 0.0;
    for (auto d: distances){ meanDist += d; }
    meanDist /= distances.size();

    //for (auto d: distances){ cout << d << endl; }
    //cout << meanDist << endl;


    Array<double> weights;
    for (auto &d: distances){
      double w = 1 / (d + meanDist * _Lambda);
      weights.push_back(w);
    }
    //normalize weights
    double weightSum = 0.0;
    for (auto w: weights){ weightSum += w; }
    for (auto &w: weights){ w /= weightSum; }

    //for (auto w: weights){ cout << w << endl; }
    //reconstruct pointdata

    Array<double> outTup(numComponents, 0.0);

    for (int i=0; i<distances.size(); i++){
      int sourceId =  sourceIds[i];
      double w = weights[i];
      Array<double> sourceTup(numComponents, 0.0);

      sourceDataArray->GetTuple(sourceId, sourceTup.data());
      for (int j=0; j<numComponents; j++){
        outTup[j] += sourceTup[j] * w;
      }
    }

    //for (auto e: outTup){ cout << e << endl; }

    if (_Norm > 0){
      double n = 0;
      for (auto e: outTup){ n += pow(e,_Norm); }
      n = pow(n, 1.0/_Norm);
      for (auto &e: outTup){ e /= n; }
    }

    //for (auto e: outTup){ cout << e << endl; }

    outDataArray->InsertNextTuple(outTup.data());

  }
  _Output->GetPointData()->AddArray(outDataArray);
}

} // mirtk namespace

