/*
 * Remesher.cc
 *
 *  Created on: 13 Feb 2015
 *      Author: rob
 */

#include "mirtk/Common.h"
#include "mirtk/MeshToSphereRemesher2.h"

#include "mirtk/PointSetUtils.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/UnorderedMap.h"
#include "mirtk/OrderedMap.h"
#include "mirtk/Pair.h"

#include <float.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkExtractEdges.h>
#include <vtkPriorityQueue.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>
#include "vtkNew.h"


namespace mirtk {

Remesher::Remesher(){

  _Input = vtkSmartPointer<vtkPolyData>::New();
  _Output = vtkSmartPointer<vtkPolyData>::New();
  _Points = vtkSmartPointer<vtkPoints>::New();
  _Cells = vtkSmartPointer<vtkCellArray>::New();


  _MaxEdgeLength = DBL_MAX;
  _MinEdgeLength = -1;
  _MaxIterations = 1;
  _Verbose = false;
  _Debug = false;
}


Remesher::~Remesher(){
}


void Remesher::Initialize(){

  vtkNew<vtkTriangleFilter> tf;
  tf->SetInputData(_Output);
  tf->Update();
  _Output = tf->GetOutput();
  _Output->BuildCells();
  _Output->BuildLinks();

  _MarkedPts.clear();
  _DeletedPts.clear();
  int numPts = _Output->GetNumberOfPoints();
  for (int ptId=0; ptId<numPts; ptId++){
    _MarkedPts.push_back(false);
    _DeletedPts.push_back(false);
  }

  _MarkedCells.clear();
  _DeletedCells.clear();
  int numCells = _Output->GetNumberOfPolys();
  for (int cellId=0; cellId<numCells; cellId++){
    _MarkedCells.push_back(false);
    _DeletedCells.push_back(false);
  }

  _Edges.Initialize(_Output);



  _NumAdjacentPts.clear();
  for (int ptId=0; ptId<numPts; ptId++){
    int adjs = _Edges.NumberOfAdjacentPoints(ptId);
    _NumAdjacentPts.push_back(adjs);
  }

  _Points->DeepCopy(_Output->GetPoints());
  _Cells->DeepCopy(_Output->GetPolys());
}

void midpoint(double p1 [],double p2 [],double mp []){
  for (int i=0;i<3;i++) mp[i]=(p1[i]+p2[i])/2;
}



double distance(double pt1 [3], double pt2 [3]){

  double dist = 0;
  for (int i=0; i<3; i++){
    dist += pow(pt2[i] - pt1[i], 2);
  }

  dist = sqrt(dist);

  return dist;
}


void Remesher::ComputeEdgeLengths(){
  EdgeIterator ei(_Edges);
  ei.InitTraversal();
  _EdgeLengths.clear();

  int pid1, pid2;
  int edge_id=0;
  while ( ei.GetNextEdge(pid1, pid2) != -1 ){
    double pt1[3], pt2[3];

    _Output->GetPoint(pid1, pt1);
    _Output->GetPoint(pid2, pt2);

    double d = distance(pt1, pt2);
    _EdgeLengths.push_back(d);
  }
}

 void Remesher::QueueEdgesByLength(){

  int splitNum = 0, collapseNum = 0;

  EdgeIterator ei(_Edges);
  ei.InitTraversal();
  _RemeshingOperationQueue.clear();

  int ptId1, ptId2;
  int edgeId = 0;
  while ( ei.GetNextEdge(ptId1, ptId2) != -1 ){

    double d = _EdgeLengths[edgeId];
    double priority;

    if (d > _MaxEdgeLength){
      RemeshingOperation ro;
      ro.collapse = false;
      ro.priority = _MaxEdgeLength / (d-_MaxEdgeLength);
      ro.ptId1 = ptId1;
      ro.ptId2 = ptId2;
      _RemeshingOperationQueue.push_back(ro);
      splitNum+=1;
    }
    else if (d < _MinEdgeLength){
      RemeshingOperation ro;
      ro.collapse = true;
      ro.priority = _MinEdgeLength / (_MinEdgeLength-d);
      ro.ptId1 = ptId1;
      ro.ptId2 = ptId2;
      _RemeshingOperationQueue.push_back(ro);

      collapseNum += 1;
    }
    edgeId++;
  }

  sort(_RemeshingOperationQueue.begin(), _RemeshingOperationQueue.end());
  if (_Verbose){
    cout << "Attempting to split " << splitNum << " edges and collapse " << collapseNum
      << " edges (" << _Edges.NumberOfEdges() << " total)." << endl;
  }
}


int Remesher::ProcessEdgeQueue(){

  int collapseCount = 0;
  int splitCount = 0;

  int edgeId;

  for (int i=0; i<_RemeshingOperationQueue.size(); i++){

    RemeshingOperation ro = _RemeshingOperationQueue[i];

    if (_MarkedPts[ro.ptId1] || _MarkedPts[ro.ptId2])
      continue;

    if (ro.collapse){
      bool success = CollapseEdge(ro.ptId1, ro.ptId2);
      collapseCount += success;
    }
    else{
      bool success = SplitEdge(ro.ptId1, ro.ptId2);
      splitCount += success;
    }

  }


  if (_Verbose){
    cout << splitCount << " edges split and " << collapseCount
      << " edges collapsed." << endl;
  }


  return collapseCount+splitCount;


 }



vtkSmartPointer<vtkIdList> Remesher::GetSharedCells(int ptId1, int ptId2){

  //get adajacencies for both and

  vtkSmartPointer<vtkIdList> cellIds1, cellIds2;
  cellIds1 = vtkSmartPointer<vtkIdList>::New();
  cellIds2 = vtkSmartPointer<vtkIdList>::New();

  _Output->GetPointCells(ptId1, cellIds1);
  _Output->GetPointCells(ptId2, cellIds2);
  cellIds1->IntersectWith(cellIds2); // common cells that both points are in

  return cellIds1;
}


vtkSmartPointer<vtkIdList> idListUnion(vtkIdList * a, vtkIdList * b){

  vtkSmartPointer<vtkIdList> c =
      vtkSmartPointer<vtkIdList>::New();

  c->DeepCopy(a);

  for (int i=0; i<b->GetNumberOfIds(); i++){
    int id = b->GetId(i);
    c->InsertUniqueId(id);
  }

  return c;
}


vtkSmartPointer<vtkIdList> Remesher::GetSharedPoints(int pId0,int pId1,
    int cellId0,
    int cellId1){

  int i,otherId;


  vtkSmartPointer<vtkIdList> cellPIds0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellPIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherPIds = vtkSmartPointer<vtkIdList>::New();
  otherPIds->Allocate(2);

  _Output->GetCellPoints(cellId0, cellPIds0);
  _Output->GetCellPoints(cellId1, cellPIds1);

  for (i=0;i<3;i++){
    otherId = cellPIds0->GetId(i);
    if (otherId != pId0 && otherId != pId1){
      otherPIds->InsertNextId(otherId);
      break;
    }
  }

  for (i=0;i<3;i++){
    otherId = cellPIds1->GetId(i);
    if (otherId != pId0 && otherId != pId1){
      otherPIds->InsertNextId(otherId);
      break;
    }
  }

  return otherPIds;
}



bool Remesher::CollapseTriad(int centerPtId){

  vtkSmartPointer<vtkIdList> cellIds =
      vtkSmartPointer<vtkIdList>::New();
  _Output->GetPointCells(centerPtId, cellIds);

  int numCellIds = cellIds->GetNumberOfIds();
  if (numCellIds != 3){
    cout << "WARNING - tried to collapse non-triad" << endl;
    return false;
  }


  //test if any of the cells have been marked
  //don't collapse edge if any cells have been marked
  for (int i=0; i<cellIds->GetNumberOfIds(); i++){
    int cellId = cellIds->GetId(i);
    if (_MarkedCells[cellId]) { return false; }
  }

  vtkSmartPointer<vtkIdList> cell, newCell;
  cell = vtkSmartPointer<vtkIdList>::New();
  newCell = vtkSmartPointer<vtkIdList>::New();

  _Output->GetCellPoints(cellIds->GetId(0), cell);
  newCell->DeepCopy(cell);

  /*
  cout << "centre " << centerPtId << endl;
  for (int i=0; i<newCell->GetNumberOfIds(); i++){
    int ptId = newCell->GetId(i);
    cout << ptId << " ";
  }
  cout << endl;
*/

  Array<int> outerPtIds;
  for (int i=0; i<cell->GetNumberOfIds(); i++){
    int ptId = cell->GetId(i);
    if (ptId != centerPtId){ outerPtIds.push_back(ptId); }
  }

  _Output->GetCellPoints(cellIds->GetId(1), cell);
  for (int i=0; i<cell->GetNumberOfIds(); i++){
    int ptId = cell->GetId(i);
    if (ptId != centerPtId && ptId != outerPtIds[0] && ptId != outerPtIds[1]){
      outerPtIds.push_back(ptId);
      break;
    }
  }

  for (int i=0; i<newCell->GetNumberOfIds(); i++){
    int ptId = newCell->GetId(i);
    if (ptId == centerPtId){
      newCell->SetId(i, outerPtIds[2]);
      break;
    }
  }
/*
  for (int i=0; i<newCell->GetNumberOfIds(); i++){
    int ptId = newCell->GetId(i);
    cout << ptId << " ";
  }
  cout << endl;
*/

  _Cells->InsertNextCell(newCell);

  for (int i=0; i<numCellIds; i++){
    int cellId = cellIds->GetId(i);
    _DeletedCells[cellId] = true;
    _MarkedCells[cellId] = true;
  }

  _DeletedPts[centerPtId] = true;
  _MarkedPts[centerPtId] = true;

  for (int i=0; i<outerPtIds.size(); i++){
    _MarkedPts[outerPtIds[i]] = true;
  }
  return true;
}


void swap(int &a, int &b){
  int tmp = a;
  a = b;
  b = tmp;
}


void Remesher::MarkPtIds(vtkIdList * cellIdList){

  for (int i=0; i<cellIdList->GetNumberOfIds(); i++){
    int cellId = cellIdList->GetId(i);
    vtkSmartPointer<vtkIdList> cellPtIds =
        vtkSmartPointer<vtkIdList>::New();
   _Output->GetCellPoints(cellId, cellPtIds);

   for (int j=0; j<cellPtIds->GetNumberOfIds(); j++){
     int ptId = cellPtIds->GetId(j);
     _MarkedPts[ptId] = true;
   }

  }

}

bool Remesher::CollapseEdge(int ptId1, int ptId2){

  if (_NumAdjacentPts[ptId1] == 3){
    return CollapseTriad(ptId1);
  }

  if (_NumAdjacentPts[ptId2] == 3){
    return CollapseTriad(ptId2);
  }

  vtkSmartPointer<vtkIdList> cellIds1, cellIds2;
  cellIds1 = vtkSmartPointer<vtkIdList>::New();
  cellIds2 = vtkSmartPointer<vtkIdList>::New();
  _Output->GetPointCells(ptId1, cellIds1);
  _Output->GetPointCells(ptId2, cellIds2);


  //test if any of the cells have been marked
  //don't collapse edge if any cells have been marked
  for (int i=0; i<cellIds1->GetNumberOfIds(); i++){
    int cellId = cellIds1->GetId(i);
    if (_MarkedCells[cellId]) { return false; }
  }
  for (int i=0; i<cellIds2->GetNumberOfIds(); i++){
    int cellId = cellIds2->GetId(i);
    if (_MarkedCells[cellId]) { return false; }
  }


  vtkSmartPointer<vtkIdList> sharedCellIds =
      vtkSmartPointer<vtkIdList>::New();
  sharedCellIds->DeepCopy(cellIds1);
  sharedCellIds->IntersectWith(cellIds2);


  if (sharedCellIds->GetNumberOfIds() != 2){
      cout << "WARNING - mesh topology error!" << endl;
      return false;
  }

  int sharedCell1Id = sharedCellIds->GetId(0);
  int sharedCell2Id = sharedCellIds->GetId(1);


  vtkSmartPointer<vtkIdList> unsharedCellIds =
      vtkSmartPointer<vtkIdList>::New();
  for (int i=0; i<cellIds1->GetNumberOfIds(); i++){
    int cellId = cellIds1->GetId(i);
    if (cellId != sharedCell1Id && cellId != sharedCell2Id)
      unsharedCellIds->InsertUniqueId(cellId);
  }
  for (int i=0; i<cellIds2->GetNumberOfIds(); i++){
    int cellId = cellIds2->GetId(i);
    if (cellId != sharedCell1Id && cellId != sharedCell2Id)
      unsharedCellIds->InsertUniqueId(cellId);
  }




  vtkSmartPointer<vtkIdList> sharedNeighbours =
      GetSharedPoints(ptId1, ptId2, sharedCell1Id, sharedCell2Id);

  if (sharedNeighbours->GetNumberOfIds() != 2)
      return false;

  int sharedNeighbourId1 = sharedNeighbours->GetId(0);
  int sharedNeighbourId2 = sharedNeighbours->GetId(1);


  //TODO handle triad in shared cell
  //if numbering point has 3 adjacencies delete all its cells
  if (_NumAdjacentPts[sharedNeighbourId1] == 3 ||
      _NumAdjacentPts[sharedNeighbourId2] == 3)
      return false;



  double pt1[3], pt2[3], m[3];


  _Output->GetPoint(ptId1, pt1);
  _Output->GetPoint(ptId2, pt2);

  midpoint(pt1,pt2,m);

  /*
  cout << pt1[0] << "," << pt1[1] << "," << pt1[2] << endl;
  cout << pt2[0] << "," << pt2[1] << "," << pt2[2] << endl;
  cout << m[0] << "," << m[1] << "," << m[2] << endl;
*/

  int newPtId = _Points->GetNumberOfPoints();



  //cout << "newPtId " << newPtId << ", total " << _Output->GetNumberOfPoints() << endl;



  //generate candidate cells
  Array< vtkSmartPointer<vtkIdList> > candidateCells;

  for (int i=0; i<unsharedCellIds->GetNumberOfIds(); i++){
    int cellId = unsharedCellIds->GetId(i);
    vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
    _Output->GetCellPoints(cellId, cell);

    for (int j=0; j<cell->GetNumberOfIds(); j++){
      int ptId = cell->GetId(j);
      if (ptId == ptId1 || ptId == ptId2){
        cell->SetId(j, newPtId);
      }
    }
    candidateCells.push_back(cell);
  }
 //std::pair<int, int>
  OrderedMap< Pair<int, int> , int > edgeCount;

  for (int i=0; i<candidateCells.size(); i++){
    vtkIdList * cell = candidateCells[i];
    int id0 = cell->GetId(0);
    int id1 = cell->GetId(1);
    int id2 = cell->GetId(2);

    int idTmp;

    //sort ids
    if (id2 < id1){ swap(id2, id1); }
    if (id2 < id0){ swap(id2, id0); }
    if (id1 < id0){ swap(id1, id0); }

    //first id is always less than second
    Pair<int,int> edge0 = MakePair(id0, id1);
    Pair<int,int> edge1 = MakePair(id0, id2);
    Pair<int,int> edge2 = MakePair(id1, id2);

    if ( edgeCount.find(edge0) == edgeCount.end() ){
      edgeCount[edge0] = 1;
    }
    else{
      edgeCount[edge0]++;
    }

    if ( edgeCount.find(edge1) == edgeCount.end() ){
      edgeCount[edge1] = 1;
    }
    else{
      edgeCount[edge1]++;
    }

    if ( edgeCount.find(edge2) == edgeCount.end() ){
      edgeCount[edge2] = 1;
    }
    else{
      edgeCount[edge2]++;
    }

  }

  //check if each half-edge appears more than twice
  for(auto iter: edgeCount){
    int count = iter.second;

    if (count > 2){
      return false; //edge collapse changes topology!
    }
  }



  //add new cells to array and mark old cells for deletion
  for (int i=0; i<candidateCells.size(); i++){
    _Cells->InsertNextCell(candidateCells[i]);
  }

  //mark old cells for deletion
  for (int i=0; i<unsharedCellIds->GetNumberOfIds(); i++){
    int cellId = unsharedCellIds->GetId(i);
    _DeletedCells[cellId] = true;
    _MarkedCells[cellId] = true;
  }

  //mark shared cells for deletion
  _DeletedCells[sharedCell1Id] = true;
  _DeletedCells[sharedCell2Id] = true;
  _MarkedCells[sharedCell1Id] = true;
  _MarkedCells[sharedCell2Id] = true;


  //Insert new point at midpoint
  _Points->InsertNextPoint(m);

  //Mark collapsed points for deletion
  _DeletedPts[ptId1] = true;
  _DeletedPts[ptId2] = true;



  MarkPtIds(unsharedCellIds);
  MarkPtIds(sharedCellIds);

  _MarkedPts[ptId1] = true;
  _MarkedPts[ptId2] = true;

  //_NumAdjacentPts[ptId1] = 0;
  //_NumAdjacentPts[ptId2] = 0;
  return true;
}





bool Remesher::SplitEdge(int ptId1, int ptId2){

  vtkSmartPointer<vtkIdList> cellIds1, cellIds2;
   cellIds1 = vtkSmartPointer<vtkIdList>::New();
   cellIds2 = vtkSmartPointer<vtkIdList>::New();
   _Output->GetPointCells(ptId1, cellIds1);
   _Output->GetPointCells(ptId2, cellIds2);


   //test if any of the cells have been marked
   //don't collapse edge if any cells have been marked
   for (int i=0; i<cellIds1->GetNumberOfIds(); i++){
     int cellId = cellIds1->GetId(i);
     if (_MarkedCells[cellId]) { return false; }
   }
   for (int i=0; i<cellIds2->GetNumberOfIds(); i++){
     int cellId = cellIds2->GetId(i);
     if (_MarkedCells[cellId]) { return false; }
   }


   vtkSmartPointer<vtkIdList> sharedCellIds =
       vtkSmartPointer<vtkIdList>::New();
   sharedCellIds->DeepCopy(cellIds1);
   sharedCellIds->IntersectWith(cellIds2);


   if (sharedCellIds->GetNumberOfIds() != 2){
       cout << sharedCellIds->GetNumberOfIds();
       cout << "WARNING - possible mesh topology error!" << endl;
       return false;
   }

   int sharedCell1Id = sharedCellIds->GetId(0);
   int sharedCell2Id = sharedCellIds->GetId(1);


   vtkSmartPointer<vtkIdList> unsharedCellIds =
       vtkSmartPointer<vtkIdList>::New();
   for (int i=0; i<cellIds1->GetNumberOfIds(); i++){
     int cellId = cellIds1->GetId(i);
     if (cellId != sharedCell1Id && cellId != sharedCell2Id)
       unsharedCellIds->InsertUniqueId(cellId);
   }
   for (int i=0; i<cellIds2->GetNumberOfIds(); i++){
     int cellId = cellIds2->GetId(i);
     if (cellId != sharedCell1Id && cellId != sharedCell2Id)
       unsharedCellIds->InsertUniqueId(cellId);
   }

   vtkSmartPointer<vtkIdList> sharedNeighbours =
       GetSharedPoints(ptId1, ptId2, sharedCell1Id, sharedCell2Id);

   if (sharedNeighbours->GetNumberOfIds() != 2)
       return false;

   int sharedNeighbourId1 = sharedNeighbours->GetId(0);
   int sharedNeighbourId2 = sharedNeighbours->GetId(1);


   double pt1[3], pt2[3], m[3];

   _Output->GetPoint(ptId1, pt1);
   _Output->GetPoint(ptId2, pt2);

   midpoint(pt1,pt2,m);

   int newPtId = _Points->GetNumberOfPoints();


   //generate new cells

   for (int i=0; i<2; i++){
     vtkSmartPointer<vtkIdList> sharedCell, newCell;
     sharedCell = vtkSmartPointer<vtkIdList>::New();
     newCell = vtkSmartPointer<vtkIdList>::New();

     int sharedCellId = sharedCellIds->GetId(i);
     _Output->GetCellPoints(sharedCellId, sharedCell);

     newCell->DeepCopy(sharedCell);

     for (int j=0; j<3; j++){
       int ptId = newCell->GetId(j);
       if (ptId == ptId1){ newCell->SetId(j, newPtId); }
     }
     _Cells->InsertNextCell(newCell);

     newCell->DeepCopy(sharedCell);
     for (int j=0; j<3; j++){
       int ptId = newCell->GetId(j);
       if (ptId == ptId2){ newCell->SetId(j, newPtId); }
     }
     _Cells->InsertNextCell(newCell);
   }

      //mark unshared cells
   for (int i=0; i<unsharedCellIds->GetNumberOfIds(); i++){
     int cellId = unsharedCellIds->GetId(i);
     _MarkedCells[cellId] = true;
   }

   //mark shared cells for deletion
   _DeletedCells[sharedCell1Id] = true;
   _DeletedCells[sharedCell2Id] = true;
   _MarkedCells[sharedCell1Id] = true;
   _MarkedCells[sharedCell2Id] = true;





   //Insert new point at midpoint
   _Points->InsertNextPoint(m);

   //Mark edge points
   _MarkedPts[ptId1] = true;
   _MarkedPts[ptId2] = true;

   return true;
}

Array<int> Remesher::UpdatePoints(){

  Array<int> mapping;

  vtkSmartPointer<vtkPoints> newPoints =
      vtkSmartPointer<vtkPoints>::New();

  int oldNumPoints = _Output->GetNumberOfPoints();
  for (int oldPtId=0; oldPtId<_Points->GetNumberOfPoints(); oldPtId++){
    if (oldPtId < oldNumPoints && _DeletedPts[oldPtId]){
      mapping.push_back(-1);
    }
    else{
      double pt [3];
      _Points->GetPoint(oldPtId, pt);
      int newPtId = newPoints->InsertNextPoint(pt);
      mapping.push_back(newPtId);
    }
  }

  _Points->DeepCopy(newPoints);

  return mapping;
}


void Remesher::UpdateCells(Array<int> &mapping){

  vtkSmartPointer<vtkCellArray> newCells =
      vtkSmartPointer<vtkCellArray>::New();

  _Cells->InitTraversal();

  vtkSmartPointer<vtkIdList> cell =
      vtkSmartPointer<vtkIdList>::New();
  int cellId=0;

  int oldNumCells =_Output->GetNumberOfPolys();

  while (_Cells->GetNextCell(cell)){

    if (cellId >= oldNumCells || !_DeletedCells[cellId]){
      for (int i=0; i<cell->GetNumberOfIds(); i++){
        int oldPtId = cell->GetId(i);
        int newPtId = mapping[oldPtId];

        if (newPtId == -1){
          cout << "old id" << oldPtId << endl;
          cout << "ERROR - Remeshe::UpdateCells invalid point mapping" << endl;
          exit(1);
        }

        cell->SetId(i, newPtId);
      }

      newCells->InsertNextCell(cell);
    }
    cellId++;
  }

  _Cells->DeepCopy(newCells);

}

void Remesher::Finalize(){
  vtkSmartPointer<vtkCellArray> newCells =
      vtkSmartPointer<vtkCellArray>::New();


  Array<int> mapping = UpdatePoints();
  UpdateCells(mapping);

  vtkSmartPointer<vtkPolyData> newOutput =
      vtkSmartPointer<vtkPolyData>::New();

  newOutput->Allocate();
  newOutput->SetPoints(_Points);
  newOutput->SetPolys(_Cells);
  newOutput->BuildCells();
  newOutput->BuildLinks();

  _Output->DeepCopy(newOutput);

}


bool Remesher::CheckTopology(){
  bool topolgyErrors = false;

  vtkSmartPointer<vtkDoubleArray> badPtIds =
      vtkSmartPointer<vtkDoubleArray>::New();

  badPtIds->SetName("badPtIds");

  int numPts = _Output->GetNumberOfPoints();
  badPtIds->SetNumberOfTuples(numPts);
  badPtIds->SetNumberOfComponents(1);
  for (int ptId=0; ptId<numPts; ptId++){
    badPtIds->SetTuple1(ptId, 0);
  }

  EdgeTable edges(_Output);
  EdgeIterator ei(edges);

  ei.InitTraversal();
  int ptId1, ptId2;



  while ( ei.GetNextEdge(ptId1, ptId2) != -1 ){

    vtkSmartPointer<vtkIdList> cellIds1, cellIds2;
    cellIds1 = vtkSmartPointer<vtkIdList>::New();
    cellIds2 = vtkSmartPointer<vtkIdList>::New();
    _Output->GetPointCells(ptId1, cellIds1);
    _Output->GetPointCells(ptId2, cellIds2);

    vtkSmartPointer<vtkIdList> sharedCellIds =
        vtkSmartPointer<vtkIdList>::New();
    sharedCellIds->DeepCopy(cellIds1);
    sharedCellIds->IntersectWith(cellIds2);


    if (sharedCellIds->GetNumberOfIds() != 2){
      badPtIds->SetTuple1(ptId1, 1);
      badPtIds->SetTuple1(ptId2, 1);
      topolgyErrors = true;
    }


  }

  _Output->GetPointData()->AddArray(badPtIds);

  return topolgyErrors;
}



void Remesher::Run(){

  if (_MinEdgeLength > _MaxEdgeLength){ return; }

  _Output = _Input;

  for (int i=0; i<_MaxIterations; i++){
    Initialize();
    ComputeEdgeLengths();
    QueueEdgesByLength();
    int numOfModifications = ProcessEdgeQueue();
    Finalize();

    if (_Debug){
      bool topologyError = CheckTopology();
      string outName = std::to_string(i) + "-out.vtk";
      WritePolyData(outName.c_str() , _Output);
      if (topologyError){
        cout << "topology error" << endl;
        exit(1);
      }

    }

    if (numOfModifications == 0) {
      cout << "No edges were modified!" << endl;
      break;
    }
  }

}




} // mirtk namespace





