/*
 * Remesher.cc
 *
 *  Created on: 13 Feb 2015
 *      Author: rob
 */

#include "mirtk/Common.h"
#include "mirtk/MeshToSphereRemesher.h"
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
#include <vtkKdTreePointLocator.h>


double Remesher::GetEdgeLength(double p1 [],double p2 []){
  double length=0;
  for (int i=0;i<3;i++) length+=pow((p2[i]-p1[i]),2);
  return sqrt(length);
}

void Remesher::Midpoint(double p1 [],double p2 [],double mp []){
  for (int i=0;i<3;i++) mp[i]=(p1[i]+p2[i])/2;
}


vtkSmartPointer<vtkPolyData> Remesher::ProcessInput(vtkPolyData * poly){

  vtkSmartPointer<vtkTriangleFilter> triangulator
      = vtkSmartPointer<vtkTriangleFilter>::New();

  triangulator->SetInputData(poly);
  triangulator->Update();

  triangulator->GetOutput()->BuildCells();
  triangulator->GetOutput()->BuildLinks();
  return triangulator->GetOutput();

}

vtkSmartPointer<vtkPolyData> Remesher::BuildEdgeStructure(vtkPolyData * poly){

  vtkSmartPointer<vtkExtractEdges> ee =
      vtkSmartPointer<vtkExtractEdges>::New();
  ee->SetInputData(poly);
  ee->Update();
  vtkSmartPointer<vtkPolyData> edges = ee->GetOutput();
  edges->BuildCells();
  return edges;
}

//Need to parallalise this!
vtkSmartPointer<vtkDoubleArray> Remesher::ComputeEdgeLengths(
    vtkPolyData * edges){
  double p0 [3], p1 [3];
  int id0,id1;
  double d;

  vtkSmartPointer<vtkDoubleArray> edgeLengths =
      vtkSmartPointer<vtkDoubleArray>::New();

  vtkSmartPointer<vtkIdList> edge =
      vtkSmartPointer<vtkIdList>::New();

  int numEdges = edges->GetNumberOfLines();
  edgeLengths->SetNumberOfComponents(1);
  edgeLengths->SetNumberOfTuples(numEdges);
  edgeLengths->SetName("Length");

  for (int edgeId=0; edgeId<numEdges;edgeId++){
      edges->GetCellPoints(edgeId,edge);
      id0 = edge->GetId(0);
      id1 = edge->GetId(1);
      edges->GetPoint(id0,p0);
      edges->GetPoint(id1,p1);
      d = GetEdgeLength(p0,p1);
      edgeLengths->InsertTuple1(edgeId,d);
  }
  return edgeLengths;
}


vtkSmartPointer<vtkPriorityQueue> Remesher::QueueEdgesByLength(
    vtkPolyData * edges, double minEdgeLength, double maxEdgeLength){

  int splitNum = 0, collapseNum = 0;
  int numEdges = edges->GetNumberOfLines();
  double priority,d;

  vtkSmartPointer<vtkPriorityQueue> edgeQueue =
      vtkSmartPointer<vtkPriorityQueue>::New();

  vtkSmartPointer<vtkDoubleArray> edgeLengths = ComputeEdgeLengths(edges);

  vtkSmartPointer<vtkIntArray> isCollapse
      = vtkSmartPointer<vtkIntArray>::New();
  isCollapse->SetNumberOfComponents(1);
  isCollapse->SetNumberOfTuples(numEdges);
  isCollapse->SetName("isCollapse");

  //now queue edges

  for (int edgeId=0; edgeId<numEdges;edgeId++){
      d = edgeLengths->GetTuple1(edgeId);
      if (d > maxEdgeLength){
          priority = maxEdgeLength / (d-maxEdgeLength);
          edgeQueue->Insert(priority,edgeId);
          splitNum+=1;
          isCollapse->SetTuple1(edgeId,0);
      }
      else if (d < minEdgeLength){
          priority = minEdgeLength / (minEdgeLength-d);
          edgeQueue->Insert(priority,edgeId);
          collapseNum += 1;
          isCollapse->SetTuple1(edgeId,1);
      }
  }
  cout << "Attempting to split " << splitNum << " edges and collapse " << collapseNum
      << " edges (" << numEdges << " total)." << endl;
  edges->GetCellData()->AddArray(isCollapse);
  return edgeQueue;
}


vtkSmartPointer<vtkIdList> Remesher::GetAdjacentCellIds(int pId0,int pId1,
    vtkPolyData * poly){

  vtkSmartPointer<vtkIdList> cellList0,cellList1;
  cellList0 = vtkSmartPointer<vtkIdList>::New();
  cellList1 = vtkSmartPointer<vtkIdList>::New();

  poly->GetPointCells(pId0,cellList0);
  poly->GetPointCells(pId1,cellList1);
  cellList0->IntersectWith(cellList1); // common cells that both points are in

  return cellList0;
}

void Remesher::DeleteEdgesAfterCollapse(int pId0,int pId1,
    vtkPolyData *  edges){
  int i,neighbourEdgeId;
  vtkSmartPointer<vtkIdList> edgeList0,edgeList1;
  edgeList0 = vtkSmartPointer<vtkIdList>::New();
  edgeList1 = vtkSmartPointer<vtkIdList>::New();

  edges->GetPointCells(pId0,edgeList0);
  edges->GetPointCells(pId1,edgeList1);
  for (i=0;i<edgeList0->GetNumberOfIds();i++){
    neighbourEdgeId = edgeList0->GetId(i);
    edges->DeleteCell(neighbourEdgeId);
  }
  for (i=0;i<edgeList1->GetNumberOfIds();i++){
    neighbourEdgeId = edgeList1->GetId(i);
    edges->DeleteCell(neighbourEdgeId);
  }
}

void DeleteEdgesAfterSplit(int pId0,int pId1,int otherPId0,
    int otherPId1, vtkPolyData *  edges){



  int i,edgeId;
  vtkSmartPointer<vtkIdList> edgeList0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> edgeList1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherEdgeList0
      = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherEdgeList1
      = vtkSmartPointer<vtkIdList>::New();


  edges->GetPointCells(pId0,edgeList0);
  edges->GetPointCells(pId1,edgeList1);
  edges->GetPointCells(otherPId0,otherEdgeList0);
  edges->GetPointCells(otherPId1,otherEdgeList1);

  //combine edgeLists

  for (i=0;i<edgeList1->GetNumberOfIds();i++){
    edgeId = edgeList1->GetId(i);
    edgeList0->InsertUniqueId(edgeId);
  }

  for (i=0;i<otherEdgeList1->GetNumberOfIds();i++){
    edgeId = otherEdgeList1->GetId(i);
    otherEdgeList0->InsertUniqueId(edgeId);
  }

  edgeList0->IntersectWith(otherEdgeList0);



  for (i=0;i<edgeList0->GetNumberOfIds();i++){
    edgeId = edgeList0->GetId(i);
    edges->DeleteCell(edgeId);
  }
}

void Remesher::CreateNewCells(int oldCellId0, int oldCellId1,int pId0,int pId1,
    int newPId, vtkCellArray * newCells, vtkPolyData * poly){

  int i;
  vtkSmartPointer<vtkIdList> oldCellPIds0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> oldCellPIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> newCell = vtkSmartPointer<vtkIdList>::New();

  poly->GetCellPoints(oldCellId0,oldCellPIds0);
  poly->GetCellPoints(oldCellId1,oldCellPIds1);

  newCell->DeepCopy(oldCellPIds0);
  i = newCell->IsId(pId0);
  newCell->SetId(i,newPId);
  newCells->InsertNextCell(newCell);

  newCell->DeepCopy(oldCellPIds0);
  i = newCell->IsId(pId1);
  newCell->SetId(i,newPId);
  newCells->InsertNextCell(newCell);

  newCell->DeepCopy(oldCellPIds1);
  i = newCell->IsId(pId0);
  newCell->SetId(i,newPId);
  newCells->InsertNextCell(newCell);

  newCell->DeepCopy(oldCellPIds1);
  i = newCell->IsId(pId1);
  newCell->SetId(i,newPId);
  newCells->InsertNextCell(newCell);
}

void Remesher::CombineCells(vtkPolyData * poly, vtkCellArray * newCells,
    vtkIntArray * toDelete){

  //add the new cells we created
  vtkSmartPointer<vtkIdList> newCell = vtkSmartPointer<vtkIdList>::New();

  for(int i = 0; i < poly->GetNumberOfCells(); i++){
    if (toDelete->GetTuple1(i) == 0){
      poly->GetCellPoints(i,newCell);
      newCells->InsertNextCell(newCell);
    }
  }
  newCells->Squeeze();
}

void Remesher::AveragePointdata(int newPId,int pId0,int pId1,
    vtkPolyData * poly){
  //TODO
}

vtkSmartPointer<vtkIdList> GetOtherPoints(int pId0,int pId1,int oldCellId0,
    int oldCellId1,vtkPolyData * poly){
  int i,otherId;


  vtkSmartPointer<vtkIdList> oldCellPIds0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> oldCellPIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherPIds = vtkSmartPointer<vtkIdList>::New();
  otherPIds->Allocate(2);

  poly->GetCellPoints(oldCellId0,oldCellPIds0);
  poly->GetCellPoints(oldCellId1,oldCellPIds1);

  for (i=0;i<3;i++){
    otherId = oldCellPIds0->GetId(i);
    if (otherId != pId0 && otherId != pId1){
      otherPIds->InsertNextId(otherId);
      break;
    }
  }

  for (i=0;i<3;i++){
    otherId = oldCellPIds1->GetId(i);
    if (otherId != pId0 && otherId != pId1){
      otherPIds->InsertNextId(otherId);
      break;
    }
  }
  return otherPIds;
}


vtkSmartPointer<vtkPolyData> Remesher::ProcessEdgeQueue(
    vtkPriorityQueue * edgeQueue, vtkPolyData * edges, vtkPolyData * poly){
  int collapseCount = 0, splitCount = 0;
  int edgeId;
  int pId0,pId1,newPId;
  int edgePId0, edgePId1, edgePId2, edgePId3; //NEW CODE
  int oldCellId0,oldCellId1;
  double p0 [3],p1 [3],m [3];

  double p2 [3], p3 [3]; //NEW CODE



  vtkSmartPointer<vtkPolyData> polyCopy = vtkSmartPointer<vtkPolyData>::New();
  polyCopy->DeepCopy(poly);

  //*NEW CODE

  vtkSmartPointer<vtkKdTreePointLocator> polyKDTree =
    vtkSmartPointer<vtkKdTreePointLocator>::New();
  polyKDTree->SetDataSet(polyCopy);
  polyKDTree->BuildLocator();

  vtkSmartPointer<vtkKdTreePointLocator> edgesKDTree =
    vtkSmartPointer<vtkKdTreePointLocator>::New();
  edgesKDTree->SetDataSet(edges);
  edgesKDTree->BuildLocator();
  cout << "finished building trees" << endl;
  //*END NEW CODE



  vtkSmartPointer<vtkIdList> edge = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> oldCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherPIds = vtkSmartPointer<vtkIdList>::New();

  vtkSmartPointer<vtkCellArray>  newCells =
      vtkSmartPointer<vtkCellArray>::New();
  newCells->Allocate(polyCopy->GetNumberOfCells() + edgeQueue->GetNumberOfItems()*4);
  //TODO polyCopy num cells change to number of cells to split * 4

  vtkIntArray * isCollapse =
      (vtkIntArray *) edges->GetCellData()->GetArray("isCollapse");

  vtkSmartPointer<vtkIntArray> cellsToDelete =
      vtkSmartPointer<vtkIntArray>::New();
  for (int i = 0; i < polyCopy->GetNumberOfCells();i++){
    cellsToDelete->InsertNextTuple1(0);
  }



  while (true){
    edgeId = edgeQueue->Pop();
    if (edgeId < 0) //break loop when we've run out of edges
        break;

    //step 1 if edge has been marked as deleted skip
    if (edges->GetCellType(edgeId) == 0)
      continue;

    //step 2 find the pids for the edge points and the cells that border it
    edges->GetCellPoints(edgeId,edge);



    //pId0 = edge->GetId(0); pId1 = edge->GetId(1); //OLD CODE



    //*NEW CODE - find the original mesh pids from the edge pids
    //cannot assume they are the same, due to the edge extraction filter code!

    edgePId0 = edge->GetId(0); edgePId1 = edge->GetId(1);

    edges->GetPoint(edgePId0,p0);
    pId0 = polyKDTree->FindClosestPoint(p0);

    edges->GetPoint(edgePId1,p1);
    pId1 = polyKDTree->FindClosestPoint(p1);

    //*END NEW CODE

    oldCellIds = GetAdjacentCellIds(pId0,pId1,polyCopy);

    //if (oldCellIds->GetNumberOfIds() < 2){
    //cout << "ERROR - found an edge with less than 2 adajacent faces" << endl;
    //continue;

    //There are cases where the two edge points have more than two edges with
    //common points. We will skip this for now, because its extremely rare and
    //much more complex.
    if (oldCellIds->GetNumberOfIds() != 2)
        continue;

    oldCellId0 = oldCellIds->GetId(0);
    oldCellId1 = oldCellIds->GetId(1);

    //step 3 calculate edge midpoint
    //polyCopy->GetPoint(pId0,p0);    polyCopy->GetPoint(pId1,p1);
    Midpoint(p0,p1,m);

    //step 4 determine if we're collapsing or splitting
    if (isCollapse->GetTuple1(edgeId) == 1){ //collapse

      // Move edge points to midpoint, effectively collapsing it!
      polyCopy->GetPoints()->SetPoint(pId0,m);
      polyCopy->GetPoints()->SetPoint(pId1,m);

      //mark neighbouring edges as "deleted" so they aren't processed in the
      //next iteration

      //DeleteEdgesAfterCollapse(pId0,pId1,edges); //OLD CODE
      DeleteEdgesAfterCollapse(edgePId0,edgePId1,edges); //NEW CODE

      collapseCount++;
    }

    else{ //split
        //Insert mid point as a new point and create new cells
        newPId = polyCopy->GetPoints()->InsertNextPoint(m);
        CreateNewCells(oldCellId0,oldCellId1,pId0,pId1,newPId,newCells,polyCopy);
        AveragePointdata(newPId,pId0,pId1,polyCopy);

        //mark neighbouring edges as "deleted" so they aren't processed in the
        //next iteration


        //*OLD CODE
        //otherPIds = GetOtherPoints(pId0,pId1,oldCellId0,oldCellId1,poly);
        //DeleteEdgesAfterSplit(pId0,pId1,otherPIds->GetId(0), otherPIds->GetId(1),edges);
        //*END OLD CODE

        //NEW CODE
        otherPIds = GetOtherPoints(pId0,pId1,oldCellId0,oldCellId1,poly);

        edges->GetPoint(otherPIds->GetId(0),p2);
        edges->GetPoint(otherPIds->GetId(1),p3);

        edgePId2 = edgesKDTree->FindClosestPoint(p2);
        edgePId3 = edgesKDTree->FindClosestPoint(p3);

        DeleteEdgesAfterSplit(edgePId0,edgePId1,edgePId2,edgePId3,edges); // OLD CODE
        //END NEW CODE

        edges->DeleteCell(edgeId);

        splitCount++;
    }

    //step 6 finally mark the old cells for deletion
    cellsToDelete->InsertTuple1(oldCellId0,1);
    cellsToDelete->InsertTuple1(oldCellId1,1);
  }

  cout << "Successfully split " << splitCount  << " edges and collapsed " << collapseCount <<
      " edges." <<endl;


  //create new polydata
  vtkSmartPointer<vtkPolyData> newPoly =
      vtkSmartPointer<vtkPolyData>::New();
  cout << polyCopy->GetNumberOfCells() << endl;

  CombineCells(polyCopy,newCells,cellsToDelete);
  newPoly->Allocate();
  newPoly->SetPoints(polyCopy->GetPoints());
  newPoly->SetPolys(newCells);

  vtkSmartPointer<vtkCleanPolyData> cleaner
      = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(newPoly);
  cleaner->PointMergingOn();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(.0);
  cleaner->GetOutput()->SetVerts(NULL);
  cleaner->GetOutput()->SetLines(NULL);
  cleaner->Update();

  vtkSmartPointer<vtkPolyData> output =
      RemoveDegenerateCells(cleaner->GetOutput());

  return output;


  /*
  poly->Modified();
  poly->BuildCells();
  AddNewCells(poly,newCells);
  poly->BuildCells();



  */
  //TODO pointdata is wiped because we
  //now have a different number of points. AveragePointdata() will address this
  //in the future.

}



vtkSmartPointer<vtkPolyData> Remesher::RemoveDegenerateCells(vtkPolyData * polydata){
    int i,j,cid,pid,pid0,pid1,pid2;

    polydata->BuildCells();
    polydata->BuildLinks();

    //build adjacencies
    int numCells = polydata->GetNumberOfCells();
    vtkIdList ** edges = new vtkIdList * [numCells];
    for (cid=0; cid<numCells; cid++)  edges[cid] = vtkIdList::New();

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    //build adjcencies
    for (cid=0; cid<numCells; cid++){
      polydata->GetCellPoints(cid, ptIds);
      pid0 = ptIds->GetId(0);
      pid1 = ptIds->GetId(1);
      pid2 = ptIds->GetId(2);

      edges[pid0]->InsertUniqueId(pid1);
      edges[pid0]->InsertUniqueId(pid2);

      edges[pid1]->InsertUniqueId(pid0);
      edges[pid1]->InsertUniqueId(pid2);

      edges[pid2]->InsertUniqueId(pid1);
      edges[pid2]->InsertUniqueId(pid0);
    }


    //initialize delete array
    vtkSmartPointer<vtkIntArray> cellsToDelete =
        vtkSmartPointer<vtkIntArray>::New();
    for (cid=0; cid<numCells; cid++) cellsToDelete->InsertNextTuple1(0);


    int numPts = polydata->GetNumberOfPoints();
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    for (pid=0; pid<numPts; pid++){
      if (edges[pid]->GetNumberOfIds() <= 2){
        //degenerate point
        polydata->GetPointCells(pid,cellIdList);
        for (j=0; j<cellIdList->GetNumberOfIds(); j++){
          //degenerate cell
          cid = cellIdList->GetId(j);
          cellsToDelete->SetTuple1(cid,1);
        }
      }
    }

    vtkSmartPointer<vtkIdList> newCell =
        vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkCellArray> newCells =
        vtkSmartPointer<vtkCellArray>::New();

    for (cid=0; cid<numCells; cid++){
      if (cellsToDelete->GetTuple1(cid) == 0){
        polydata->GetCellPoints(cid,newCell);
        newCells->InsertNextCell(newCell);
      }
    }
    newCells->Squeeze();
    vtkSmartPointer<vtkPolyData> newPolydata =
        vtkSmartPointer<vtkPolyData>::New();
    newPolydata->Allocate();
    newPolydata->SetPoints(polydata->GetPoints());
    newPolydata->SetPolys(newCells);

    //run cleaner to remove lone points
    vtkSmartPointer<vtkCleanPolyData> cleaner =
        vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(newPolydata);
    cleaner->ConvertLinesToPointsOn();
    cleaner->ConvertPolysToLinesOn();
    cleaner->ConvertStripsToPolysOn();
    cleaner->PointMergingOn();
    cleaner->SetAbsoluteTolerance(0);
    cleaner->Update();
    cleaner->GetOutput()->SetVerts(NULL);
    cleaner->GetOutput()->SetLines(NULL);

    for (cid=0; cid<numCells; cid++)  edges[cid]->Delete();
    delete [] edges;


    return cleaner->GetOutput();
}




Remesher::Remesher(){
  _maxEdgeLength=DBL_MAX;
  _minEdgeLength=-1;
  _maxIterations=1;
}

Remesher::~Remesher(){
}

void Remesher::SetMaxEdgeLength(double maxEdgeLength){
  _maxEdgeLength=maxEdgeLength;
}
void Remesher::SetMinEdgeLength(double minEdgeLength){
  _minEdgeLength=minEdgeLength;
}

void Remesher::SetMaxIterations(int iterations){
  _maxIterations=iterations;
}

vtkSmartPointer<vtkPolyData> Remesher::Remesh(vtkPolyData * poly){
  cout << "Remeshing" << endl;


  vtkSmartPointer<vtkPolyData> remeshed = ProcessInput(poly);


  for (int i=0;i<_maxIterations;i++){
    remeshed->GetPointData()->Initialize();
    remeshed->GetCellData()->Initialize();
    vtkSmartPointer<vtkPolyData> edges = BuildEdgeStructure(remeshed);
    vtkSmartPointer<vtkPriorityQueue> edgeQueue = QueueEdgesByLength(edges,
        _minEdgeLength,_maxEdgeLength);

    if (edgeQueue->GetNumberOfItems() == 0) break;

    remeshed = ProcessEdgeQueue(edgeQueue,edges,remeshed);
  }


  return remeshed;
 }







