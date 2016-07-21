//my mirtk includes
#include "mirtk/M2SDiffuser.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointSetIO.h"

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================

M2SDiffuser::M2SDiffuser(){
  _AdaptiveTermination = true;
  _MaxIterations = 1000;

  _Input = vtkSmartPointer<vtkPolyData>::New();
  _Output = vtkSmartPointer<vtkPolyData>::New();

  _Debug = false;
}


// -----------------------------------------------------------------------------


M2SDiffuser::M2SDiffuser(const M2SDiffuser &other)
{
  *this = other;
}

// -----------------------------------------------------------------------------

M2SDiffuser &M2SDiffuser::operator =(const M2SDiffuser &other)
{
  if (this != &other) {
    M2SDiffuser::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
M2SDiffuser::~M2SDiffuser(){
}

void M2SDiffuser::CopyAttributes(const M2SDiffuser &other){
  cout << "TODO copy attributes" << endl;
}

// =============================================================================
// Class methods
// =============================================================================


void M2SDiffuser::Run(){
  //compute edge lengths

  _Output->DeepCopy(_Input);

  int numPts = _Output->GetNumberOfPoints();

  int numSourcePts = 0;

  for (int ptId=0; ptId< numPts; ptId++){
    numSourcePts += _SourceMask[ptId];
  }

  EdgeTable edges(_Output);
  int numEdges = edges.NumberOfEdges();


  //TODO implement different weighting schemes
  Array<double> edgeWeights;
  for (int edgeId=0; edgeId<numEdges; edgeId++){
    edgeWeights.push_back(1.0);
  }






  Array<double> totalWeighting;


  Array<bool> previouslyVisitedMask = _SourceMask;
  Array<bool> visitedMask = _SourceMask;
  Array<double> weightsum(numPts, 0.0);


  //zero non active points
   //oldpt[3];
  vtkSmartPointer<vtkPoints> newPts =
      vtkSmartPointer<vtkPoints>::New();

  for (int i=0; i<numPts; i++){
    double p[3] = {0,0,0};
    newPts->InsertNextPoint(p);
  }


  //TODO remove this loop
  //zeroing non source points for visualisation when outputting mesh files

  for (int ptId=0; ptId< numPts; ptId++){
    if(!_SourceMask[ptId]){
      double p[3] = {0,0,0};
      _Output->GetPoints()->SetPoint(ptId,p);
    }
  }



  bool visitedAllPts = false;

  EdgeIterator ei(edges);

  int maxIterations = _MaxIterations;

  for (int iteration=0; iteration<maxIterations; iteration++){
    ei.InitTraversal();
    int edgeId = 0;
    int pid1, pid2;
    double pt[3], oldpt [3];
    while ( ei.GetNextEdge(pid1, pid2) != -1 ){
      //cout << "processing edge " << edgeId << endl;

      double w = edgeWeights[edgeId];


      if (w == 0.0){
        cout << "zero weight" << endl;
        exit(1);
      }

      //modify pt 1 if it's not a source pt (only if pt 2 has been visited already)
      if (! _SourceMask[pid1] && previouslyVisitedMask[pid2]){
        _Output->GetPoint(pid2, oldpt);
        newPts->GetPoint(pid1, pt);

         for (int i=0; i<3; i++){ pt[i] += oldpt[i]*w; }
         newPts->SetPoint(pid1, pt);
         weightsum[pid1] += w;

         if(!visitedMask[pid1])  visitedMask[pid1] = true;
      }

      //modify pt 2 if it's not a source pt (only if pt 2 has been visited already)
      if (! _SourceMask[pid2] && previouslyVisitedMask[pid1]){
        _Output->GetPoint(pid1, oldpt);
         newPts->GetPoint(pid2, pt);

         for (int i=0; i<3; i++){ pt[i] += oldpt[i]*w; }
         newPts->SetPoint(pid2, pt);
         weightsum[pid2] += w;

         if(!visitedMask[pid2])  visitedMask[pid2] = true;
      }
      edgeId++;
    }

    //update non-source points
    int totalVisited = 0;
    for (int ptId=0; ptId<numPts; ptId++){
      if (! _SourceMask[ptId] && visitedMask[ptId]){

        newPts->GetPoint(ptId, pt);




        if (weightsum[ptId] == 0){  cout << "something bad happened" << endl; exit(1); }

        for (int i=0; i<3; i++)  { pt[i] /= weightsum[ptId]; }

        weightsum[ptId] = 0;
        previouslyVisitedMask[ptId] = true;

        //update


        //_Output->GetPoint(ptId, oldpt);
        //for (int i=0; i<3; i++)  { pt[i] = (pt[i] + oldpt[i])/2; }

        _Output->GetPoints()->SetPoint(ptId, pt);

        if (pt[0] == 0 || pt[1] == 0 || pt[2] == 0){
          cout << "ptId" << ptId << endl;
          cout << "zero coorindate!" << endl;
          exit(1);
        }

        if (pt[0] > 2 || pt[1] > 2 || pt[2] > 2){
          cout << "ptId" << ptId << endl;
          cout << "coordinate > 2!" << endl;
          exit(1);
        }

        if (pt[0] < -2 || pt[1] < -2 || pt[2] < -2){
          cout << "ptId" << ptId << endl;
          cout << "coordinate < -2!" << endl;
          exit(1);
        }


        //reset point for next iteration
        for (int i=0; i<3; i++) { pt[i]=0; }
        newPts->SetPoint(ptId, pt);
        totalVisited++;
        }
      }


    if (!visitedAllPts  && _AdaptiveTermination && totalVisited + numSourcePts == numPts){
      //cout << "all points visited!" << endl;
      visitedAllPts = true;
      //modify max iterations based on time to visit all points
      int potentialMaxIterations = iteration*20;
      if (_MaxIterations > potentialMaxIterations)
        maxIterations = potentialMaxIterations;
    }

    if (_Debug) {
      string outName = "diffuser-" + std::to_string(iteration) + ".vtk";
      WritePolyData(outName.c_str(), _Output);
    }

  }



}


} // mirtk namespace



