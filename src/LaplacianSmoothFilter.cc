#include "mirtk/LaplacianSmoothFilter.h"
#include "mirtk/Common.h"

#include "math.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"


LaplacianSmoothFilter::LaplacianSmoothFilter(){
  _iterations=1;
  _lambda=1;
  _sigma=0;
  _pointWeightingName=NULL;
}
LaplacianSmoothFilter::~LaplacianSmoothFilter(){}


void LaplacianSmoothFilter::SetInput(vtkSmartPointer<vtkPolyData> input){
  _input=input;
}

vtkSmartPointer<vtkPolyData> LaplacianSmoothFilter::GetOutput(){
  return _output;
}


void LaplacianSmoothFilter::SetIterations(int iterations){
 _iterations=iterations;
}


void LaplacianSmoothFilter::SetLambda(double lambda){
  _lambda=lambda;
}
void LaplacianSmoothFilter::SetSigma(double sigma){
  _sigma=sigma;
}

void LaplacianSmoothFilter::SetPointWeighting(char * name){
  _pointWeightingName=name;
}

void LaplacianSmoothFilter::Update(){
  int i,j,k, pIdx0,pIdx1,pIdx2,numPoints,numFaces;
  double p [3], p0 [3], p1 [3],p2 [3],np [3],np0 [3], np1 [3],np2 [3];
  double e0x, e0y, e0z, e0l, e1x, e1y, e1z, e1l, e2x, e2y, e2z, e2l; //, ux, uy, uz;
  double pw;
  vtkSmartPointer<vtkPoints> oldPoints, newPoints, newPoints2;
  vtkSmartPointer<vtkCellArray> faces;
  vtkSmartPointer<vtkIdList> pts;

  oldPoints=_input->GetPoints();
  numPoints=oldPoints->GetNumberOfPoints();

  double weights [numPoints];



  newPoints=vtkSmartPointer<vtkPoints>::New();
  newPoints2=vtkSmartPointer<vtkPoints>::New();
  newPoints->DeepCopy(oldPoints);
  newPoints2->DeepCopy(oldPoints);

  faces=_input->GetPolys();
  numFaces=faces->GetNumberOfCells();


  vtkDoubleArray * pointWeighting = NULL;
  if (_pointWeightingName != NULL){
    if (_input->GetPointData()->HasArray(_pointWeightingName)){
      pointWeighting = (vtkDoubleArray *)
          _input->GetPointData()->GetArray(_pointWeightingName);
    }
  }



   pts=vtkSmartPointer<vtkIdList>::New();


   for (j=0; j<_iterations; j++) {

       // Clean the weights
       for (i=0; i<numPoints; i++) {
         weights[i] =0;
         newPoints2->SetPoint(i,0,0,0);
       }
       faces->InitTraversal();
       /* Calculate all face normals and angles */


       for (i=0; i<numFaces; i++){

           /* Get indices of face vertices */
         faces->GetNextCell(pts);

           pIdx0=pts->GetId(0);
           pIdx1=pts->GetId(1);
           pIdx2=pts->GetId(2);

           /* Calculate edge lengths */
           newPoints->GetPoint(pIdx0,p0);
           newPoints->GetPoint(pIdx1,p1);
           newPoints->GetPoint(pIdx2,p2);

           e0x=p0[0]-p1[0];
           e0y=p0[1]-p1[1];
           e0z=p0[2]-p1[2];
           e1x=p1[0]-p2[0];
           e1y=p1[1]-p2[1];
           e1z=p1[2]-p2[2];
           e2x=p2[0]-p0[0];
           e2y=p2[1]-p0[1];
           e2z=p2[2]-p0[2];
           e0l=1 / (sqrt(e0x*e0x + e0y*e0y + e0z*e0z)+_sigma);
           e1l=1 / (sqrt(e1x*e1x + e1y*e1y + e1z*e1z)+_sigma);
           e2l=1 / (sqrt(e2x*e2x + e2y*e2y + e2z*e2z)+_sigma);

           newPoints2->GetPoint(pIdx0,np0);
           newPoints2->GetPoint(pIdx1,np1);
           newPoints2->GetPoint(pIdx2,np2);

           for (k=0; k<3; k++){
             np0[k]+=p1[k]*e0l;
             np1[k]+=p0[k]*e0l;
             np1[k]+=p2[k]*e1l;
             np2[k]+=p1[k]*e1l;
             np2[k]+=p0[k]*e2l;
             np0[k]+=p2[k]*e2l;
           }

           weights[pIdx0]+=e0l;
           weights[pIdx1]+=e0l;
           weights[pIdx1]+=e1l;
           weights[pIdx2]+=e1l;
           weights[pIdx2]+=e2l;
           weights[pIdx0]+=e2l;

           newPoints2->SetPoint(pIdx0,np0);
           newPoints2->SetPoint(pIdx1,np1);
           newPoints2->SetPoint(pIdx2,np2);

       }

       // Normalize the Vertices */
       for (i=0; i<numPoints; i++){
         if (weights[i] == 0){
           cout << "WARNING: point " << i << " has no adjacencies!" << endl;
           for (k=0; k<3; k++)
             np[k]=p[k];
         }
         else{
           newPoints2->GetPoint(i,np);
           newPoints->GetPoint(i,p);
           /*
           ux=0; uy=0; uz=0;
           ux=np[0]/weights[i];
           uy=np[1]/weights[i];
           uz=np[2]/weights[i];

           ux=ux-p[0];
           uy=uy-p[1];
           uz=uz-p[2];

           np[0]=p[0]+ux*_lambda;
           np[1]=p[1]+uy*_lambda;
           np[2]=p[2]+uz*_lambda;


           */

           if (pointWeighting == NULL){
             for (k=0; k<3; k++){
               np[k]=p[k]+((np[k]/weights[i])-p[k])*_lambda;
             }
           }
           else{
             pw = pointWeighting->GetTuple1(i);
             for (k=0; k<3; k++){
               np[k]=p[k]+((np[k]/weights[i])-p[k])*_lambda*pw;
             }
           }




           newPoints2->SetPoint(i,np);
         }
       }

       newPoints->DeepCopy(newPoints2);

   }
   _output=vtkSmartPointer<vtkPolyData>::New();
   _output->DeepCopy(_input);
   _output->SetPoints(newPoints);






}


