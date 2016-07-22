/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/MeshToSphere.h"

#include "mirtk/Common.h"
#include "mirtk/Options.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointSetIO.h"




using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "Usage: " << name << " <input> <output> -parin <parameter_file>" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Robert Wright's tool for embedding a mesh onto a sphere" << endl;
  cout << endl;
  cout << "Input options:" << endl;

  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================


// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  verbose = 1; // default verbosity level
  EXPECTS_POSARGS(2);

  // Read input point set
  FileOption fopt = FO_ASCII;
  vtkSmartPointer<vtkPolyData> input = ReadPolyData(POSARG(1));

  // Optional arguments
  string parin = "";
  bool debug = false;

  for (ALL_OPTIONS) {
    // Input
    if (OPTION("-parin") || OPTION("-pi")) parin = ARGUMENT;
    else if (OPTION("-debug") || OPTION("-d")) debug = true;      
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  ParameterParser parser;
  if (parin != ""){
    parser.ParseFile(parin);
  }
  else {
    //TODO implement default parameter set based on input mesh
    cout << "Default parameter not implemented yet!" << endl;
    cout << "Please specify parameter file with -parin flag." << endl;
    exit(1);
	//parameters.Initialize();
  }

  cout << "Embedding mesh with the following parameters: " << endl;
  cout << endl;
  cout << parser.ToString();
  cout << endl;

/*
  MeshToSphere::Params p0;
  p0.minEdgeLength = 4;
  p0.maxEdgeLength = p0.minEdgeLength*2.2;
  p0.connectivityDistanceThreshold = -1.0;
  p0.maxIterations = 100;

  MeshToSphere::Params p1;
  p1.minEdgeLength = 2;
  p1.maxEdgeLength = p1.minEdgeLength*2.2;
  p1.connectivityDistanceThreshold = 20;
  p1.maxIterations = 50;

  MeshToSphere::Params p2;
  p2.minEdgeLength = 1;
  p2.maxEdgeLength = p2.minEdgeLength*2.2;
  p2.connectivityDistanceThreshold = 10;
  p2.maxIterations = 25;

  Array<MeshToSphere::Params> params;
  //params.push_back(p0);
  params.push_back(p1);
  //params.push_back(p2);
*/

  MeshToSphere m2s;
  m2s.Input(input);
  m2s.Parameters(parser.GetParameters());
  m2s.Debug(debug);
  m2s.Run();
  vtkSmartPointer<vtkPolyData> output = m2s.Output();

  // Center output point set
  //if (center_output)

  // Scale output surface to match input area
  //if (match_area) Scale(output, sqrt(Area(input) / Area(output)));

  //Write output surface
  return WritePolyData(POSARG(2), output, fopt);
}
