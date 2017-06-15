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


// #include "mirtk/DeformableSurfaceModel.h"
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
  cout << "  -inflated <inflated_mesh>    use inflated mesh to aid decimation (must have same topology as input mesh):" << endl;

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
  string parin = "", inflated_name = "";
  bool debug = false;

  for (ALL_OPTIONS) {
    // Input
    if (OPTION("-parin") || OPTION("-pi")) parin = ARGUMENT;
    else if (OPTION("-debug") || OPTION("-d")) debug = true;
    else if (OPTION("-inflated") || OPTION("-i")) inflated_name = ARGUMENT;

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

  //Read inflated mesh if given
  vtkSmartPointer<vtkPolyData> inflated;
  if (inflated_name != "")
    inflated = ReadPolyData(inflated_name.c_str());


  MeshToSphere m2s;
  m2s.Input(input);
  if (inflated_name != "")
    m2s.Inflated(inflated);
  m2s.Parameters(parser.GetParameters());
  m2s.Debug(debug);
  m2s.Run();
  vtkSmartPointer<vtkPolyData> output = m2s.Output();

  return WritePolyData(POSARG(2), output, fopt);
}
