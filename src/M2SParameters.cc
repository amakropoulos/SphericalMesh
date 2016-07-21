

//mirtk includes
#include "mirtk/M2SParameters.h"

#include "mirtk/Common.h"

//my mirtk includes
#include "mirtk/EdgeTable.h"
#include <boost/regex.hpp>
//vtk includes

namespace mirtk {
// =============================================================================
// Construction/Destruction
// =============================================================================

template<>
void Parameter<bool>::setValue(string valueString) {
  Array<string> trueStrings = { "True", "true", "1", "On", "on" };
  Array<string> falseStrings = { "False", "false", "0", "Off", "off" };

  for (auto s: trueStrings){
    if (valueString == s){ value = true;  return; }
  }

  for (auto s: falseStrings){
    if (valueString == s){ value = false;  return; }
  }

  cout << "Unrecognised bool value " << valueString << endl;
  exit(1);
}


template<>
void Parameter<int>::setValue(string valueString) {
  value = atoi(valueString.c_str());
}

template<>
void Parameter<double>::setValue(string valueString) {
  value = atof(valueString.c_str());
}

template<>
void Parameter<float>::setValue(string valueString) {
  value = atof(valueString.c_str());
}
template<>
void Parameter<string>::setValue(string valueString) {
  value = valueString;
}

template<>
string Parameter<float>::getValue() {
  string valStr = std::to_string(value);
  return formatFloatString(valStr);
}

template<>
string Parameter<double>::getValue() {
  string valStr = std::to_string(value);
  return formatFloatString(valStr);
}



template<class T>
string Parameter<T>::formatFloatString(string str){
  int strEnd = str.find_last_not_of('0');
  if (strEnd != std::string::npos){
    int strRange = strEnd + 1;
    str = str.substr(0, strRange);
  }

  int decimalPtId = str.find('.');
  if (decimalPtId == str.size() - 1){
    str = str.substr(0, decimalPtId);
  }
  return str;
}


ParameterParser::ParameterParser(){
  _CurrentParameterSet = NULL;

	//string preprocessing
	//remove leading/trailing whitespace
	//and replace consecutive whitespace chars with single spaces
	/*
	for (auto &e: attributeStringMap){
		e.second = SingleSpace(e.second);
		e.second = RemoveLeadingAndTrailingWhiteSpace(e.second);
	}

	for (auto &e: categoryStringMap){
		e.second = SingleSpace(e.second);
		e.second = RemoveLeadingAndTrailingWhiteSpace(e.second);
	}
	*/

}

// -----------------------------------------------------------------------------


ParameterParser::ParameterParser(const ParameterParser &other)
{
  *this = other;
}

// -----------------------------------------------------------------------------

ParameterParser &ParameterParser::operator =(const ParameterParser &other)
{
  if (this != &other) {
    ParameterParser::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ParameterParser::~ParameterParser(){
}

void ParameterParser::CopyAttributes(const ParameterParser &other){
  cout << "TODO copy attributes" << endl;
}


// =============================================================================
// Class methods
// =============================================================================

const string whitespace = " _\t";
const Array<string> commentStrings = {"//", "#"};


Array<string> &split(const string &s, char delim, Array<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }

string ParameterParser::SingleSpace(string str) {
  std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
  str.erase(new_end, str.end());
  return str;
}

string RemoveChars(string str, string chars){
   for (auto c:chars){
      str.erase (std::remove(str.begin(), str.end(), c), str.end());
   }
   return str;
}


string ParameterParser::RemoveLeadingAndTrailingWhiteSpace(string str){
  str = RemoveTrailingWhiteSpace(str);
  return RemoveLeadingWhiteSpace(str);
}

string ParameterParser::RemoveTrailingWhiteSpace(string str){
  int strEnd = str.find_last_not_of(whitespace);

  if (strEnd != std::string::npos){
    int strRange = strEnd + 1;
    return str.substr(0, strRange);
  }
  else
    return "";
}

string ParameterParser::RemoveLeadingWhiteSpace(string str){

  int strBegin = str.find_first_not_of(whitespace);

  if (strBegin != std::string::npos){
    int strRange = str.size() - strBegin + 1;
    return str.substr(strBegin, strRange);
  }
  else
    return "";

}

string ParameterParser::ToLowerCase(string str){
	string strCopy = str;
    std::transform(strCopy.begin(), strCopy.end(), strCopy.begin(), ::tolower);
    return strCopy;
}




Array<string> split(const string &s, char delim) {
    Array<string>  elems;
    split(s, delim, elems);
    return elems;
}

string RemoveComment(string line){

  for (auto commentStr : commentStrings){

    int commentBegin = line.find_first_of(commentStr,0);
    if (commentBegin !=std::string::npos)
      line = line.substr(0, commentBegin);
  }

    return line;
}


bool ParameterParser::ParseLevel(string str){
  //cout << "parsing level" << endl;

  Array<string> tokens = split(str,' ');
  if (tokens.size() != 2){ return false;}
  if (ToLowerCase(tokens[0]) != "level") { return false;}

  //cout << "I've found a level statement!" << endl;

  string levelNumStr = tokens[1];
  for (auto c: levelNumStr){
    if (! isdigit(c)) { return false;}
  }

  //cout << "and it has a number! " << levelNumStr << endl;


  LevelParameters p;
  _ParameterHierarchy.push_back(p);

  int level = atoi(levelNumStr.c_str());

  if (level != _ParameterHierarchy.size()){
    cout << "M2SParameters::ParseLevel invalid level number." << endl;
    cout << str << endl;
    cout << "Expected Level " << _ParameterHierarchy.size() << "!" << endl;
    cout << level << endl;
    exit(1);
  }

  //cout << "so level statement parsed, winning" << endl;
  _CurrentParameterSet = NULL;

  return true;
}

bool ParameterParser::ParseParameterSet(string str){
  //cout << "Parsing parameter set"<< endl;

  if (_ParameterHierarchy.empty()){
    cout << "M2SParameters::ParseCategory" << endl;
    cout << str << endl;
    cout << "Expected level before category !" << endl;
    exit(1);
  }

  bool active = true;
  Array<string> tokens = split(str,' ');
  if (tokens.back() == "on"){
    str = str.substr(0, str.size() - 3); // remove "on" and preceding space
  }
  else if (tokens.back() == "off"){
    active = false;
    str = str.substr(0, str.size() - 4); // remove "off" and preceding space
  }

  Array<ParameterSet *> parSets = _ParameterHierarchy.back().parameterSets;

  for (auto parameterSet: parSets ){
    if (str == ToLowerCase(parameterSet->name)){
      parameterSet->active = active;
      _CurrentParameterSet = parameterSet;
      //cout << "Found parameter set " << parameterSet->name << endl;
      return true;
    }
  }


  return false;
}


bool ParameterParser::ParseAttribute(string str){

  if (_CurrentParameterSet == NULL){
    cout << "M2SParameters::ParseAttribute" << endl;
    cout << str << endl;
    cout << "Expected parameterset heading before parameter ! " << endl;
    exit(1);
  }

  Array<string> tokens = split(str,'=');
  if (tokens.size() != 2){ return false;}

  //cout << "I've found something that looks like a parameter!" << endl;


  string name = tokens[0];
  //format string for matching
  name = RemoveLeadingAndTrailingWhiteSpace(name);
  name = SingleSpace(name);
  name = ToLowerCase(name);

  string valStr = RemoveLeadingAndTrailingWhiteSpace(tokens[1]).c_str();


  for (auto &p: _CurrentParameterSet->parameters){
    string n = p->str;
    n = RemoveLeadingAndTrailingWhiteSpace(n);
    n = SingleSpace(n);
    n = ToLowerCase(n);
    if (name == n){
      p->setValue(valStr);
      return true;
    }

  }


  return false;
}


bool ParameterParser::ParseLine(string line){
  string str = RemoveComment(line);
  str = RemoveChars(str,"[]_");

  str = RemoveLeadingAndTrailingWhiteSpace(str);
  str = ToLowerCase(str);
  str = SingleSpace(str);

  if (str == ""){ return true; } //line is whitespace
  if (ParseLevel(str)) { return true; }
  if (ParseParameterSet(str)){ return true; }
  return ParseAttribute(str);

/*
  if (str == ""){ return; } //line is whitespace






  int attribute = ParseAttribute(attributeStr);
  const char * val = RemoveLeadingAndTrailingWhiteSpace(tokens[1]).c_str();



  switch (attribute){


     case Attribute::UseMeshEdgeLinks:
       p->UseMeshEdgeLinks = atoi(val);
       break;

  }*/
}

void ParameterParser::ParseFile(string fileName){
  ifstream infile(fileName);

  string line;
  while (std::getline(infile, line)){
    bool success = ParseLine(line);




    if (!success) {
      cout << "failed to parse line" << line << endl;
      exit(1);
    }
  }

}
/*
void M2SParameters::DefaultInitialize(vtkPolyData * poly){

  double maxEdgeLength = 4;
  double minEdgeLength = 1;


  int numPoints = poly->GetNumberOfPoints();
  double meanEdgeLength = 0.0;

  EdgeTable edges(poly);
  EdgeIterator ei(EdgeTable);
  ei.InitTraversal();
  int pid1, pid2;
  while ( ei.GetNextEdge(pid1, pid2) != -1 ){
    double pt1[3], pt2[3];

    poly->GetPoint(pid1, pt1);
    poly->GetPoint(pid2, pt2);

    double d = 0;
    for (int i=0; i<3; i++){ d += pow(pt2[i] - pt1[i], 2); }
    meanEdgeLength += sqrt(d);
  }

  meanEdgeLength /= edges.NumberOfEdges();

  double minEdgeLengthByMesh = pow(2, ceil(log2(meanEdgeLength)));
  if (minEdgeLengthByMesh >  minEdgeLength) { minEdgeLength = minEdgeLengthByMesh; }

  // 0.5 //
  // 8.8 // 8

  if (maxEdgeLength > minEdgeLengthByMesh){
    LevelParameters p;
    p.smds.maxIterations = 100;
    _ParameterHierarchy.push_back(p);
  }

  for (double length = maxEdgeLength; length >= minEdgeLength; length/=2){

    LevelParameters p;
    p.smds.maxIterations = (int) 25 * length;
    p.interp.k = 100;
    p.interp.lambda = 0.01;
    _ParameterHierarchy.push_back(p);
    _ParameterHierarchy.push_back(p);
  }

  double edgeLength = 2;

  LevelParameters p;
  p.remesh.maxEdgeLength = -1;
  p.remesh.minEdgeLength = -1;
  p.smds.maxIterations = 0;
  p.interp.k = 100;
  p.interp.lambda = 0.01;
  _ParameterHierarchy.push_back(p);






}
*/



string ParameterParser::ToString(){


  ostringstream os;

  int l = 0;
  for (auto levelParameters: _ParameterHierarchy){
    l++;
    os << "\n";
    os << "Level " << std::to_string(l) << "\n";
    os << "_______";
    int digits = 1 + (l / 10);
    for (int i=0; i<digits; i++){ os << "_"; }
    os << "\n\n\n";



    for (auto parameterSet: levelParameters.parameterSets){
      if (parameterSet->active){
        os << "  [ " << parameterSet->name <<  " ]" << "\n" << "\n";
        for (auto p: parameterSet->parameters){
          os << "    " << p->toString() << "\n";
        }
      }
      else{
        os << "  [ " << parameterSet->name <<  " Off ]" << "\n" << "\n";
      }
      os << "\n" << "\n";
    }

  }

  string str = os.str();

  if (! _ParameterHierarchy.empty())
    str = str.substr(0, str.size() - 2);

  return str;
}

Array<LevelParameters> ParameterParser::GetParameters(){
  return _ParameterHierarchy;
}





} // mirtk namespace




