
#ifndef M2SPARAMETERS_H
#define M2SPARAMETERS_H

#include "mirtk/Object.h" // for attribute macros
#include "mirtk/EdgeTable.h"

//vtk headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <map>

namespace mirtk {

class AbstractParameter {
  public:
    string str = "Unknown Parameter";
    virtual void setValue(string) = 0;
    virtual string getValue() { return ""; }

    virtual void setValue(void *) = 0;
    virtual void *getValuePtr() = 0;

    virtual string toString() { return ""; }
    virtual ~AbstractParameter() {};
};

template<class T>
class Parameter : public AbstractParameter{

  public:
    T value;
    Parameter(){};
    ~Parameter() {};
    Parameter(T v) { value = v; };
    Parameter(string s, T v) { str = s; value = v;}


    virtual void setValue(string);
    virtual string getValue() { return std::to_string(value); };

    virtual void setValue(void * v) { value = *(T *)v; };
    virtual void *getValuePtr() { return (void*)&value; };

    virtual string toString() { return str + " = " + std::to_string(value); };
    string formatFloatString(string value);

};

struct ParameterSet {
  public:
    bool active = false;
    string name = "Unknown Parameter Set";
    Array<AbstractParameter *> parameters;

    ParameterSet(){
      UpdateLinks();
    };

    virtual ~ParameterSet() {};

    ParameterSet(const ParameterSet &other){
      UpdateLinks();
      CopyAttributes(other);
      UpdateLinks();
    };

    ParameterSet &operator =(const ParameterSet &other){
      if (this != &other) {
        UpdateLinks();
        CopyAttributes(other);
        UpdateLinks();
      }
      return *this;
    }


    virtual void UpdateLinks() {};

    void CopyAttributes(const ParameterSet &other){
      name = other.name;
      active = other.active;
      for (int i=0; i<parameters.size(); i++){
        parameters[i]->str = other.parameters[i]->str;
        void * val = other.parameters[i]->getValuePtr();
        parameters[i]->setValue(val);
      }
    };
};

struct RemeshingParameterSet : public ParameterSet{
   Parameter<double> minEdgeLength =
       Parameter<double>("Minimum Edge Length", -1.0);

   Parameter<double> maxEdgeLength =
       Parameter<double>("Maximum Edge Length", -1.0);

   RemeshingParameterSet(){
     name = "Remeshing";
     UpdateLinks();
   };

   void UpdateLinks(){
     parameters.clear();
     parameters.push_back(&minEdgeLength);
     parameters.push_back(&maxEdgeLength);
   }


 };

 struct SMDSParameterSet : public ParameterSet{
   Parameter<double> distThreshold =
       Parameter<double>("Distance Threshold", -1.0);

   Parameter<double> epsilon =
       Parameter<double>("Epsilon", 0.0);

   Parameter<int> maxIterations =
       Parameter<int>("Maximum Iterations", 100);

   SMDSParameterSet(){
     name = "SMDS";
     UpdateLinks();
   };

   void UpdateLinks(){
     parameters.clear();
     parameters.push_back(&distThreshold);
     parameters.push_back(&epsilon);
     parameters.push_back(&maxIterations);
   }

 };

 struct GeodesicInterpolationParameterSet : public ParameterSet{

   Parameter<double> lambda =
       Parameter<double>("Lambda", 0.01);

   Parameter<int> k =
       Parameter<int>("K", 100);

   Parameter<int> iterations =
       Parameter<int>("Iterations", 1);

   GeodesicInterpolationParameterSet(){
     name = "Geodesic Initialisation";
     UpdateLinks();
   };

   void UpdateLinks(){
     parameters.clear();
     parameters.push_back(&lambda);
     parameters.push_back(&k);
     parameters.push_back(&iterations);
   }
 };

 struct PostProcessingParameterSet : public ParameterSet{

   Parameter<int> iterations =
       Parameter<int>("Iterations", 10);

   Parameter<int> diffusionIterations =
       Parameter<int>("Diffusion Iterations", 1000);

   Parameter<double> percentile =
       Parameter<double>("Percentile", 99);

   Parameter<double> threshold =
       Parameter<double>("Threshold", -1.0);


   PostProcessingParameterSet(){
     name = "Post Processing";
     UpdateLinks();
   };

   void UpdateLinks(){
     parameters.clear();
     parameters.push_back(& iterations);
     parameters.push_back(& diffusionIterations);
     parameters.push_back(& percentile);
     parameters.push_back(& threshold);
   }

 };

 struct DiffusionParameterSet : public ParameterSet{

   Parameter<bool> adaptive =
       Parameter<bool>("Adaptive", 1);

   Parameter<int> maxIterations =
       Parameter<int>("Maximum Iterations", 1000);

   DiffusionParameterSet(){
     name = "Diffusion Initialisation";
     UpdateLinks();
   };

   void UpdateLinks(){
     parameters.clear();
     parameters.push_back(& adaptive);
     parameters.push_back(& maxIterations);
   }

 };



 struct LevelParameters{
   RemeshingParameterSet remesh;
   GeodesicInterpolationParameterSet interp;
   DiffusionParameterSet diffusion;
   PostProcessingParameterSet postprocess;
   SMDSParameterSet smds;


   Array<ParameterSet *> parameterSets;

   LevelParameters(){
     UpdateLinks();
   };

   LevelParameters(const LevelParameters &other){
     UpdateLinks();
     CopyAttributes(other);
     UpdateLinks();
   };

   LevelParameters &operator =(const LevelParameters &other)
   {
     if (this != &other) {
       UpdateLinks();
       CopyAttributes(other);
       UpdateLinks();
     }
     return *this;
   }

   void CopyAttributes(const LevelParameters &other){
     for (int i=0; i<parameterSets.size(); i++){
       *parameterSets[i] = *other.parameterSets[i];
     }
   }

   void UpdateLinks(){
     parameterSets.clear();
     parameterSets.push_back(& remesh);
     parameterSets.push_back(& interp);
     parameterSets.push_back(& diffusion);
     parameterSets.push_back(& postprocess);
     parameterSets.push_back(& smds);
   }


 };


class ParameterParser {

	// ---------------------------------------------------------------------------
	// Attributes

  private:
    ParameterSet * _CurrentParameterSet;
    Array<LevelParameters> _ParameterHierarchy;

	 // Construction/Destruction
	public:

	  /// Constructor
    ParameterParser();

	  /// Copy constructor
    ParameterParser(const ParameterParser &);

	  /// Assignment operator
    ParameterParser &operator =(const ParameterParser &);

	  /// Copy attributes of this class from another instance
	  void CopyAttributes(const ParameterParser &);

	  /// Destructor
	  virtual ~ParameterParser();


	  //Methods

    void ParseFile(string);

    void Initialize(vtkPolyData *);

    string ToString();

    Array<LevelParameters> GetParameters();


	private:
    bool ParseLine(string);
    bool ParseLevel(string);
    bool ParseParameterSet(string);
    bool ParseAttribute(string);

    string SingleSpace(string str);
    string RemoveLeadingWhiteSpace(string str);
    string RemoveTrailingWhiteSpace(string str);
    string RemoveLeadingAndTrailingWhiteSpace(string str);
    string ToLowerCase(string str);

}; // end class

} // namespace mirtk

#endif // M2SPARAMETERS_H
