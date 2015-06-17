#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"
#include "Python.h"

using namespace Gyoto;
using namespace Gyoto::Spectrum;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Spectrum::Python)
GYOTO_PROPERTY_STRING(Python, Module, module)
GYOTO_PROPERTY_STRING(Python, Function, function)
GYOTO_PROPERTY_VECTOR_DOUBLE(Spectrum::Python, parameters, parameters,
			     "List of parameters")
GYOTO_PROPERTY_END(Gyoto::Spectrum::Python,
		   Gyoto::Spectrum::Generic::properties)

Spectrum::Python::Python()
: Generic("Python"),
  module_(""), pModule_(NULL),
  function_(""), pFunc_(NULL),
  parameters_(), expression_("")
{}

Spectrum::Python::~Python(){
  Py_XDECREF(pFunc_);
  Py_XDECREF(pModule_);
}

Spectrum::Python* Spectrum::Python::clone() const {return new Python(*this);}

std::string Python::module() const { return module_; }
void Python::module(const std::string &m) {
  module_=m;
  PyObject *pName=PyString_FromString(m.c_str());
  if (!pName) throwError("Failed translating string to Python");
  Py_XDECREF(pModule_);
  pModule_ = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule_) throwError("Failed loading Python module");
  if (function_ != "") function(function_);
}

std::string Python::function() const { return function_; }
void Python::function(const std::string &f) {
  function_=f;
  if (!pModule_) return;
  Py_XDECREF(pFunc_);
  pFunc_ = PyObject_GetAttrString(pModule_, function_.c_str());
  if (!pFunc_) throwError("Could not find function in module");
  if (!PyCallable_Check(pFunc_)) {
    Py_DECREF(pFunc_);
    pFunc_ = NULL;
    throwError("Function is not callable");
  }
}

std::vector<double> Spectrum::Python::parameters() const {return parameters_;}
void Spectrum::Python::parameters(const std::vector<double> &p){parameters_=p;}

double Spectrum::Python::operator()(double nu) const {
  if (!pFunc_) throwError("Python function not loaded yet");
  PyObject * pArgs  = PyTuple_New(1);
  PyObject * pValue = PyFloat_FromDouble(nu);
  /* pValue reference stolen here */
  PyTuple_SetItem(pArgs, 0, pValue);
  pValue = PyObject_CallObject(pFunc_, pArgs);
  Py_DECREF(pArgs);
  if (!pValue) {
    PyErr_Print();
    throwError("Python function call failed");
  }
  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  return res;
}
