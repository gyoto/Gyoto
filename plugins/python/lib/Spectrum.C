#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"
#include "Python.h"

using namespace Gyoto;
using namespace Gyoto::Spectrum;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Spectrum::Python)
GYOTO_PROPERTY_STRING(Python, Module, module)
GYOTO_PROPERTY_STRING(Python, Class, klass)
GYOTO_PROPERTY_VECTOR_DOUBLE(Spectrum::Python, Parameters, parameters,
			     "List of parameters")
GYOTO_PROPERTY_END(Gyoto::Spectrum::Python,
		   Gyoto::Spectrum::Generic::properties)

Spectrum::Python::Python()
: Generic("Python"),
  module_(""), class_(""),
  pModule_(NULL), pClass_(NULL), pInstance_(NULL),
  parameters_()
{}

Spectrum::Python::~Python(){
  Py_XDECREF(pInstance_);
  Py_XDECREF(pClass_);
  Py_XDECREF(pModule_);
}

Spectrum::Python* Spectrum::Python::clone() const {return new Python(*this);}

std::string Python::module() const { return module_; }
void Python::module(const std::string &m) {
  module_=m;
  PyObject *pName=PyString_FromString(m.c_str());
  if (!pName) {
    PyErr_Print();
    throwError("Failed translating string to Python");
  }
  Py_XDECREF(pModule_);
  pModule_ = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule_) {
    PyErr_Print();
    throwError("Failed loading Python module");
  }
  if (class_ != "") klass(class_);
}

std::string Python::klass() const { return class_; }
void Python::klass(const std::string &f) {
  class_=f;
  if (!pModule_) return;
  Py_XDECREF(pInstance_); pInstance_=NULL;
  Py_XDECREF(pClass_); pClass_=NULL;
  pClass_ = PyObject_GetAttrString(pModule_, class_.c_str());
  if (!pClass_) {
    PyErr_Print();
    throwError("Could not find class in module");
  }
  if (!PyCallable_Check(pClass_)) {
    Py_DECREF(pClass_);
    pClass_ = NULL;
    PyErr_Print();
    throwError("Class is not callable");
  }
  pInstance_ = PyObject_CallObject(pClass_, NULL);
  if (!pInstance_) {
    PyErr_Print();
    throwError("Failed instanciating Python class");
  }
  if (parameters_.size()) parameters(parameters_);
}

std::vector<double> Python::parameters() const {return parameters_;}
void Python::parameters(const std::vector<double> &p){
  parameters_=p;
  if (!pInstance_) return;
  PyObject * pArgs  = PyTuple_New(1);
  PyObject * pValue;
  for (size_t i=0; i<parameters_.size(); ++i) {
    pValue = PyFloat_FromDouble(parameters_[i]);
    if (!pValue) throwError("Failed converting parameter to Python Value");
    /* pValue reference stolen here */
    PyTuple_SetItem(pArgs, i, pValue);
  }
  pValue = PyObject_CallMethod(pInstance_, "setParameters", "O", pArgs);
  Py_DECREF(pArgs);
  if (!pValue) {
    PyErr_Print();
    throwError("Failed calling Python method setParameters");
  }
  Py_DECREF(pValue);
}

double Spectrum::Python::operator()(double nu) const {
  if (!pInstance_) throwError("Python class not loaded yet");
  PyObject * pValue =
    PyObject_CallMethod(pInstance_, "__call__", "(d)", nu);
  if (!pValue) {
    PyErr_Print();
    throwError("Python class call failed");
  }
  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  return res;
}
