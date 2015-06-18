#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"
#include "Python.h"

using namespace Gyoto;
using namespace Gyoto::Spectrum;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Spectrum::Python,
      "Python-based Spectrum class")
GYOTO_PROPERTY_STRING(Gyoto::Spectrum::Python, Module, module,
      "Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Spectrum::Python, Class, klass,
      "Python class (in Module) implementing the Spectrum.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Gyoto::Spectrum::Python, Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_END(Gyoto::Spectrum::Python,
		   Gyoto::Spectrum::Generic::properties)

Spectrum::Python::Python()
: Generic("Python"), Base(),
  pCall_(NULL), pIntegrate_(NULL), pCall_overloaded_(false)
{}

Spectrum::Python::Python(const Python&o)
  : Generic(o), Base(o),
    pCall_(o.pCall_), pIntegrate_(o.pIntegrate_),
    pCall_overloaded_(o.pCall_overloaded_)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XINCREF(pCall_);
  Py_XINCREF(pIntegrate_);
  PyGILState_Release(gstate);
}

Spectrum::Python::~Python(){
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pIntegrate_);
  Py_XDECREF(pCall_);
  PyGILState_Release(gstate);
}

Spectrum::Python* Spectrum::Python::clone() const {return new Python(*this);}

/* Gyoto::Python::Base API: reimplement ::klass method */

std::vector<double> Spectrum::Python::parameters() const
{return Python::Base::parameters();}
void Spectrum::Python::parameters(const std::vector<double>& m)
{Python::Base::parameters(m);}
std::string Spectrum::Python::module() const {return Python::Base::module();}
void Spectrum::Python::module(const std::string& m){Python::Base::module(m);}
std::string Spectrum::Python::klass() const {return Python::Base::klass();}
void Spectrum::Python::klass(const std::string &f) {

  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pIntegrate_); pIntegrate_=NULL;
  Py_XDECREF(pCall_); pCall_=NULL;
  PyGILState_Release(gstate);

  Python::Base::klass(f);
  if (!pModule_) return;

  gstate = PyGILState_Ensure();
  GYOTO_DEBUG << "Checking Python class methods" << f << endl;

  pCall_=PyObject_GetAttrString(pInstance_, "__call__");
  if (PyErr_Occurred() || !pCall_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("This class does not seem to implement __call__");
  }

  if (!PyCallable_Check(pCall_)) {
    throwError("Member \"__call__\" present but not callable\n");
  }

  PyObject * pName = PyUnicode_FromString("inspect");
  PyObject * pModule = PyImport_Import(pName);
  Py_XDECREF(pName); pName=NULL;
  PyObject * pFunc = PyObject_GetAttrString(pModule, "getargspec");
  PyObject * pArgSpec = PyObject_CallFunctionObjArgs(pFunc, pCall_, NULL);
  Py_XDECREF(pFunc);
  if (PyErr_Occurred() || !pArgSpec) {
    PyErr_Print();
    Py_XDECREF(pArgSpec);
    PyGILState_Release(gstate);
    throwError("Error checking __call__ arguments");
  }

  pCall_overloaded_ = (PyTuple_GetItem(pArgSpec, 1) != Py_None);
  Py_XDECREF(pArgSpec);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error preparing Python string");
  }

  pName = PyUnicode_FromString("integrate");
  if (PyErr_Occurred() || !pName) {
    PyErr_Print();
    Py_XDECREF(pName);
    PyGILState_Release(gstate);
    throwError("Error preparing Python string");
  }

  if (PyObject_HasAttr(pInstance_, pName)) {
    pIntegrate_ = PyObject_GetAttr(pInstance_, pName);
    Py_DECREF(pName);
    if (!PyCallable_Check(pIntegrate_)) {
      GYOTO_WARNING << "Member \"integrate\" present but not callable\n";
      Py_DECREF(pIntegrate_);
    }
  } else Py_XDECREF(pName);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error looking for integrate method");
  }

  PyGILState_Release(gstate);
  if (parameters_.size()) parameters(parameters_);
  GYOTO_DEBUG << "Done checking Python class methods" << f << endl;
}

/* Gyoto::Spectrum::Generic API */

double Spectrum::Python::operator()(double nu) const {
  if (!pCall_) throwError("Python class not loaded yet");
  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  PyObject * pArgs = Py_BuildValue("(d)", nu);
  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    Py_XDECREF(pArgs);
    PyGILState_Release(gstate);
    throwError("Failed building argument list");
  }

  PyObject * pValue = PyObject_CallObject(pCall_, pArgs);
  Py_DECREF(pArgs);
  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    Py_XDECREF(pValue);
    PyGILState_Release(gstate);
    throwError("Failed calling Python method __call__");
  }

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error interpreting result as double");
  }

  PyGILState_Release(gstate);

  return res;
}

double Spectrum::Python::operator()(double nu, double opacity, double ds) const {
  if (!pCall_overloaded_) return Generic::operator()(nu, opacity, ds);

  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  PyObject * pArgs = Py_BuildValue("(ddd)", nu, opacity, ds);
  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    Py_XDECREF(pArgs);
    PyGILState_Release(gstate);
    throwError("Failed building argument list");
  }

  PyObject * pValue = PyObject_CallObject(pCall_, pArgs);
  Py_DECREF(pArgs);
  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    Py_XDECREF(pValue);
    PyGILState_Release(gstate);
    throwError("Failed calling Python method __call__");
  }

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error interpreting result as double");
  }

  PyGILState_Release(gstate);

  return res;

}

double Spectrum::Python::integrate(double nu1, double nu2) {
  if (!pIntegrate_) return Generic::integrate(nu1, nu2);

  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject * pArgs = Py_BuildValue("dd", nu1, nu2);
  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    Py_XDECREF(pArgs);
    PyGILState_Release(gstate);
    throwError("Failed building argument list");
  }

  PyObject * pValue = PyObject_CallObject(pIntegrate_, pArgs);
  Py_DECREF(pArgs);
  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    Py_XDECREF(pValue);
    PyGILState_Release(gstate);
    throwError("Failed calling Python method integrate");
  }

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error interpreting result as double");
  }

  PyGILState_Release(gstate);

  return res;
}
