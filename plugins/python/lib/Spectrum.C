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
      "Name of Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Spectrum::Python, InlineModule, inlineModule,
      "Inline code of Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Spectrum::Python, Class, klass,
      "Python class (in Module) implementing the Spectrum.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Gyoto::Spectrum::Python, Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_END(Gyoto::Spectrum::Python,
		   Gyoto::Spectrum::Generic::properties)

GYOTO_PROPERTY_THREAD_UNSAFE(Gyoto::Spectrum::Python)

Spectrum::Python::Python()
: Python::Object<Generic>(),
  pCall_(NULL), pIntegrate_(NULL), pCall_overloaded_(false)
{
  kind("Python");
}

Spectrum::Python::Python(const Python&o)
  : Python::Object<Generic>(o),
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
std::string Spectrum::Python::inlineModule() const {return Python::Base::inlineModule();}
void Spectrum::Python::inlineModule(const std::string& m){Python::Base::inlineModule(m);}
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

  pCall_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "__call__");
  pIntegrate_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "integrate");

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error while retrieving methods");
  }

  if (!pCall_) {
    PyGILState_Release(gstate);
    GYOTO_ERROR("Object does not implement required method \"__call__\"");
  }

  pCall_overloaded_ = Gyoto::Python::PyCallable_HasVarArg(pCall_);

  Gyoto::Python::PyInstance_SetThis(pInstance_,
				    Gyoto::Python::pGyotoSpectrum(),
				    this);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error while setting this");
  }

  PyGILState_Release(gstate);
  if (parameters_.size()) parameters(parameters_);
  GYOTO_DEBUG << "Done checking Python class methods" << f << endl;
}

/* Gyoto::Spectrum::Generic API */

double Spectrum::Python::operator()(double nu) const {
  if (!pCall_) GYOTO_ERROR("Python class not loaded yet");
  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  PyObject * pArgs = Py_BuildValue("(d)", nu);
  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    Py_XDECREF(pArgs);
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed building argument list");
  }

  PyObject * pValue = PyObject_CallObject(pCall_, pArgs);
  Py_DECREF(pArgs);
  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    Py_XDECREF(pValue);
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed calling Python method __call__");
  }

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error interpreting result as double");
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
    GYOTO_ERROR("Failed building argument list");
  }

  PyObject * pValue = PyObject_CallObject(pCall_, pArgs);
  Py_DECREF(pArgs);
  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    Py_XDECREF(pValue);
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed calling Python method __call__");
  }

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error interpreting result as double");
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
    GYOTO_ERROR("Failed building argument list");
  }

  PyObject * pValue = PyObject_CallObject(pIntegrate_, pArgs);
  Py_DECREF(pArgs);
  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    Py_XDECREF(pValue);
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed calling Python method integrate");
  }

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error interpreting result as double");
  }

  PyGILState_Release(gstate);

  return res;
}
