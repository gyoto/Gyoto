#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"
#include "Python.h"

using namespace Gyoto;
using namespace Gyoto::Spectrum;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Spectrum::Python,
      "Python-based Spectrum class")
GYOTO_PROPERTY_STRING(Python, Module, module,
      "Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Python, Class, klass,
      "Python class (in Module) implementing the Spectrum.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Spectrum::Python, Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_END(Gyoto::Spectrum::Python,
		   Gyoto::Spectrum::Generic::properties)

Spectrum::Python::Python()
: Generic("Python"),
  module_(""), class_(""),
  pModule_(NULL), pClass_(NULL), pInstance_(NULL), pCall_(NULL), pIntegrate_(NULL),
  parameters_()
{}

Python::Python(const Python&o)
  : Generic(o), module_(o.module_), class_(o.class_),
    pModule_(o.pModule_), pClass_(o.pClass_), pInstance_(o.pInstance_),
    pCall_(o.pCall_), pIntegrate_(o.pIntegrate_),
    parameters_(o.parameters_)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XINCREF(pModule_);
  Py_XINCREF(pClass_);
  Py_XINCREF(pInstance_);
  Py_XINCREF(pCall_);
  Py_XINCREF(pIntegrate_);
  PyGILState_Release(gstate);
}

Spectrum::Python::~Python(){
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pIntegrate_);
  Py_XDECREF(pCall_);
  Py_XDECREF(pInstance_);
  Py_XDECREF(pClass_);
  Py_XDECREF(pModule_);
  PyGILState_Release(gstate);
}

Spectrum::Python* Spectrum::Python::clone() const {return new Python(*this);}

std::string Python::module() const { return module_; }
void Python::module(const std::string &m) {
  GYOTO_DEBUG << "Loading Python module " << m << endl;
  PyGILState_STATE gstate;
  module_=m;

  gstate = PyGILState_Ensure();
  PyObject *pName=PyUnicode_FromString(m.c_str());
  if (!pName) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Failed translating string to Python");
  }
  Py_XDECREF(pModule_);
  pModule_ = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Failed loading Python module");
  }
  PyGILState_Release(gstate);
  if (class_ != "") klass(class_);
  GYOTO_DEBUG << "Done loading Python module " << m << endl;
}

std::string Python::klass() const { return class_; }
void Python::klass(const std::string &f) {
  GYOTO_DEBUG << "Instanciating Python class " << f << endl;
  PyGILState_STATE gstate;
  class_=f;
  if (!pModule_) return;
  gstate = PyGILState_Ensure();
  Py_XDECREF(pIntegrate_); pIntegrate_=NULL;
  Py_XDECREF(pCall_); pCall_=NULL;
  Py_XDECREF(pInstance_); pInstance_=NULL;
  Py_XDECREF(pClass_); pClass_=NULL;
  pClass_ = PyObject_GetAttrString(pModule_, class_.c_str());
  if (!pClass_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Could not find class in module");
  }
  if (!PyCallable_Check(pClass_)) {
    Py_DECREF(pClass_);
    pClass_ = NULL;
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Class is not callable");
  }
  pInstance_ = PyObject_CallObject(pClass_, NULL);
  if (!pInstance_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Failed instanciating Python class");
  }

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error instanciating python class");
  }

  pCall_=PyObject_GetAttrString(pInstance_, "__call__");
  if (PyErr_Occurred() || !pCall_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("This class does not seem to implement __call__");
  }

  if (!PyCallable_Check(pCall_)) {
    throwError("Member \"__call__\" present but not callable\n");
  }

  PyObject * pName = PyUnicode_FromString("integrate");
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
  }

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Error looking for integrate method");
  }

  PyGILState_Release(gstate);
  if (parameters_.size()) parameters(parameters_);
  GYOTO_DEBUG << "Done instanciating Python class " << f << endl;
}

std::vector<double> Python::parameters() const {return parameters_;}
void Python::parameters(const std::vector<double> &p){
  parameters_=p;
  if (!pInstance_ || p.size()==0) return;

  PyGILState_STATE gstate = PyGILState_Ensure();

  for (size_t i=0; i<p.size(); ++i) {
    Py_XDECREF(PyObject_CallMethod(pInstance_,
				   "__setitem__",
				   "id", i, p[i]));
    if (PyErr_Occurred()) {
      PyErr_Print();
      PyGILState_Release(gstate);
      throwError("Failed constructing args");
    }

  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "done.\n";
}

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
