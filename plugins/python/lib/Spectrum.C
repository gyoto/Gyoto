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
  pModule_(NULL), pClass_(NULL), pInstance_(NULL), pIntegrate_(NULL),
  parameters_()
{}

Python::Python(const Python&o)
  : Generic(o), module_(o.module_), class_(o.class_),
    pModule_(o.pModule_), pClass_(o.pClass_), pInstance_(o.pInstance_), pIntegrate_(o.pIntegrate_),
    parameters_(o.parameters_)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XINCREF(pModule_);
  Py_XINCREF(pClass_);
  Py_XINCREF(pInstance_);
  Py_XINCREF(pIntegrate_);
  PyGILState_Release(gstate);
}

Spectrum::Python::~Python(){
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pIntegrate_);
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
  GYOTO_DEBUG << "Setting parameters to vector of size " << p.size() << endl;
  parameters_=p;
  if (!pInstance_) return;

  /* Multi-threading guard */
  PyGILState_STATE gstate = PyGILState_Ensure();

  /* Find setPArameter method in insance. The object returns knows 'self'. */
  GYOTO_DEBUG << "Looking for method \"setParameters\" in instance" << endl;
  PyObject * pMethod =  PyObject_GetAttrString(pInstance_, "setParameters");
  if (!pMethod) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Could not find method in class");
  }
  if (!PyCallable_Check(pMethod)) {
    Py_DECREF(pMethod);
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Method is not callable");
  }

  /* Transform parameter vector into Python tuple */
  PyObject * pArgs  = PyTuple_New(p.size());
  PyObject * pValue;
  for (size_t i=0; i<parameters_.size(); ++i) {
    GYOTO_DEBUG << "Adding parameter #" << i << " with value " << p[i] << " into parameter list\n";
    pValue = PyFloat_FromDouble(parameters_[i]);
    if (!pValue) throwError("Failed converting parameter to Python Value");
    /* pValue reference stolen here */
    PyTuple_SetItem(pArgs, i, pValue);
    if (PyErr_Occurred()) {
      PyErr_Print();
      Py_DECREF(pArgs);
      PyGILState_Release(gstate);
      throwError("Failed constructing args");
    }
  }

  /* Actually call method */
  GYOTO_DEBUG << "Calling method setParameters\n";
  pValue = PyObject_CallObject(pMethod, pArgs);

  /* Check for error, release memory, release lock */
  GYOTO_DEBUG << "Doing house-keeping tasks\n";
  Py_DECREF(pArgs);
  Py_DECREF(pMethod);
  Py_XDECREF(pValue);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Failed calling Python method setParameters");
  }
  PyGILState_Release(gstate);
  GYOTO_DEBUG << "done.\n";
}

double Spectrum::Python::operator()(double nu) const {
  if (!pInstance_) throwError("Python class not loaded yet");
  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  PyObject * pValue =
    PyObject_CallMethod(pInstance_, "__call__", "(d)", nu);
  if (!pValue) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Python class call failed");
  }
  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pValue);
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
