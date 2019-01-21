#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL GyotoPython_ARRAY_API
#include <numpy/arrayobject.h>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Metric::Python,
      "Python-based Metric class")
GYOTO_PROPERTY_STRING(Gyoto::Metric::Python, Module, module,
      "Python module containing the Metric implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Metric::Python, InlineModule, inlineModule,
      "Inline code of Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Metric::Python, Class, klass,
      "Python class (in Module) implementing the Metric.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Gyoto::Metric::Python, Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_BOOL(Metric::Python, Spherical, Cartesian, spherical,
      "Whether the coordinate system is Spherical or (default) Cartesian.")
GYOTO_PROPERTY_END(Metric::Python, Generic::properties)

GYOTO_PROPERTY_THREAD_UNSAFE(Metric::Python)

// Birth and death
Gyoto::Metric::Python::Python()
: Generic(GYOTO_COORDKIND_CARTESIAN, "Python"),
  Base(),
  pGmunu_(NULL), pChristoffel_(NULL)
{}

Gyoto::Metric::Python::Python(const Python& o)
: Generic(o),
  Base(o),
  pGmunu_(o.pGmunu_), pChristoffel_(o.pChristoffel_)
{
  Py_XINCREF(pGmunu_);
  Py_XINCREF(pChristoffel_);
}

Gyoto::Metric::Python::~Python() {
  Py_XDECREF(pChristoffel_);
  Py_XDECREF(pGmunu_);
}

Metric::Python* Gyoto::Metric::Python::clone() const {return new Python(*this);}

// Property accessors
void Gyoto::Metric::Python::spherical(bool t) {
  coordKind(t?GYOTO_COORDKIND_SPHERICAL:GYOTO_COORDKIND_CARTESIAN);
  
  if (!pInstance_) return;

  GYOTO_DEBUG << "Set \"spherical\"\n";
  PyGILState_STATE gstate = PyGILState_Ensure();

  int res = PyObject_SetAttrString(pInstance_, "spherical", t?Py_True:Py_False);

  if (PyErr_Occurred() || res == -1) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed setting \"spherical\" using __setattr__");
  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "done.\n";

}

bool Gyoto::Metric::Python::spherical() const {
  if (coordKind() == GYOTO_COORDKIND_UNSPECIFIED)
    GYOTO_ERROR("coordKind unspecified");
  return coordKind() == GYOTO_COORDKIND_SPHERICAL;
}

void Gyoto::Metric::Python::mass(double m) {
  Generic::mass(m);
  
  if (!pInstance_) return;

  GYOTO_DEBUG << "Setting \"mass\"\n";
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject * pM = PyFloat_FromDouble(mass());

  int res = PyObject_SetAttrString(pInstance_, "mass", pM);

  Py_DECREF(pM);

  if (PyErr_Occurred() || res == -1) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed setting \"mass\" using __setattr__");
  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "done.\n";
  
}

std::vector<double> Metric::Python::parameters() const
{return Python::Base::parameters();}
void Metric::Python::parameters(const std::vector<double>& m)
{Python::Base::parameters(m);}
std::string Metric::Python::module() const {return Python::Base::module();}
void Metric::Python::module(const std::string& m){Python::Base::module(m);}
std::string Metric::Python::inlineModule() const {return Python::Base::inlineModule();}
void Metric::Python::inlineModule(const std::string& m){Python::Base::inlineModule(m);}
std::string Metric::Python::klass() const {return Python::Base::klass();}
void Gyoto::Metric::Python::klass(const std::string &f) {

  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pChristoffel_); pChristoffel_=NULL;
  Py_XDECREF(pGmunu_); pGmunu_=NULL;
  PyGILState_Release(gstate);
  
  Python::Base::klass(f);
  if (!pModule_) return;

  gstate = PyGILState_Ensure();
  GYOTO_DEBUG << "Checking Python class methods" << f << endl;

  pGmunu_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "gmunu");
  pChristoffel_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "christoffel");

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error while retrieving methods");
  }

  if (!pGmunu_) {
    PyGILState_Release(gstate);
    GYOTO_ERROR("Object does not implement required method \"__call__\"");
  }

  if (!pChristoffel_) {
    PyGILState_Release(gstate);
    GYOTO_ERROR("Object does not implement required method \"getVelocity\"");
  }

  Gyoto::Python::PyInstance_SetThis(pInstance_,
				    Gyoto::Python::pGyotoMetric(),
				    this);

  PyGILState_Release(gstate);
  if (parameters_.size()) parameters(parameters_);
  if (coordKind()) spherical(spherical());
  mass(mass());
  GYOTO_DEBUG << "Done checking Python class methods" << f << endl;
}

void Metric::Python::gmunu(double g[4][4], const double * x) const {
  if (!pGmunu_) GYOTO_ERROR("gmunu method not loaded yet");
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp g_dims[] = {4, 4};
  
  PyObject * pG = PyArray_SimpleNewFromData(2, g_dims, NPY_DOUBLE, g);
  PyObject * pX = PyArray_SimpleNewFromData(1, g_dims, NPY_DOUBLE, const_cast<double*>(x));
  PyObject * pR = PyObject_CallFunctionObjArgs(pGmunu_, pG, pX, NULL);

  Py_XDECREF(pR);
  Py_XDECREF(pX);
  Py_XDECREF(pG);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::Python::gmunu");
  }
   
  PyGILState_Release(gstate);
}

int Metric::Python::christoffel(double dst[4][4][4], const double * x) const {
  if (!pChristoffel_) GYOTO_ERROR("christoffel method not loaded yet");
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp d_dims[] = {4, 4, 4};
  
  PyObject * pD = PyArray_SimpleNewFromData(3, d_dims, NPY_DOUBLE, dst);
  PyObject * pX = PyArray_SimpleNewFromData(1, d_dims, NPY_DOUBLE, const_cast<double*>(x));
  PyObject * pR = PyObject_CallFunctionObjArgs(pChristoffel_, pD, pX, NULL);

  double r = PyFloat_AsDouble(pR);
  
  Py_XDECREF(pR);
  Py_XDECREF(pX);
  Py_XDECREF(pD);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::Python::gmunu");
  }
   
  PyGILState_Release(gstate);

  return r;
}
