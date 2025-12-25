#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"
#include "GyotoValue.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"
#include "GyotoFactoryMessenger.h"

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL GyotoPython_ARRAY_API
#include <numpy/arrayobject.h>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace Gyoto::Python;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Metric::Python,
      "Python-based Metric class")
GYOTO_PROPERTY_STRING(Gyoto::Metric::Python, Module, module,
      "Python module containing the Metric implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Metric::Python, InlineModule, inlineModule,
      "Inline code of Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Metric::Python, Class, klass,
      "Python class (in Module) implementing the Metric.")
GYOTO_PROPERTY_SIZE_T(Gyoto::Metric::Python, Instance, instance,
      "ID of pre-constructed Class instance as per gyoto.core.gyotoid().")
GYOTO_PROPERTY_VECTOR_DOUBLE(Gyoto::Metric::Python, Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_BOOL(Metric::Python, Spherical, Cartesian, spherical,
      "Whether the coordinate system is Spherical or (default) Cartesian.")
GYOTO_PROPERTY_END(Metric::Python, Generic::properties)

GYOTO_PROPERTY_THREAD_UNSAFE(Metric::Python)

// Birth and death
Gyoto::Metric::Python::Python()
: Python::Object<Generic>(),
  pGmunu_(NULL), pChristoffel_(NULL), pGetRmb_(NULL), pGetRms_(NULL),
  pGetSpecificAngularMomentum_(NULL), pGetPotential_(NULL),
  pIsStopCondition_(NULL), pCircularVelocity_(NULL)
{
  kind("Python");
  coordKind(GYOTO_COORDKIND_CARTESIAN);
}

Gyoto::Metric::Python::Python(const Python& o)
  :
  Python::Object<Generic>(o),
  pGmunu_(o.pGmunu_), pChristoffel_(o.pChristoffel_), pGetRmb_(o.pGetRmb_),
  pGetRms_(o.pGetRms_),
  pGetSpecificAngularMomentum_(o.pGetSpecificAngularMomentum_),
  pGetPotential_(o.pGetPotential_),
  pIsStopCondition_(o.pIsStopCondition_),
  pCircularVelocity_(o.pCircularVelocity_)

{
  Py_XINCREF(pGmunu_);
  Py_XINCREF(pChristoffel_);
  Py_XINCREF(pGetRmb_);
  Py_XINCREF(pGetRms_);
  Py_XINCREF(pGetSpecificAngularMomentum_);
  Py_XINCREF(pGetPotential_);
  Py_XINCREF(pIsStopCondition_);
  Py_XINCREF(pCircularVelocity_);
}

Gyoto::Metric::Python::~Python() {
  Py_XDECREF(pCircularVelocity_);
  Py_XDECREF(pIsStopCondition_);
  Py_XDECREF(pGetPotential_);
  Py_XDECREF(pGetSpecificAngularMomentum_);
  Py_XDECREF(pGetRms_);
  Py_XDECREF(pGetRmb_);
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
void Gyoto::Metric::Python::klass(const std::string &c)
{Gyoto::Python::Base::klass(c);}

void Metric::Python::instance(size_t i)
{Gyoto::Python::Base::instance(reinterpret_cast<PyObject*>(i));}
size_t Metric::Python::instance() const
{return reinterpret_cast<size_t>(pInstance_);}

void Metric::Python::detachInstance() {
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  Py_CLEAR(pGetPotential_);
  Py_CLEAR(pGetSpecificAngularMomentum_);
  Py_CLEAR(pGetRms_);
  Py_CLEAR(pGetRmb_);
  Py_CLEAR(pChristoffel_);
  Py_CLEAR(pGmunu_);
  Gyoto::Python::Base::detachInstance();
}

void Gyoto::Metric::Python::attachInstance (PyObject * instance) {
  Gyoto::Python::Base::attachInstance(instance);
  GYOTO_DEBUG << "Checking methodes for Python class " << class_ << endl;

  pGmunu_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "gmunu");
  pChristoffel_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "christoffel");
  pGetRmb_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "getRmb");
  pGetRms_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "getRms");
  pGetSpecificAngularMomentum_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "getSpecificAngularMomentum");
  pGetPotential_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "getPotential");
  pIsStopCondition_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "isStopCondition");
  pCircularVelocity_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "circularVelocity");

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error while retrieving methods");
  }

  if (!pGmunu_) {
    GYOTO_ERROR("Object does not implement required method \"gmunu\"");
  }

  if (!pChristoffel_) {
    GYOTO_ERROR("Object does not implement required method \"christoffel\"");
  }

  GYOTO_DEBUG << "setting 'this' in the instance" << endl;
  Gyoto::Python::PyInstance_SetThis(pInstance_,
				    Gyoto::Python::pGyotoMetric(),
				    this);

  if (parameters_.size()) parameters(parameters_);
  if (coordKind()) spherical(spherical());
  mass(mass());
  GYOTO_DEBUG << "Done checking methods for Python class " << class_ << endl;
}

void Metric::Python::gmunu(double g[4][4], const double * x) const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pGmunu_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pGmunu_);
  }

  GYOTO_DEBUG_EXPR(pGmunu_);
  if (!pGmunu_) GYOTO_ERROR("gmunu method not loaded yet");
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp g_dims[] = {4, 4};

  PyObject * pG = PyArray_SimpleNewFromData(2, g_dims, NPY_DOUBLE, g);
  PyObject * pX = PyArray_SimpleNewFromData(1, g_dims, NPY_DOUBLE, const_cast<double*>(x));

  isRecursing = true;
  PyObject * pR = PyObject_CallFunctionObjArgs(pGmunu_, pG, pX, NULL);
  isRecursing = false;

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
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pChristoffel_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pChristoffel_);
  }

  GYOTO_DEBUG_EXPR(pChristoffel_);
  if (!pChristoffel_) GYOTO_ERROR("christoffel method not loaded yet");
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp d_dims[] = {4, 4, 4};

  PyObject * pD = PyArray_SimpleNewFromData(3, d_dims, NPY_DOUBLE, dst);
  PyObject * pX = PyArray_SimpleNewFromData(1, d_dims, NPY_DOUBLE, const_cast<double*>(x));

  isRecursing=true;
  PyObject * pR = PyObject_CallFunctionObjArgs(pChristoffel_, pD, pX, NULL);
  isRecursing=false;

  Py_XDECREF(pX);
  Py_XDECREF(pD);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::Python::christoffel");
  }

  double r = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);

  PyGILState_Release(gstate);

  return r;
}

double Metric::Python::getRmb
()
  const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pGetRmb_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pGetRmb_);
  }

  GYOTO_DEBUG_EXPR(pGetRmb_);
  if (!pGetRmb_)
    return Metric::Generic::getRmb();

  PyGILState_STATE gstate = PyGILState_Ensure();

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pGetRmb_, NULL);
  isRecursing=false;

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::getRmb()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

double Metric::Python::getRms
()
  const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pGetRms_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pGetRms_);
  }

  GYOTO_DEBUG_EXPR(pGetRms_);
  if (!pGetRms_)
    return Metric::Generic::getRms();

  PyGILState_STATE gstate = PyGILState_Ensure();

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pGetRms_, NULL);
  isRecursing=false;

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::getRms()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

double Metric::Python::getSpecificAngularMomentum
(double rr)
  const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pGetSpecificAngularMomentum_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pGetSpecificAngularMomentum_);
  }

  GYOTO_DEBUG_EXPR(pGetSpecificAngularMomentum_);
  if (!pGetSpecificAngularMomentum_)
    return Metric::Generic::getSpecificAngularMomentum(rr);

  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject * pRr = PyFloat_FromDouble(rr);

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pGetSpecificAngularMomentum_, pRr, NULL);
  isRecursing=false;

  Py_XDECREF(pRr);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::getSpecificAngularMomentum()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

double Metric::Python::getPotential
(double const pos[4], double l_cst)
  const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pGetPotential_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pGetPotential_);
  }

  GYOTO_DEBUG_EXPR(pGetPotential_);
  if (!pGetPotential_)
    return Metric::Generic::getPotential(pos, l_cst);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims_pos[] = {4};

  PyObject * pPo = PyArray_SimpleNewFromData(1, dims_pos, NPY_DOUBLE, const_cast<double*>(pos));
  PyObject * pCs = PyFloat_FromDouble(l_cst);

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pGetPotential_, pPo, pCs, NULL);
  isRecursing=false;

  Py_XDECREF(pCs);
  Py_XDECREF(pPo);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::getPotential()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

int Metric::Python::isStopCondition
(double const coord[8])
  const {
  /* This method, part of the Gyoto::Metric::Generic API, wraps around the
     Python method of same same in a Python class.
   */

  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pIsStopCondition_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pIsStopCondition_);
  }

  GYOTO_DEBUG_EXPR(pIsStopCondition_);
  // If this method was not found is the Python class, then the
  // pIsStopCondition_ pointer is null.
  if (!pIsStopCondition_)
    return Metric::Generic::isStopCondition(coord);

  // All functions need to claim the GIL state.
  PyGILState_STATE gstate = PyGILState_Ensure();

  // Wrap the coord input argument as a Python array.
  npy_intp dims_coord[] = {8};
  PyObject * pPo = PyArray_SimpleNewFromData(1, dims_coord, NPY_DOUBLE, const_cast<double*>(coord));

  // Execute the Python isStopCondition method. This creates a new
  // Python object *pR to hold tge result.
  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pIsStopCondition_, pPo, NULL);
  isRecursing=false;

  // Release the Python array holding the coord input argument.
  Py_XDECREF(pPo);

  // Check for errors.
  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::isStopCondition()");
  }

  // Interpret the result depending on its Python type.
  int res = 0;
  if (PyLong_Check(pR)) {
    GYOTO_DEBUG << "Python code returned a Long, return it unchanged.\n";
    res =  PyLong_AsLong(pR);
  } else {
    GYOTO_DEBUG << "Python code did not return a Long, return truth value.\n";
    res = PyObject_IsTrue(pR);
  }

  // Check for errors.
  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::isStopCondition()");
  }

  // Release all temporary Python objects as well as the GIL.
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  // Return the result.
  return res;
}

void Metric::Python::circularVelocity(double const pos[4], double vel[4],
				      double dir) const
{
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pCircularVelocity_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Metric::Python*>(this)->pCircularVelocity_);
  }

  GYOTO_DEBUG_EXPR(pCircularVelocity_);
  if (!pCircularVelocity_ || keplerian_) {
    Metric::Generic::circularVelocity(pos, vel, dir);
    return;
  }

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims[] = {4};

  PyObject * pP = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, const_cast<double*>(pos));
  PyObject * pV = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, vel);
  PyObject * pD = PyFloat_FromDouble(dir);

  isRecursing=true;
  PyObject * pR = PyObject_CallFunctionObjArgs(pCircularVelocity_, pP, pV, pD, NULL);
  isRecursing=false;

  Py_XDECREF(pR);
  Py_XDECREF(pD);
  Py_XDECREF(pV);
  Py_XDECREF(pP);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Metric::Python::circularVelocity");
  }

  PyGILState_Release(gstate);
}
