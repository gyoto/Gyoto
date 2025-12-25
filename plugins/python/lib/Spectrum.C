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
GYOTO_PROPERTY_SIZE_T(Gyoto::Spectrum::Python, Instance, instance,
      "ID of pre-constructed Class instance as per gyoto.core.gyotoid().")
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
void Spectrum::Python::klass(const std::string &c)
{Gyoto::Python::Base::klass(c);}

void Spectrum::Python::instance(size_t i)
{Gyoto::Python::Base::instance(reinterpret_cast<PyObject*>(i));}
size_t Spectrum::Python::instance() const
{return reinterpret_cast<size_t>(pInstance_);}

void Spectrum::Python::detachInstance() {
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  Py_CLEAR(pCall_);
  Py_CLEAR(pIntegrate_);
  Gyoto::Python::Base::detachInstance();
}


void Gyoto::Spectrum::Python::attachInstance (PyObject * instance) {
  Gyoto::Python::Base::attachInstance(instance);
  GYOTO_DEBUG << "Checking methodes for Python class " << class_ << endl;

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  pCall_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "__call__");
  pIntegrate_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "integrate");

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error while retrieving methods");
  }

  pCall_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "__call__");
  pIntegrate_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "integrate");

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error while retrieving methods");
  }

  if (!pCall_) {
    GYOTO_ERROR("Object does not implement required method \"__call__\"");
  }

  pCall_overloaded_ = Gyoto::Python::PyCallable_HasVarArg(pCall_);

  Gyoto::Python::PyInstance_SetThis(pInstance_,
				    Gyoto::Python::pGyotoSpectrum(),
				    this);
  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error while setting this");
  }

  if (parameters_.size()) parameters(parameters_);
  GYOTO_DEBUG << "Done checking Python class methods" << class_ << endl;
}

/* Gyoto::Spectrum::Generic API */

double Spectrum::Python::operator()(double nu) const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pCall_"
		<< std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Spectrum::Python*>(this)->pCall_);
  }

  GYOTO_DEBUG_EXPR(pCall_);
  if (!pCall_) GYOTO_ERROR("class does not implement required method '__call__'");
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyObject * pArgs = Py_BuildValue("(d)", nu);
  guardian.track(pArgs);

  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    GYOTO_ERROR("Failed building argument list");
  }

  isRecursing=true;
  PyObject * pValue = PyObject_CallObject(pCall_, pArgs);
  isRecursing=false;
  guardian.track(pValue);

  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    GYOTO_ERROR("Failed calling Python method __call__");
  }

  double res = PyFloat_AsDouble(pValue);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error interpreting result as double");
  }

  return res;
}

double Spectrum::Python::operator()(double nu, double opacity, double ds) const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pCall_ and pCall_overloaded_"
		<< std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Spectrum::Python*>(this)->pCall_);
    const_cast<Gyoto::Spectrum::Python*>(this)->pCall_overloaded_=false;
  }

  if (!pCall_) GYOTO_ERROR("class does not implement required method '__call__'");
  if (!pCall_overloaded_) return Generic::operator()(nu, opacity, ds);

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyObject * pArgs = Py_BuildValue("(ddd)", nu, opacity, ds);
  guardian.track(pArgs);

  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    GYOTO_ERROR("Failed building argument list");
  }

  isRecursing=true;
  PyObject * pValue = PyObject_CallObject(pCall_, pArgs);
  isRecursing=false;
  guardian.track(pValue);

  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    GYOTO_ERROR("Failed calling Python method __call__");
  }

  double res = PyFloat_AsDouble(pValue);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error interpreting result as double");
  }

  return res;

}

double Spectrum::Python::integrate(double nu1, double nu2) {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pIntegrate_"
		<< std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Spectrum::Python*>(this)->pIntegrate_);
  }

  if (!pIntegrate_) return Generic::integrate(nu1, nu2);

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;

  PyObject * pArgs = Py_BuildValue("dd", nu1, nu2);
  guardian.track(pArgs);

  if (PyErr_Occurred() || !pArgs) {
    PyErr_Print();
    GYOTO_ERROR("Failed building argument list");
  }

  isRecursing=true;
  PyObject * pValue = PyObject_CallObject(pIntegrate_, pArgs);
  isRecursing=true;
  guardian.track(pValue);

  if (PyErr_Occurred() || !pValue) {
    PyErr_Print();
    GYOTO_ERROR("Failed calling Python method integrate");
  }

  double res = PyFloat_AsDouble(pValue);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error interpreting result as double");
  }

  return res;
}
