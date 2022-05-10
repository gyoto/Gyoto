#include "GyotoPython.h"
#include "GyotoProperty.h"
#include "GyotoError.h"

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL GyotoPython_ARRAY_API
#include <numpy/arrayobject.h>

using namespace Gyoto;
using namespace std;

GYOTO_PROPERTY_START(Gyoto::Astrobj::Python::Standard,
      "Python-based Astrobj class")
GYOTO_PROPERTY_STRING(Gyoto::Astrobj::Python::Standard, Module, module,
      "Python module containing the Astrobj implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Astrobj::Python::Standard, InlineModule, inlineModule,
      "Inline code of Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Astrobj::Python::Standard, Class, klass,
      "Python class (in Module) implementing the Astrobj.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Gyoto::Astrobj::Python::Standard,
			     Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_DOUBLE(Gyoto::Astrobj::Python::Standard,
		       CriticalValue, criticalValue,
      "The object is defined by __call__ < this value")
GYOTO_PROPERTY_END(Astrobj::Python::Standard, Astrobj::Standard::properties)

GYOTO_PROPERTY_THREAD_UNSAFE(Astrobj::Python::Standard)

// Birth and death
Gyoto::Astrobj::Python::Standard::Standard()
: Gyoto::Python::Object<Astrobj::Standard>(),
  pEmission_(NULL), pIntegrateEmission_(NULL), pTransmission_(NULL),
  pCall_(NULL), pGetVelocity_(NULL), pGiveDelta_(NULL),
  pEmission_overloaded_(false), pIntegrateEmission_overloaded_(false)
{
  kind("Python::Standard");
}

Gyoto::Astrobj::Python::Standard::Standard(const Standard& o)
  : Gyoto::Python::Object<Astrobj::Standard>(o),
  pEmission_(o.pEmission_), pIntegrateEmission_(o.pIntegrateEmission_),
  pTransmission_(o.pTransmission_), pCall_(o.pCall_),
  pGetVelocity_(o.pGetVelocity_), pGiveDelta_(o.pGiveDelta_),
  pEmission_overloaded_(o.pEmission_overloaded_),
  pIntegrateEmission_overloaded_(o.pIntegrateEmission_overloaded_)
{
  Py_XINCREF(pEmission_);
  Py_XINCREF(pIntegrateEmission_);
  Py_XINCREF(pTransmission_);
  Py_XINCREF(pCall_);
  Py_XINCREF(pGetVelocity_);
  Py_XINCREF(pGiveDelta_);
}

Gyoto::Astrobj::Python::Standard::~Standard() {
  Py_XDECREF(pEmission_);
  Py_XDECREF(pIntegrateEmission_);
  Py_XDECREF(pTransmission_);
  Py_XDECREF(pCall_);
  Py_XDECREF(pGetVelocity_);
  Py_XDECREF(pGiveDelta_);
}

Astrobj::Python::Standard* Gyoto::Astrobj::Python::Standard::clone() const
{return new Standard(*this);}

double Astrobj::Python::Standard::criticalValue() const
{return critical_value_;}
void Astrobj::Python::Standard::criticalValue(double v)
{critical_value_=v;}
  
std::vector<double> Astrobj::Python::Standard::parameters() const
{return Gyoto::Python::Base::parameters();}
void Astrobj::Python::Standard::parameters(const std::vector<double>& m)
{Gyoto::Python::Base::parameters(m);}
std::string Astrobj::Python::Standard::module() const
{return Gyoto::Python::Base::module();}
std::string Astrobj::Python::Standard::inlineModule() const
{return Gyoto::Python::Base::inlineModule();}
void Astrobj::Python::Standard::inlineModule(const std::string& m)
{Gyoto::Python::Base::inlineModule(m);}
void Astrobj::Python::Standard::module(const std::string& m)
{Gyoto::Python::Base::module(m);}
std::string Astrobj::Python::Standard::klass() const
{return Gyoto::Python::Base::klass();}

void Gyoto::Astrobj::Python::Standard::klass(const std::string &f) {

  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pEmission_);
  Py_XDECREF(pIntegrateEmission_);
  Py_XDECREF(pTransmission_);
  Py_XDECREF(pCall_);
  Py_XDECREF(pGetVelocity_);
  Py_XDECREF(pGiveDelta_);
  PyGILState_Release(gstate);

  pEmission_overloaded_ = false;
  pIntegrateEmission_overloaded_ = false;
  
  Gyoto::Python::Base::klass(f);
  if (!pModule_) return;

  gstate = PyGILState_Ensure();
  GYOTO_DEBUG << "Checking Python class methods" << f << endl;

  pEmission_          =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "emission");
  pIntegrateEmission_ =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "integrateEmission");
  pTransmission_      =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "transmission");
  pCall_              =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "__call__");
  pGetVelocity_       =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "getVelocity");
  pGiveDelta_         =
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "giveDelta");
  
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error while retrieving methods");
  }

  if (!pCall_) {
    PyGILState_Release(gstate);
    GYOTO_ERROR("Object does not implement required method \"__call__\"");
  }

  if (!pGetVelocity_) {
    PyGILState_Release(gstate);
    GYOTO_ERROR("Object does not implement required method \"getVelocity\"");
  }

  pEmission_overloaded_ = pEmission_ &&
    Gyoto::Python::PyCallable_HasVarArg(pEmission_);
  
  pIntegrateEmission_overloaded_ = pIntegrateEmission_ &&
    Gyoto::Python::PyCallable_HasVarArg(pIntegrateEmission_);

  Gyoto::Python::PyInstance_SetThis(pInstance_,
				    Gyoto::Python::pGyotoStandardAstrobj(),
				    this);

  PyGILState_Release(gstate);
  if (parameters_.size()) parameters(parameters_);
  GYOTO_DEBUG << "Done checking Python class methods" << f << endl;
}

double Gyoto::Astrobj::Python::Standard::operator()(double const coord[4]) {
  if (!pCall_) GYOTO_ERROR("__call__ not loaded yet");
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims[] = {4};
  
  PyObject * pCoord = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE,
						const_cast<double*>(coord));
  PyObject * pR = PyObject_CallFunctionObjArgs(pCall_, pCoord, NULL);
  double res = PyFloat_AsDouble(pR);
  
  Py_XDECREF(pR);
  Py_XDECREF(pCoord);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::operator()()");
  }
   
  PyGILState_Release(gstate);
  return res;
}

void Gyoto::Astrobj::Python::Standard::getVelocity
(double const coord[4], double vel[4]) {
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims[] = {4};
  
  PyObject * pCoord = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE,
						const_cast<double*>(coord));
  PyObject * pVel = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, vel);
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pGetVelocity_, pCoord, pVel, NULL);
  
  Py_XDECREF(pR);
  Py_XDECREF(pCoord);
  Py_XDECREF(pVel);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::getVelocity()");
  }
   
  PyGILState_Release(gstate);
}

double Gyoto::Astrobj::Python::Standard::giveDelta(double coord[8]) {
  if (!pGiveDelta_) return Astrobj::Standard::giveDelta(coord);
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims[] = {8};
  
  PyObject * pCoord = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, coord);
  PyObject * pR = PyObject_CallFunctionObjArgs(pGiveDelta_, pCoord, NULL);
  double res = PyFloat_AsDouble(pR);
  
  Py_XDECREF(pR);
  Py_XDECREF(pCoord);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::giveDelta()");
  }
   
  PyGILState_Release(gstate);
  return res;
}

double Gyoto::Astrobj::Python::Standard::emission
(double nu_em, double dsem, state_t const &coord_ph, double const coord_obj[8])
  const {
  if (!pEmission_)
    return Astrobj::Standard::emission(nu_em, dsem, coord_ph, coord_obj);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims_co[] = {8};
  npy_intp dims_ph[] = {npy_intp(coord_ph.size())};
  
  PyObject * pNu = PyFloat_FromDouble(nu_em);
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, dims_ph, NPY_DOUBLE, const_cast<double*>(&coord_ph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, dims_co, NPY_DOUBLE, const_cast<double*>(coord_obj));
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pEmission_, pNu, pDs, pCp, pCo, NULL);

  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pNu);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::emission()");
  }
 
  double res = PyFloat_AsDouble(pR);  
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

void Gyoto::Astrobj::Python::Standard::emission
(double Inu[], double const nu_em[], size_t nbnu, double dsem, state_t const &coord_ph,
 double const coord_obj[8]) const {
  if (!pEmission_ || !pEmission_overloaded_) {
    Astrobj::Standard::emission(Inu, nu_em, nbnu, dsem,
				coord_ph, coord_obj);
    return;
  }

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp I_dims[] = {static_cast<npy_intp>(nbnu)};
  npy_intp dims_co[] = {8};
  npy_intp dims_cp[] = {npy_intp(coord_ph.size())};
  
  PyObject * pIn = PyArray_SimpleNewFromData(1, I_dims, NPY_DOUBLE, Inu);
  PyObject * pNu = PyArray_SimpleNewFromData(1, I_dims, NPY_DOUBLE, const_cast<double*>(nu_em));
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, dims_cp, NPY_DOUBLE, const_cast<double*>(&coord_ph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, dims_co, NPY_DOUBLE, const_cast<double*>(coord_obj));
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pEmission_, pIn, pNu, pDs, pCp, pCo, NULL);
  
  Py_XDECREF(pR);
  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pNu);
  Py_XDECREF(pIn);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::emission()");
  }
   
  PyGILState_Release(gstate);
}

double Gyoto::Astrobj::Python::Standard::integrateEmission
(double nu1, double nu2, double dsem, state_t const &c_ph, double const c_obj[8])
  const {
  if (!pIntegrateEmission_)
    return Astrobj::Standard::integrateEmission(nu1, nu2, dsem, c_ph,c_obj);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims_co[] = {8};
  npy_intp dims_cp[] = {npy_intp(c_ph.size())};
  
  PyObject * pN1 = PyFloat_FromDouble(nu1);
  PyObject * pN2 = PyFloat_FromDouble(nu2);
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, dims_cp, NPY_DOUBLE, const_cast<double*>(&c_ph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, dims_co, NPY_DOUBLE, const_cast<double*>(c_obj));
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pIntegrateEmission_,
				 pN1, pN2, pDs, pCp, pCo, NULL);

  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pN2);
  Py_XDECREF(pN1);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::integrateEmission()");
  }
   
  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

void Gyoto::Astrobj::Python::Standard::integrateEmission
(double * I, double const * boundaries, size_t const * chaninds,
 size_t nbnu, double dsem, state_t const &cph, const double *co) const{
  if (!pIntegrateEmission_ || !pIntegrateEmission_overloaded_) {
    Gyoto::Astrobj::Standard::integrateEmission(I, boundaries, chaninds,
						nbnu, dsem, cph, co);
    return;
  }

  PyGILState_STATE gstate = PyGILState_Ensure();

  size_t nbo=0;
  for (size_t i=0; i<2*nbnu; ++i)
    if (nbo < chaninds[i]) nbo=chaninds[i];
  npy_intp nNu = nbnu, nBo = nbo, nCh = 2*nbnu, nCo = 8, nCp = npy_intp(cph.size());

  PyObject * pI  = PyArray_SimpleNewFromData(1, &nNu, NPY_DOUBLE, I);
  PyObject * pBo = PyArray_SimpleNewFromData(1, &nBo, NPY_DOUBLE,
					     const_cast<double*>(boundaries));
  PyObject * pCh = PyArray_SimpleNewFromData(1, &nCh, NPY_UINTP,
					     const_cast<size_t *>(chaninds));
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, &nCp, NPY_DOUBLE, const_cast<double*>(&cph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, &nCo, NPY_DOUBLE, const_cast<double*>(co));
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pIntegrateEmission_,
				 pI, pBo, pCh, pDs, pCp, pCo, NULL);
  
  Py_XDECREF(pR);
  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pCh);
  Py_XDECREF(pBo);
  Py_XDECREF(pI);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::integrateEmission()");
  }
   
  PyGILState_Release(gstate);
}

double Gyoto::Astrobj::Python::Standard::transmission
(double nuem, double dsem, state_t const &cph, double const *co) const {
  if (!pTransmission_)
    return Astrobj::Standard::transmission(nuem, dsem, cph, co);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp pdims[] = {npy_intp(cph.size())};
  npy_intp odims[] = {8};
  
  PyObject * pNu = PyFloat_FromDouble(nuem);
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, pdims, NPY_DOUBLE, const_cast<double*>(&cph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, odims, NPY_DOUBLE, const_cast<double*>(co));
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pTransmission_, pNu, pDs, pCp, pCo, NULL);

  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pNu);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in Standard::transmission()");
  }
   
  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

