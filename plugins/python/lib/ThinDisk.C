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

GYOTO_PROPERTY_START(Gyoto::Astrobj::Python::ThinDisk,
      "Python-based Astrobj class")
GYOTO_PROPERTY_STRING(Gyoto::Astrobj::Python::ThinDisk, Module, module,
      "Python module containing the Astrobj implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Astrobj::Python::ThinDisk, InlineModule, inlineModule,
      "Inline code of Python module containing the Spectrum implementation.")
GYOTO_PROPERTY_STRING(Gyoto::Astrobj::Python::ThinDisk, Class, klass,
      "Python class (in Module) implementing the Astrobj.")
GYOTO_PROPERTY_SIZE_T(Gyoto::Astrobj::Python::ThinDisk, Instance, instance,
      "Memory address of Python instance implementing the Astrobj as returned by id().")
GYOTO_PROPERTY_VECTOR_DOUBLE(Gyoto::Astrobj::Python::ThinDisk,
			     Parameters, parameters,
      "Parameters for the class instance.")
GYOTO_PROPERTY_END(Astrobj::Python::ThinDisk, Astrobj::ThinDisk::properties)

GYOTO_PROPERTY_THREAD_UNSAFE(Astrobj::Python::ThinDisk)

// Birth and death
Gyoto::Astrobj::Python::ThinDisk::ThinDisk()
: Gyoto::Python::Object<Astrobj::ThinDisk>(),
  pEmission_(NULL), pIntegrateEmission_(NULL), pTransmission_(NULL),
  pCall_(NULL), pGetVelocity_(NULL),
  pEmission_overloaded_(false), pIntegrateEmission_overloaded_(false)
{
  kind("Python::ThinDisk");
}

Gyoto::Astrobj::Python::ThinDisk::ThinDisk(const ThinDisk& o)
  :
  Gyoto::Python::Object<Astrobj::ThinDisk>(o),
  pEmission_(o.pEmission_), pIntegrateEmission_(o.pIntegrateEmission_),
  pTransmission_(o.pTransmission_), pCall_(o.pCall_),
  pGetVelocity_(o.pGetVelocity_),
  pEmission_overloaded_(o.pEmission_overloaded_),
  pIntegrateEmission_overloaded_(o.pIntegrateEmission_overloaded_)
{
  Py_XINCREF(pEmission_);
  Py_XINCREF(pIntegrateEmission_);
  Py_XINCREF(pTransmission_);
  Py_XINCREF(pCall_);
  Py_XINCREF(pGetVelocity_);
}

Gyoto::Astrobj::Python::ThinDisk::~ThinDisk() {
  Py_CLEAR(pEmission_);
  Py_CLEAR(pIntegrateEmission_);
  Py_CLEAR(pTransmission_);
  Py_CLEAR(pCall_);
  Py_CLEAR(pGetVelocity_);
}

Astrobj::Python::ThinDisk* Gyoto::Astrobj::Python::ThinDisk::clone() const
{return new ThinDisk(*this);}

std::vector<double> Astrobj::Python::ThinDisk::parameters() const
{return Gyoto::Python::Base::parameters();}
void Astrobj::Python::ThinDisk::parameters(const std::vector<double>& m)
{Gyoto::Python::Base::parameters(m);}
std::string Astrobj::Python::ThinDisk::module() const
{return Gyoto::Python::Base::module();}
void Astrobj::Python::ThinDisk::module(const std::string& m)
{Gyoto::Python::Base::module(m);}
std::string Astrobj::Python::ThinDisk::inlineModule() const
{return Gyoto::Python::Base::inlineModule();}
void Astrobj::Python::ThinDisk::inlineModule(const std::string& m)
{Gyoto::Python::Base::inlineModule(m);}
std::string Astrobj::Python::ThinDisk::klass() const
{return Gyoto::Python::Base::klass();}
void Astrobj::Python::ThinDisk::klass(const std::string& c)
{Gyoto::Python::Base::klass(c);}

void Astrobj::Python::ThinDisk::instance(size_t i)
{Gyoto::Python::Base::instance(reinterpret_cast<PyObject*>(i));}
size_t Astrobj::Python::ThinDisk::instance() const
{return reinterpret_cast<size_t>(pInstance_);}

void Gyoto::Astrobj::Python::ThinDisk::detachInstance() {
  pEmission_overloaded_ = false;
  pIntegrateEmission_overloaded_ = false;
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  Py_CLEAR(pEmission_);
  Py_CLEAR(pIntegrateEmission_);
  Py_CLEAR(pTransmission_);
  Py_CLEAR(pCall_);
  Py_CLEAR(pGetVelocity_);
  Gyoto::Python::Base::detachInstance();
}

void Gyoto::Astrobj::Python::ThinDisk::attachInstance (PyObject * instance) {
  Gyoto::Python::Base::attachInstance(instance);
  GYOTO_DEBUG << "Checking methodes for Python class " << class_ << endl;

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

  GYOTO_DEBUG_EXPR(pEmission_);
  GYOTO_DEBUG_EXPR(pIntegrateEmission_);
  GYOTO_DEBUG_EXPR(pTransmission_);
  GYOTO_DEBUG_EXPR(pCall_);
  GYOTO_DEBUG_EXPR(pGetVelocity_);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error while retrieving methods");
  }

  pEmission_overloaded_ = pEmission_ &&
    Gyoto::Python::PyCallable_HasVarArg(pEmission_);

  pIntegrateEmission_overloaded_ = pIntegrateEmission_ &&
    Gyoto::Python::PyCallable_HasVarArg(pIntegrateEmission_);

  GYOTO_DEBUG << "setting 'this' in the instance" << endl;
  Gyoto::Python::PyInstance_SetThis(pInstance_,
				    Gyoto::Python::pGyotoThinDisk(),
				    this);

  if (parameters_.size()) parameters(parameters_);
  GYOTO_DEBUG << "Done checking methods for Python class " << class_ << endl;
}

double Gyoto::Astrobj::Python::ThinDisk::operator()(double const coord[4]) {
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
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pCall_);
  }

  GYOTO_DEBUG_EXPR(pCall_);
  if (!pCall_) return Gyoto::Astrobj::ThinDisk::operator()(coord);
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims[] = {4};

  PyObject * pCoord = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE,
						const_cast<double*>(coord));
  isRecursing=true;
  PyObject * pR = PyObject_CallFunctionObjArgs(pCall_, pCoord, NULL);
  isRecursing=false;

  Py_XDECREF(pCoord);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in ThinDisk::operator()()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

void Gyoto::Astrobj::Python::ThinDisk::getVelocity
(double const coord[4], double vel[4]) {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pGetVelocity_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pGetVelocity_);
  }

  GYOTO_DEBUG_EXPR(pGetVelocity_);
  if (!pGetVelocity_) return Gyoto::Astrobj::ThinDisk::getVelocity(coord, vel);
  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims[] = {4};

  PyObject * pCoord = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE,
						const_cast<double*>(coord));
  PyObject * pVel = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, vel);

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pGetVelocity_, pCoord, pVel, NULL);
  isRecursing=false;

  Py_XDECREF(pR);
  Py_XDECREF(pCoord);
  Py_XDECREF(pVel);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in ThinDisk::getVelocity()");
  }

  PyGILState_Release(gstate);
}

double Gyoto::Astrobj::Python::ThinDisk::emission
(double nu_em, double dsem, state_t const &coord_ph, double const coord_obj[8])
  const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pEmission_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pEmission_);
  }

  GYOTO_DEBUG_EXPR(pEmission_);
  if (!pEmission_)
    return Astrobj::ThinDisk::emission(nu_em, dsem, coord_ph, coord_obj);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims_co[] = {8};
  npy_intp dims_cp[] = {npy_intp(coord_ph.size())};

  PyObject * pNu = PyFloat_FromDouble(nu_em);
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, dims_cp, NPY_DOUBLE, const_cast<double*>(&coord_ph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, dims_co, NPY_DOUBLE, const_cast<double*>(coord_obj));

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pEmission_, pNu, pDs, pCp, pCo, NULL);
  isRecursing=false;

  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pNu);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in ThinDisk::emission()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

void Gyoto::Astrobj::Python::ThinDisk::emission
(double Inu[], double const nu_em[], size_t nbnu, double dsem, state_t const &coord_ph,
 double const coord_obj[8]) const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pEmission_ "
		     << "and pEmission_overloaded_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pEmission_);
    const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pEmission_overloaded_ = false;
  }

  GYOTO_DEBUG_EXPR(pEmission_);
  GYOTO_DEBUG_EXPR(pEmission_overloaded_);
  if (!pEmission_ || !pEmission_overloaded_) {
    Astrobj::ThinDisk::emission(Inu, nu_em, nbnu, dsem,
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

  isRecursing = true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pEmission_, pIn, pNu, pDs, pCp, pCo, NULL);
  isRecursing = false;

  Py_XDECREF(pR);
  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pNu);
  Py_XDECREF(pIn);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in ThinDisk::emission()");
  }

  PyGILState_Release(gstate);
}

double Gyoto::Astrobj::Python::ThinDisk::integrateEmission
(double nu1, double nu2, double dsem, state_t const &c_ph, double const c_obj[8])
  const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pIntegrateEmission_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pIntegrateEmission_);
  }

  GYOTO_DEBUG_EXPR(pIntegrateEmission_);
  if (!pIntegrateEmission_)
    return Astrobj::ThinDisk::integrateEmission(nu1, nu2, dsem, c_ph,c_obj);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp dims_co[] = {8};
  npy_intp dims_cp[] = {npy_intp(c_ph.size())};

  PyObject * pN1 = PyFloat_FromDouble(nu1);
  PyObject * pN2 = PyFloat_FromDouble(nu2);
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, dims_cp, NPY_DOUBLE, const_cast<double*>(&c_ph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, dims_co, NPY_DOUBLE, const_cast<double*>(c_obj));

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pIntegrateEmission_,
				 pN1, pN2, pDs, pCp, pCo, NULL);
  isRecursing=false;

  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pN2);
  Py_XDECREF(pN1);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in ThinDisk::integrateEmission()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}

void Gyoto::Astrobj::Python::ThinDisk::integrateEmission
(double * I, double const * boundaries, size_t const * chaninds,
 size_t nbnu, double dsem, state_t const &cph, double const *co) const{
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pIntegrateEmission_ "
		     << "and pIntegrateEmission_overloaded_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pIntegrateEmission_);
    const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pIntegrateEmission_overloaded_ = false;
  }

  GYOTO_DEBUG_EXPR(pIntegrateEmission_);
  GYOTO_DEBUG_EXPR(pIntegrateEmission_overloaded_);
  if (!pIntegrateEmission_ || !pIntegrateEmission_overloaded_) {
    Gyoto::Astrobj::ThinDisk::integrateEmission(I, boundaries, chaninds,
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
  PyObject * pCp = PyArray_SimpleNewFromData(1, &nCp, NPY_DOUBLE,
					     const_cast<double*>(&cph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, &nCo, NPY_DOUBLE,
					     const_cast<double*>(co));

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pIntegrateEmission_,
				 pI, pBo, pCh, pDs, pCp, pCo, NULL);
  isRecursing=false;

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
    GYOTO_ERROR("Error occurred in ThinDisk::integrateEmission()");
  }

  PyGILState_Release(gstate);
}

double Gyoto::Astrobj::Python::ThinDisk::transmission
(double nuem, double dsem, state_t const &cph, double const * co) const {
  /* Recursion check:

     When a class derives from StandardBase, it inherits several
     methods from PythonStandard. In that case, if the class doesn't
     reimplement it, these methods will recurse indefinitely. Here we
     break this loop.
   */
  thread_local bool isRecursing = false;
  if (isRecursing) {
    GYOTO_DEBUG << "recursion detected, resetting pTransmission_"
		     << std::endl;
    // ugly but right: bypass constness
    Py_CLEAR(const_cast<Gyoto::Astrobj::Python::ThinDisk*>(this)->pTransmission_);
  }

  GYOTO_DEBUG_EXPR(pTransmission_);
  if (!pTransmission_)
    return Astrobj::ThinDisk::transmission(nuem, dsem, cph, co);

  PyGILState_STATE gstate = PyGILState_Ensure();

  npy_intp pdims[] = {npy_intp(cph.size())};
  npy_intp odims[] = {8};

  PyObject * pNu = PyFloat_FromDouble(nuem);
  PyObject * pDs = PyFloat_FromDouble(dsem);
  PyObject * pCp = PyArray_SimpleNewFromData(1, pdims, NPY_DOUBLE, const_cast<double*>(&cph[0]));
  PyObject * pCo = PyArray_SimpleNewFromData(1, odims, NPY_DOUBLE, const_cast<double*>(co));

  isRecursing=true;
  PyObject * pR =
    PyObject_CallFunctionObjArgs(pTransmission_, pNu, pDs, pCp, pCo, NULL);
  isRecursing=false;

  Py_XDECREF(pCo);
  Py_XDECREF(pCp);
  Py_XDECREF(pDs);
  Py_XDECREF(pNu);

  if (PyErr_Occurred()) {
    Py_XDECREF(pR);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred in ThinDisk::transmission()");
  }

  double res = PyFloat_AsDouble(pR);
  Py_XDECREF(pR);
  PyGILState_Release(gstate);

  return res;
}
