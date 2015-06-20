#include "GyotoPython.h"

#include <Python.h>

using namespace Gyoto;
using namespace Gyoto::Python;
using namespace std;

PyObject * Gyoto::Python::PyInstance_GetMethod
(PyObject* pInstance, const char *name) {
  PyObject * pName = PyUnicode_FromString(name);
  if (!pName) return NULL;
  if (!PyObject_HasAttr(pInstance, pName)) {
    Py_XDECREF(pName);
    return NULL;
  }
  PyObject * pMethod = PyObject_GetAttr(pInstance, pName);
  Py_DECREF(pName);
  if (!pMethod) return NULL;
  if (!PyCallable_Check(pMethod)) {
    Py_DECREF(pMethod);
    return NULL;
  }
  return pMethod;
}

bool Gyoto::Python::PyCallable_HasVarArg(PyObject * pMethod) {
  static PyObject * pGetArgSpec = NULL;

  if (!pGetArgSpec) {
    PyObject * pName = PyUnicode_FromString("inspect");
    PyObject * pModule = PyImport_Import(pName);
    Py_XDECREF(pName); pName=NULL;
    pGetArgSpec = PyObject_GetAttrString(pModule, "getargspec");
  }

  PyObject * pArgSpec =
    PyObject_CallFunctionObjArgs(pGetArgSpec, pMethod, NULL);
  bool answer = (PyTuple_GetItem(pArgSpec, 1) != Py_None);
  Py_XDECREF(pArgSpec);

  return answer;
}

PyObject * Gyoto::Python::PyImport_Gyoto() {
  static PyObject * pModule = NULL;
  static bool need_load = true;

  if (need_load) {
    need_load=false;
    pModule = PyImport_ImportModule("gyoto");
    if (PyErr_Occurred()) {
      GYOTO_WARNING << "";
      PyErr_Print();
    }
  }

  return pModule;
}

PyObject * Gyoto::Python::pGyotoSpectrum() {
  PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "Spectrum");
  }
  return res;
}

PyObject * Gyoto::Python::pGyotoThinDisk() {
  PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "ThinDisk");
  }
  return res;
}

PyObject * Gyoto::Python::pGyotoMetric() {
  PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "Metric");
  }
  return res;
}

PyObject * Gyoto::Python::pGyotoStandardAstrobj() {
  PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "StandardAstrobj");
  }
  return res;
}

void Gyoto::Python::PyInstance_SetThis(PyObject * pInstance,
				       PyObject * pNew,
				       void * ptr) {
  PyObject * pThis = NULL;
  if (pNew) pThis = PyObject_CallFunction(pNew, "l", (long)ptr);
  else {
    pThis = Py_None;
    Py_INCREF(pThis);
  }
  PyObject_SetAttrString(pInstance, "this", pThis);
  Py_XDECREF(pThis);
}


// Birth and death
Base::Base()
: module_(""), class_(""),   parameters_(),
  pModule_(NULL), pInstance_(NULL)
{}

Base::Base(const Base& o)
: module_(o.module_), class_(o.class_), parameters_(o.parameters_),
  pModule_(o.pModule_), pInstance_(o.pInstance_)
{
  Py_XINCREF(pModule_);
  Py_XINCREF(pInstance_);
}

Base::~Base() {
  Py_XDECREF(pInstance_);
  Py_XDECREF(pModule_);
}

/* Accessors */

std::string Base::module() const { return module_; }
void Base::module(const std::string &m) {
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
  if (PyErr_Occurred() || !pModule_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    throwError("Failed loading Python module");
  }
  PyGILState_Release(gstate);
  if (class_ != "") klass(class_);
  GYOTO_DEBUG << "Done loading Python module " << m << endl;
}

std::string Base::klass() const { return class_; }
void Base::klass(const std::string &f) {
  class_=f;
  if (!pModule_) return;

  GYOTO_DEBUG << "Instanciating Python class " << f << endl;
  
  PyGILState_STATE gstate = PyGILState_Ensure();
  
  Py_XDECREF(pInstance_); pInstance_=NULL;

  PyObject * pClass = PyObject_GetAttrString(pModule_, class_.c_str());
  if (PyErr_Occurred() || !pClass) {
    PyErr_Print();
    Py_XDECREF(pClass);
    PyGILState_Release(gstate);
    throwError("Could not find class in module");
  }
  if (!PyCallable_Check(pClass)) {
    Py_DECREF(pClass);
    PyGILState_Release(gstate);
    throwError("Class is not callable");
  }

  pInstance_ = PyObject_CallObject(pClass, NULL);
  Py_DECREF(pClass); pClass=NULL;
  if (PyErr_Occurred() || !pInstance_) {
    PyErr_Print();
    Py_XDECREF(pInstance_); pInstance_=NULL;
    PyGILState_Release(gstate);
    throwError("Failed instanciating Python class");
  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "Done instanciating Python class " << f << endl;
}

std::vector<double> Base::parameters() const {return parameters_;}
void Base::parameters(const std::vector<double> &p){
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
      throwError("Failed calling __setitem__");
    }

  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "done.\n";
}
