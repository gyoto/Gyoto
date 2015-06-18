#include "GyotoPython.h"

#include <Python.h>

using namespace Gyoto;
using namespace Gyoto::Python;
using namespace std;

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
