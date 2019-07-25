/*
    Copyright © 2015-2017, 2019 Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.

 */

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
    pModule = PyImport_ImportModule("gyoto.core");
    if (PyErr_Occurred()) {
      GYOTO_WARNING << "";
      PyErr_Print();
    }
  }

  return pModule;
}

PyObject * Gyoto::Python::pGyotoSpectrum() {
  static PyObject * res = NULL;
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
  static PyObject * res = NULL;
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
  static PyObject * res = NULL;
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
  static PyObject * res = NULL;
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

PyObject * Gyoto::Python::PyModule_NewFromPythonCode(const char * source_code) {
  PyObject * pDedent = NULL;
  if (!pDedent) {
    GYOTO_DEBUG << "importing textwrap.dedent\n";
    PyObject * mod = PyImport_ImportModule("textwrap");
    if (PyErr_Occurred() || !mod) {
      Py_XDECREF(mod);
      return NULL;
    }
    pDedent = PyObject_GetAttrString(mod, "dedent");
    Py_XDECREF(mod);
    if (PyErr_Occurred() || !pDedent) {
      return NULL;
    }
    GYOTO_DEBUG << "done importing textwrap.dedent\n";
  }

  GYOTO_DEBUG << "dedenting source code... \n";
  PyObject * pCode = PyObject_CallFunction(pDedent, "s", source_code);
  if (PyErr_Occurred() || !pCode) {
    GYOTO_DEBUG << "failed dedenting source code!\n";
    Py_XDECREF(pCode);
    return NULL;
  }

  const char * new_source_code = NULL;

  if (PyUnicode_Check(pCode)) {
    PyObject * tmp = PyUnicode_AsUTF8String(pCode);
    Py_DECREF(pCode);
    pCode = tmp;
  }

  if (!PyBytes_Check(pCode)) {
    GYOTO_DEBUG << "not a PyBytes string\n";
    Py_DECREF(pCode);
    return NULL;
  }

  new_source_code = PyBytes_AsString(pCode);

  GYOTO_DEBUG << "compiling inline code...\n";
  PyObject * object_code=Py_CompileString(new_source_code, "<inline>", Py_file_input);
  Py_DECREF(pCode);
  if (PyErr_Occurred() || !object_code) {
  GYOTO_DEBUG << "failed compiling inline code!\n";
    Py_XDECREF(object_code);
    return NULL;
  }

  GYOTO_DEBUG << "importing object code as module...\n";
  PyObject * mod = PyImport_ExecCodeModule("gyoto_inline", object_code);
  Py_XDECREF(object_code);
  if (PyErr_Occurred() || !mod) {
    GYOTO_DEBUG << "failed importing object code as module!\n";
    Py_XDECREF(mod);
    PyErr_Print();
    return NULL;
  }

  return mod;
}


// Birth and death
Base::Base()
  : module_(""), inline_module_(""), class_(""),   parameters_(),
  pModule_(NULL), pInstance_(NULL)
{}

Base::Base(const Base& o)
: module_(o.module_), inline_module_(o.inline_module_),
  class_(o.class_), parameters_(o.parameters_),
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
  if (m=="") return;
  inline_module_="";

  gstate = PyGILState_Ensure();
  PyObject *pName=PyUnicode_FromString(m.c_str());
  if (!pName) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed translating string to Python");
  }
  Py_XDECREF(pModule_);
  pModule_ = PyImport_Import(pName);
  Py_DECREF(pName);
  if (PyErr_Occurred() || !pModule_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed loading Python module");
  }
  PyGILState_Release(gstate);
  if (class_ != "") klass(class_);
  GYOTO_DEBUG << "Done loading Python module " << m << endl;
}

std::string Base::inlineModule() const { return inline_module_; }
void Base::inlineModule(const std::string &m) {
  inline_module_=m;
  if (m=="") return;
  module_="";

  GYOTO_DEBUG << "Loading inline Python module :" << m << endl;
  PyGILState_STATE gstate = PyGILState_Ensure();
  Py_XDECREF(pModule_);
  pModule_ = Gyoto::Python::PyModule_NewFromPythonCode(m.c_str());
  if (PyErr_Occurred() || !pModule_) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed loading inline Python module");
  }
  PyGILState_Release(gstate);
  if (class_ != "") klass(class_);
  GYOTO_DEBUG << "Done loading Python module " << m << endl;
}

std::string Base::klass() const { return class_; }
void Base::klass(const std::string &f) {
  class_=f;
  if (!pModule_) return;

  GYOTO_DEBUG << "Instantiating Python class " << f << endl;
  
  PyGILState_STATE gstate = PyGILState_Ensure();
  
  Py_XDECREF(pInstance_); pInstance_=NULL;

  if (class_ == ""){
    GYOTO_DEBUG << "class_ is empty: check whether there is a single class in module...\n";

    PyObject * dict = PyModule_GetDict(pModule_);
    PyObject *key, *value, *tmp;
    Py_ssize_t pos = 0, nclass=0;

    while (PyDict_Next(dict, &pos, &key, &value)) {
      if (
#if PY_VERSION_HEX < 0x03000000
	  PyClass_Check(value) ||
#endif
	  PyObject_TypeCheck(value, &PyType_Type)) {
	++nclass;
	if (PyUnicode_Check(key)) {
	  tmp = PyUnicode_AsUTF8String(key);
	} else {
	  tmp = key;
	  Py_INCREF(tmp);
	}

	if (!PyBytes_Check(tmp)) {
	  Py_DECREF(tmp);
	  PyGILState_Release(gstate);
	  GYOTO_ERROR("not a PyBytes string");
	}

	class_= PyBytes_AsString(tmp);

	Py_DECREF(tmp);
      }
    }
    if (nclass>1) {
      GYOTO_DEBUG << "several classes in module" << endl;
      class_ = "";
    } else if (nclass == 1) GYOTO_DEBUG << "single class in module: " << class_ << endl;
    else if (nclass == 0) {
      PyGILState_Release(gstate);
      GYOTO_ERROR("no class in Python module\n");
    }

  }

  PyObject * pClass = PyObject_GetAttrString(pModule_, class_.c_str());
  if (PyErr_Occurred() || !pClass) {
    PyErr_Print();
    Py_XDECREF(pClass);
    PyGILState_Release(gstate);
    GYOTO_ERROR("Could not find class in module");
  }
  if (!PyCallable_Check(pClass)) {
    Py_DECREF(pClass);
    PyGILState_Release(gstate);
    GYOTO_ERROR("Class is not callable");
  }

  pInstance_ = PyObject_CallObject(pClass, NULL);
  Py_DECREF(pClass); pClass=NULL;
  if (PyErr_Occurred() || !pInstance_) {
    PyErr_Print();
    Py_XDECREF(pInstance_); pInstance_=NULL;
    PyGILState_Release(gstate);
    GYOTO_ERROR("Failed instantiating Python class");
  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "Done instantiating Python class " << f << endl;
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
      GYOTO_ERROR("Failed calling __setitem__");
    }

  }

  PyGILState_Release(gstate);
  GYOTO_DEBUG << "done.\n";
}
