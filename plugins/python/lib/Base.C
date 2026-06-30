/*
    Copyright © 2015-2017, 2019, 2022 Thibaut Paumard

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
#include "GyotoValue.h"
#include "GyotoProperty.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL GyotoPython_ARRAY_API
#include <numpy/arrayobject.h>

using namespace Gyoto;
using namespace Gyoto::Python;
using namespace std;

class Proxy:
  public Gyoto::Object {
public:
  void set_double(double);
};


PyObject * Gyoto::Python::PyObject_FromGyotoValue(const Gyoto::Value& val){
  PyObject * pVal = NULL;
  switch (val.type) {
  case Property::double_t:
    pVal = PyFloat_FromDouble(val);
    break;
  case Property::long_t:
    pVal = PyLong_FromLong(val);
    break;
  case Property::unsigned_long_t:
    pVal = PyLong_FromUnsignedLong(val);
    break;
  case Property::size_t_t:
    pVal = PyLong_FromSize_t(val);
    break;
  case Property::bool_t:
    pVal = PyBool_FromLong(val);
    break;
  case Property::filename_t:
  case Property::string_t:
    pVal = PyUnicode_FromString(std::string(val).c_str());
    break;
  case Property::vector_double_t:
    {
      std::vector<double> vval = val;
      npy_intp dims[] = {npy_intp(vval.size())};

      pVal = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      for (npy_intp k=0; k<dims[0]; ++k) *(double*)PyArray_GetPtr((PyArrayObject*)pVal, &k)=vval[k];
    }
    break;
  case Property::vector_unsigned_long_t:
    {
      std::vector<unsigned long> vval = val;
      npy_intp dims[] = {npy_intp(vval.size())};

      pVal = PyArray_SimpleNew(1, dims, NPY_ULONG);
      for (npy_intp k=0; k<dims[0]; ++k) *(unsigned long*)PyArray_GetPtr((PyArrayObject*)pVal, &k)=vval[k];
    }
    break;
  case Property::metric_t:
    {
      GYOTO_DEBUG_EXPR(val.type);
      pVal = PyObject_CallFunction(pGyotoMetric(), "l", (long)(Gyoto::Metric::Generic*)Gyoto::SmartPointer<Gyoto::Metric::Generic>(val));
    }
    break;
  case Property::screen_t:
    {
      GYOTO_DEBUG_EXPR(val.type);
      pVal = PyObject_CallFunction(pGyotoScreen(), "l", (long)(Gyoto::Screen*)Gyoto::SmartPointer<Gyoto::Screen>(val));
    }
    break;
  case Property::astrobj_t:
    {
      GYOTO_DEBUG_EXPR(val.type);
      pVal = PyObject_CallFunction(pGyotoAstrobj(), "l", (long)(Gyoto::Astrobj::Generic*)Gyoto::SmartPointer<Gyoto::Astrobj::Generic>(val));
    }
    break;
  case Property::spectrum_t:
    {
      GYOTO_DEBUG_EXPR(val.type);
      pVal = PyObject_CallFunction(pGyotoSpectrum(), "l", (long)(Gyoto::Spectrum::Generic*)Gyoto::SmartPointer<Gyoto::Spectrum::Generic>(val));
    }
    break;
  case Property::spectrometer_t:
    {
      GYOTO_DEBUG_EXPR(val.type);
      pVal = PyObject_CallFunction(pGyotoSpectrometer(), "l", (long)(Gyoto::Spectrometer::Generic*)Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>(val));
    }
    break;
  case Property::empty_t:
    pVal = Py_None;
    break;
  default:
    GYOTO_ERROR("Type not implemented in Python::Metric::PyObject_FromGyotoValue()");
  }
  return pVal;
}

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

  if (!pGetArgSpec) { // search for getfullargspec
    PyObject * pName = PyUnicode_FromString("inspect");
    PyObject * pModule = PyImport_Import(pName);
    Py_XDECREF(pName); pName=NULL;
    pGetArgSpec = PyObject_GetAttrString(pModule, "getfullargspec");
  }

  if (!pGetArgSpec) { // getargspec deprecated since Python 3.0, removed in 3.11
    PyObject * pName = PyUnicode_FromString("inspect");
    PyObject * pModule = PyImport_Import(pName);
    Py_XDECREF(pName); pName=NULL;
    pGetArgSpec = PyObject_GetAttrString(pModule, "getargspec");
  }

  if (!pGetArgSpec) {
    PyErr_Print();
    GYOTO_ERROR("Failed finding method getargspec or getfullargspec in module inspect");
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

PyObject * Gyoto::Python::pGyotoScreen() {
  static PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "Screen");
  }
  return res;
}

PyObject * Gyoto::Python::pGyotoAstrobj() {
  static PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "Astrobj");
  }
  return res;
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

PyObject * Gyoto::Python::pGyotoSpectrometer() {
  static PyObject * res = NULL;
  static bool need_load = true;
  if (need_load) {
    need_load=false;
    PyObject* pGyoto=Gyoto::Python::PyImport_Gyoto();
    if (pGyoto)
      res = PyObject_GetAttrString(pGyoto, "Spectrometer");
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
    PyErr_Print();
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
    pModule_(NULL), pInstance_(NULL),
    pProperties_(NULL), pSet_(NULL), pGet_(NULL)
{}

Base::Base(const Base& o)
: module_(o.module_), inline_module_(o.inline_module_),
  class_(o.class_), parameters_(o.parameters_),
  pModule_(o.pModule_), pInstance_(o.pInstance_),
  pProperties_(o.pProperties_), pSet_(o.pSet_), pGet_(o.pGet_)
{
  Py_XINCREF(pModule_);
  Py_XINCREF(pInstance_);
  Py_XINCREF(pProperties_);
  Py_XINCREF(pSet_);
  Py_XINCREF(pGet_);
}

Base::~Base() {
  Py_XDECREF(pGet_);
  Py_XDECREF(pSet_);
  Py_XDECREF(pProperties_);
  Py_XDECREF(pInstance_);
  Py_XDECREF(pModule_);
}


/* Accessors */

std::string Base::module() const { return module_; }
void Base::module(const std::string &m) {
  GYOTO_DEBUG << "Loading Python module " << m << endl;
  module_=m;
  inline_module_="";
  detachInstance();
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  Py_XDECREF(pModule_); pModule_=nullptr;
  if (m=="") {
    class_ = "";
    return;
  }

  PyObject *pName=PyUnicode_FromString(m.c_str());
  PyObject *tmp = nullptr;
  guardian.track(tmp);

  if (!pName) {
    PyErr_Print();
    GYOTO_ERROR("Failed translating string to Python");
  }
  pModule_ = PyImport_Import(pName);
  Py_DECREF(pName);
  if (PyErr_Occurred() || !pModule_) {
    PyErr_Print();
    Py_XDECREF(pModule_); pModule_=nullptr;
    GYOTO_ERROR("Failed loading Python module");
  }

  // if class_ is the empty string, check whether there is a single
  // class in module_. In this case set class_ to its name.
  checkModuleForSingleClass();

  if (class_ != "") {
    GYOTO_DEBUG << "class_ is set already (to '" << class_ << "')."
      "Instantiating it now." << std::endl;
    klass(class_);
  }
  GYOTO_DEBUG << "Done loading Python module " << m << endl;
}

std::string Base::inlineModule() const { return inline_module_; }
void Base::inlineModule(const std::string &m) {
  inline_module_=m;
  module_="";
  detachInstance();

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  Py_XDECREF(pModule_);

  if (m=="") {
    class_="";
    return;
  }

  GYOTO_DEBUG << "Loading inline Python module :" << m << endl;
  pModule_ = Gyoto::Python::PyModule_NewFromPythonCode(m.c_str());
  if (PyErr_Occurred() || !pModule_) {
    PyErr_Print();
    GYOTO_ERROR("Failed loading inline Python module");
  }

  // if class_ is the empty string, check whether there is a single
  // class in module_. In this case set class_ to its name.
  checkModuleForSingleClass();

  if (class_ != "") klass(class_);
  GYOTO_DEBUG << "Done loading Python module " << m << endl;
}

void Base::checkModuleForSingleClass() {
  if (!pModule_ || class_ != "") return;
  GYOTO_DEBUG << "class_ is empty: check whether there is a single class in module...\n";
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;

  PyObject * dict = PyModule_GetDict(pModule_);
  PyObject *key, *value, *tmp=nullptr;
  Py_ssize_t pos = 0, nclass=0;
  guardian.track(tmp);

  // Loop on all objects in module_. Store in class_ the name of the
  // first class found. As soon as 2 classes have been found, reset
  // class_ to "" and break out. Therefore, the net effect is that
  // class_ is set to the name of the (first) class in the module if
  // and only if there is a single class in the module.
  while (PyDict_Next(dict, &pos, &key, &value)) {
    if (
#if PY_VERSION_HEX < 0x03000000
	PyClass_Check(value) ||
#endif
	PyObject_TypeCheck(value, &PyType_Type)) {
      ++nclass;
      if (nclass >= 2) {
	class_ = "";
	GYOTO_DEBUG << "several classes in module" << endl;
	return;
      }

      if (PyUnicode_Check(key)) {
	tmp = PyUnicode_AsUTF8String(key);
      } else {
	tmp = key;
	Py_INCREF(tmp);
      }

      if (!PyBytes_Check(tmp))
	GYOTO_ERROR("not a PyBytes string");

      class_= PyBytes_AsString(tmp);

      Py_DECREF(tmp); tmp=nullptr;
    }
  }
  if (nclass == 1) {
    GYOTO_DEBUG << "single class in module: " << class_ << endl;
  } else if (nclass == 0) {
    GYOTO_ERROR("no class in module\n");
  }
}

void Base::detachInstance() {
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  Py_XDECREF(pProperties_); pProperties_=NULL;
  Py_XDECREF(pGet_); pGet_=NULL;
  Py_XDECREF(pSet_); pSet_=NULL;
  Py_XDECREF(pInstance_); pInstance_=NULL;
}

PyObject * Base::instantiateClass(std::string &klass) const {
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyObject * pClass = nullptr;
  guardian.track(pClass);

  if (klass == "") GYOTO_ERROR ("klass is the empty string");

  pClass = PyObject_GetAttrString(pModule_, class_.c_str());
  if (PyErr_Occurred() || !pClass) {
    PyErr_Print();
    GYOTO_ERROR("Could not find class in module");
  }
  if (!PyCallable_Check(pClass)) {
    GYOTO_ERROR("Class is not callable");
  }

  PyObject * retval = PyObject_CallObject(pClass, NULL);
  // will be done automaticallu upon destruction of guardian
  //Py_DECREF(pClass); pClass=NULL;
  if (PyErr_Occurred() || !retval) {
    PyErr_Print();
    Py_XDECREF(retval);
    GYOTO_ERROR("Could not instantiate class");
  }

  return retval;
}


void Base::attachInstance (PyObject * instance) {
  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyTypeObject * type = nullptr;
  PyObject * inspect_module = nullptr;
  guardian.track(inspect_module);
  PyObject * getsource_func = nullptr;
  guardian.track(getsource_func);
  PyObject * source_str = nullptr;
  guardian.track(source_str);

  // in any case, free the previous instance
  Py_XDECREF(pInstance_); pInstance_=nullptr;

  // if instance is null, return
  if (instance == nullptr)
    return;

  // else check that instance is (likely) a class instance
  if (PyType_Check(instance) || PyModule_Check(instance) ||
      PyFunction_Check(instance) || PyCFunction_Check(instance) ||
      !PyObject_TypeCheck(instance, &PyBaseObject_Type))
    GYOTO_ERROR("address doesn't seem to point to a class instance");

  pInstance_ = instance;
  Py_INCREF(pInstance_);

  // is class_ is empty, set it to the instance's class
  if (class_=="") {
    // Get the type object of the instance
    type = Py_TYPE(instance);
    // Get the tp_name field (const char*)
    class_ = (type && type->tp_name) ? type->tp_name : "UNKNOWN";

    // if pModule_ is null, get the module that contains the class
    if ( module_ == "" && inline_module_ == "" && (type = Py_TYPE(instance))) {
      PyObject* module_name_obj = PyObject_GetAttrString((PyObject*)type, "__module__");
      if (!module_name_obj) {
	PyErr_Print();
	GYOTO_ERROR("Error getting module name");
      }

      // Ensure it's a string
      if (!PyUnicode_Check(module_name_obj)) {
	Py_DECREF(module_name_obj);
	PyErr_Print();
	GYOTO_ERROR("Module name is not a string");
      }

      // Convert to UTF-8
      const char* module_name = PyUnicode_AsUTF8(module_name_obj);
      if (!module_name) {
        Py_DECREF(module_name_obj);
	PyErr_Print();
	GYOTO_ERROR("Could not convert module name to UTF8");
      }

      // set module_
      module_  = module_name;

      // set pModule_
      if (module_ == "__main__") {
	// Try to store source code in InlineModule
	// Import the inspect module
	inspect_module = PyImport_ImportModule("inspect");
	if (!inspect_module) {
	  PyErr_Print();
	  Py_CLEAR(pInstance_);
	  GYOTO_ERROR("failed importing inspect module");
	}
	// Get inspect.getsource
	getsource_func = PyObject_GetAttrString(inspect_module, "getsource");
	if (!getsource_func || PyErr_Occurred()) {
	  PyErr_Print();
	  Py_CLEAR(pInstance_);
	  GYOTO_ERROR("failed getting function inpect.getsource");
	}
	// Call inspect.getsource.  If source cannot be retrieved or
	// converted tu UTF8, module_ remains set to __main__
	source_str = PyObject_CallFunction(getsource_func, "O", type);
	if (PyErr_Occurred() || ! source_str) {
	  // ignore this error
	  PyErr_Clear();
	} else {
	  // successfully got source: store is as inline_module_
	  const char* source_utf8 = PyUnicode_AsUTF8(source_str);
	  if (PyErr_Occurred() || ! source_utf8) {
	    // ignore this error
	    PyErr_Clear();
	  } else {
	    inline_module_ = source_utf8;
	    module_ = "";
	  }
	}
	// in any case, load the actual module in pModule_
	// this is a bit special for __main__
	PyObject* sys_modules = PyImport_GetModuleDict();
	pModule_ = PyDict_GetItemString(sys_modules, "__main__");
	if (PyErr_Occurred() || !pModule_) {
	  PyErr_Print();
	  Py_CLEAR(pInstance_);
	  GYOTO_ERROR("Error getting __main__ module");
	}
	Py_INCREF(pModule_);
      } else {
	// set pModule_
	pModule_ = PyImport_ImportModule(module_name);
	if (PyErr_Occurred() || !pModule_) {
	  PyErr_Print();
	  Py_CLEAR(pInstance_);
	  GYOTO_ERROR("Error getting instance class module");
	}
      }
      Py_DECREF(module_name_obj);
    }
  }

  pSet_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "set");

  if (PyErr_Occurred()) {
    PyErr_Print();
    Py_XDECREF(pSet_); pSet_=NULL;
    Py_XDECREF(pInstance_); pInstance_=NULL;
    GYOTO_ERROR("Error getting set method");
  }

  pGet_=
    Gyoto::Python::PyInstance_GetMethod(pInstance_, "get");

  if (PyErr_Occurred()) {
    PyErr_Print();
    Py_XDECREF(pSet_); pSet_=NULL;
    Py_XDECREF(pGet_); pGet_=NULL;
    Py_XDECREF(pInstance_); pInstance_=NULL;
    GYOTO_ERROR("Error getting set method");
  }

  pProperties_ = NULL;
  if (PyObject_HasAttrString(pInstance_, "properties"))
    pProperties_ = PyObject_GetAttrString(pInstance_, "properties");
  if (PyErr_Occurred()) {
    PyErr_Print();
    Py_XDECREF(pProperties_); pInstance_=NULL;
    Py_XDECREF(pGet_); pGet_=NULL;
    Py_XDECREF(pSet_); pSet_=NULL;
    Py_XDECREF(pInstance_); pInstance_=NULL;
    GYOTO_ERROR("Error getting properties member");
  }
  if (pProperties_) {
    if (!PyDict_Check(pProperties_)) {
      Py_XDECREF(pProperties_); pInstance_=NULL;
      Py_XDECREF(pGet_); pGet_=NULL;
      Py_XDECREF(pSet_); pSet_=NULL;
      Py_XDECREF(pInstance_); pInstance_=NULL;
      GYOTO_ERROR("'properties' is not a dict'");
    }
  }
}

void Base::instance(PyObject * instance) {
  detachInstance();
  Base::module("");
  class_="";
  if (!instance) return;
  attachInstance(instance);
}

std::string Base::klass() const { return class_; }
void Base::klass(const std::string &f) {
  Gyoto::Python::GILGuard guardian;
  PyObject * tmpinstance = nullptr;
  guardian.track(tmpinstance);

  // always start by storing the value
  class_=f;

  // always detach the previously attached instance and methods
  detachInstance();

  // if the string is empty, that's a signal to also free the module
  if (class_ == "") {
    module_ = "";
    inline_module_ = "";
    Py_XDECREF(pModule_); pModule_ = nullptr;
  }
  if (!pModule_) return;

  GYOTO_DEBUG << "Instantiating Python class " << f << endl;

  // actually instantiate the class.
  tmpinstance = instantiateClass(class_);
  if (!tmpinstance)
    GYOTO_ERROR("Failed instantiating Python class");

  // attach the instance
  attachInstance(tmpinstance);

  // will be done when guardian is destroyed
  // Py_XDECREF(tmpinstance); tmpinstance=NULL;

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

bool Base::hasPythonProperty(std::string const &key) const {
  if (!pProperties_) return false;
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject * pKey = PyUnicode_FromString(key.c_str());

  GYOTO_DEBUG_EXPR(key);
  GYOTO_DEBUG_EXPR(pKey);
  GYOTO_DEBUG_EXPR(pProperties_);

  int key_in_props = PyDict_Contains(pProperties_, pKey);

  Py_XDECREF(pKey);
  PyGILState_Release(gstate);

  GYOTO_DEBUG_EXPR(key_in_props);

  if (key_in_props==-1) GYOTO_ERROR("Error checking for key in Python properties");
  return key_in_props;

}

Property::type_e Base::pythonPropertyType(std::string const &key) const {
  GYOTO_DEBUG_EXPR(key);
  if (!pProperties_) GYOTO_ERROR("no properties");
  if (!hasPythonProperty(key)) GYOTO_ERROR("no such property");

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyObject * pKey = PyUnicode_FromString(key.c_str());
  guardian.track(pKey);

  GYOTO_DEBUG_EXPR(pKey);

  if (!pKey) {
    GYOTO_ERROR("Could not convert '"+key+"' to PyUnicode string");
  }

  GYOTO_DEBUG_EXPR(pProperties_);

  PyObject * pType = PyDict_GetItem(pProperties_, pKey);
  if (PyDict_Check(pType)) {
    GYOTO_DEBUG << "property[" << key <<"] is a dict" << std::endl;
    pType = PyDict_GetItemString(pType, "type");
    if (!pType) { GYOTO_ERROR("Please specify type for property: " + key); }
  }

  if (!PyUnicode_Check(pType)) {
    GYOTO_ERROR("Type of '"+key+"' is not a UNICODE string");
  }

  std::string stype=PyUnicode_AsUTF8(pType);
  GYOTO_DEBUG_EXPR(stype);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("could not convert type to UTF8");
  }

  return Property::typeFromString(stype);
}

void Base::setPythonProperty(std::string const &key, Value val) {

  if (!pSet_) GYOTO_ERROR("self(self, key, val) method not implemented");

  GYOTO_DEBUG_EXPR(key);
  GYOTO_DEBUG_EXPR(val.type);

  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject * pKey = PyUnicode_FromString(key.c_str());

  GYOTO_DEBUG_EXPR(pKey);
  GYOTO_DEBUG_EXPR(pProperties_);

  PyObject * pVal = PyObject_FromGyotoValue(val);

  if (PyErr_Occurred()) {
    Py_XDECREF(pKey);
    Py_XDECREF(pVal);
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred while setting property");
  }

  PyObject * pR =
    PyObject_CallFunctionObjArgs(pSet_, pKey, pVal, NULL);

    Py_XDECREF(pR);
    Py_XDECREF(pKey);
    Py_XDECREF(pVal);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyGILState_Release(gstate);
    GYOTO_ERROR("Error occurred while setting property");
  }

  PyGILState_Release(gstate);
}

void Base::setPythonProperty(std::string const &key, Value val, std::string const &unit) {

  if (!pSet_) GYOTO_ERROR("self(self, key, val, unit) method not implemented");

  GYOTO_DEBUG_EXPR(key);
  GYOTO_DEBUG_EXPR(val.type);
  GYOTO_DEBUG_EXPR(unit);

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyObject * pKey = PyUnicode_FromString(key.c_str());
  guardian.track(pKey);

  GYOTO_DEBUG_EXPR(pKey);
  GYOTO_DEBUG_EXPR(pProperties_);

  PyObject * pVal = PyObject_FromGyotoValue(val);
  guardian.track(pVal);

  PyObject * pUnit = PyUnicode_FromString(unit.c_str());
  guardian.track(pUnit);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error occurred while setting property");
  }

  PyObject * pR =
    PyObject_CallFunctionObjArgs(pSet_, pKey, pVal, pUnit, NULL);

  Py_XDECREF(pR);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error occurred while setting property");
  }

}

Value Base::getPythonProperty(std::string const &key) const {
  return getPythonProperty(key, "");
}

Value Base::getPythonProperty(std::string const &key, std::string const &unit) const {
  GYOTO_DEBUG_EXPR(key);
  if (!pProperties_) GYOTO_ERROR("no properties");
  if (!hasPythonProperty(key)) GYOTO_ERROR("no such property");
  if (!pGet_) GYOTO_ERROR("get(self, key) method not implemented");

  [[maybe_unused]] Gyoto::Python::GILGuard guardian;
  PyObject * pKey = PyUnicode_FromString(key.c_str());
  guardian.track(pKey);

  GYOTO_DEBUG_EXPR(pKey);
  if (!pKey) {
    GYOTO_ERROR("Could not convert '"+key+"' to PyUnicode string");
  }

  PyObject * pUnit = PyUnicode_FromString(unit.c_str());
  guardian.track(pUnit);

  GYOTO_DEBUG_EXPR(pUnit);
  if (!pUnit) {
    GYOTO_ERROR("Could not convert '"+unit+"' to PyUnicode string");
  }

  GYOTO_DEBUG_EXPR(pProperties_);

  PyObject * pVal = nullptr;
  if (unit=="") pVal = PyObject_CallFunctionObjArgs(pGet_, pKey, NULL);
  else          pVal = PyObject_CallFunctionObjArgs(pGet_, pKey, pUnit, NULL);

  guardian.track(pVal);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error occurred while calling get() in getPythonProperty()");
  }

  Value val;
  Property::type_e type=pythonPropertyType(key);

  GYOTO_DEBUG_EXPR(type);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error getting property type");
  }

  switch (type) {
  case Property::double_t:
    val = PyFloat_AsDouble(pVal);
    break;
  case Property::long_t:
    val = PyLong_AsLong(pVal);
    break;
  case Property::unsigned_long_t:
    val = PyLong_AsUnsignedLong(pVal);
    break;
  case Property::size_t_t:
    val = PyLong_AsSize_t(pVal);
    break;
  case Property::bool_t:
    {
      int result = PyObject_IsTrue(pVal);
      if (result == -1 || PyErr_Occurred()) {
        PyErr_Print();
	GYOTO_ERROR("this value can't be converted to boolean");
      }
      val =  bool(result);
    }
    break;
  case Property::string_t:
  case Property::filename_t:
    {
      if (!PyUnicode_Check(pVal))
	GYOTO_ERROR("this doesn't look like a string");

      const char* cstring = PyUnicode_AsUTF8(pVal);
      if (!cstring || PyErr_Occurred()) {
	PyErr_Print();
	GYOTO_ERROR("error converting Python string to C string");
      }
      std::string sstring = cstring;
      val = sstring;
    }
    break;
  case Property::vector_double_t:
    {
    PyArray_Descr* pd = PyArray_DescrFromType(NPY_DOUBLE);
    PyObject * pArr = PyArray_FromAny(pVal, pd, 0, 0, NPY_ARRAY_CARRAY, NULL);
    //    Py_XDECREF(pd); // PyArray_FromAny steals a reference to *pd

    if (PyErr_Occurred()) {
      Py_XDECREF(pArr);

      PyErr_Print();
      GYOTO_ERROR("Error occurred while calling get() in getPythonProperty()");
    }

    double *buffer=(double*)PyArray_DATA((PyArrayObject*)pArr);
    npy_intp sz = PyArray_Size(pArr);

    std::vector<double> vec(sz);
    for (npy_intp k=0; k<sz; ++k) vec[k]=buffer[k];
    val=vec;

    Py_XDECREF(pArr);
    }
    break;
  case Property::vector_unsigned_long_t:
    {
    PyArray_Descr* pd = PyArray_DescrFromType(NPY_ULONG);
    PyObject * pArr = PyArray_FromAny(pVal, pd, 0, 0, NPY_ARRAY_CARRAY, NULL);
    //    Py_XDECREF(pd); // PyArray_FromAny steals a reference to *pd

    if (PyErr_Occurred()) {
      Py_XDECREF(pArr);

      PyErr_Print();
      GYOTO_ERROR("Error occurred while calling get() in getPythonProperty()");
    }

    unsigned long *buffer=(unsigned long*)PyArray_DATA((PyArrayObject*)pArr);
    npy_intp sz = PyArray_Size(pArr);

    std::vector<unsigned long> vec(sz);
    for (npy_intp k=0; k<sz; ++k) vec[k]=buffer[k];
    val=vec;

    Py_XDECREF(pArr);
    }
    break;
  case Property::metric_t:
    {
    PyObject * pAddr = PyObject_CallMethod(pVal, "getPointer", NULL);
    if (PyErr_Occurred()) {
      Py_XDECREF(pAddr);
      PyErr_Print();
      GYOTO_ERROR("Error occurred in getPythonProperty()");
    }
    long address = PyLong_AsLong(pAddr);
    Py_CLEAR(pAddr);

    val = Gyoto::SmartPointer<Gyoto::Metric::Generic>((Gyoto::Metric::Generic*) (address));
    }
    break;
  case Property::screen_t:
    {
    PyObject * pAddr = PyObject_CallMethod(pVal, "getPointer", NULL);
    if (PyErr_Occurred()) {
      Py_XDECREF(pAddr);
      PyErr_Print();
      GYOTO_ERROR("Error occurred in getPythonProperty()");
    }
    long address = PyLong_AsLong(pAddr);
    Py_CLEAR(pAddr);

    val = Gyoto::SmartPointer<Gyoto::Screen>((Gyoto::Screen*) (address));
    }
    break;
  case Property::astrobj_t:
    {
    PyObject * pAddr = PyObject_CallMethod(pVal, "getPointer", NULL);
    if (PyErr_Occurred()) {
      Py_XDECREF(pAddr);
      PyErr_Print();
      GYOTO_ERROR("Error occurred in getPythonProperty()");
    }
    long address = PyLong_AsLong(pAddr);
    Py_CLEAR(pAddr);

    val = Gyoto::SmartPointer<Gyoto::Astrobj::Generic>((Gyoto::Astrobj::Generic*) (address));
    }
    break;
  case Property::spectrum_t:
    {
    PyObject * pAddr = PyObject_CallMethod(pVal, "getPointer", NULL);
    if (PyErr_Occurred()) {
      Py_XDECREF(pAddr);
      PyErr_Print();
      GYOTO_ERROR("Error occurred in getPythonProperty()");
    }
    long address = PyLong_AsLong(pAddr);
    Py_CLEAR(pAddr);

    val = Gyoto::SmartPointer<Gyoto::Spectrum::Generic>((Gyoto::Spectrum::Generic*) (address));
    }
    break;
  case Property::spectrometer_t:
    {
    PyObject * pAddr = PyObject_CallMethod(pVal, "getPointer", NULL);
    if (PyErr_Occurred()) {
      Py_XDECREF(pAddr);
      PyErr_Print();
      GYOTO_ERROR("Error occurred in getPythonProperty()");
    }
    long address = PyLong_AsLong(pAddr);
    Py_CLEAR(pAddr);

    val = Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>((Gyoto::Spectrometer::Generic*) (address));
    }
    break;
  case Property::empty_t:
    GYOTO_ERROR("cant' get an empty property");
    val = false;
    break;
  default:
    GYOTO_ERROR("unknown property type");
  }

  Py_CLEAR(pVal);

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Error occurred in getPythonProperty()");
  }

  GYOTO_DEBUG << "exiting" << std::endl;
  return val;

}

// GILGuard class
Gyoto::Python::GILGuard::GILGuard() : gstate_(PyGILState_Ensure()) {}

Gyoto::Python::GILGuard::~GILGuard() {
  // Release the GIL
  PyGILState_Release(gstate_);

  // Decrement references for all tracked PyObject* pointers
  for (PyObject** obj_ptr : tracked_objects_) {
    if (obj_ptr && *obj_ptr) {
      Py_XDECREF(*obj_ptr);
      *obj_ptr = nullptr;  // Set to nullptr to avoid dangling pointers
    }
  }
}

void Gyoto::Python::GILGuard::track(PyObject*& obj) {
  tracked_objects_.push_back(&obj);
}
