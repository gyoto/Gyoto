/*
    Copyright 2015, 2022 Thibaut Paumard

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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL GyotoPython_ARRAY_API
#include <numpy/arrayobject.h>

using namespace Gyoto;

#if PY_VERSION_HEX < 0x03070000
static PyThreadState* mainPyThread=NULL;
#endif

namespace Gyoto {
  // import_array is actually a MACRO which returns a value.
  // We want to eat this return.
  bool eat_import_array() { import_array1(false); return true;}
}

extern "C" void __GyotoPluginInit() {
  Spectrum::Register("Python",
		     &(Spectrum::Subcontractor<Spectrum::Python>));
  Metric::Register("Python",
		     &(Metric::Subcontractor<Metric::Python>));
  Astrobj::Register("Python::Standard",
		    &(Astrobj::Subcontractor<Astrobj::Python::Standard>));
  Astrobj::Register("Python::ThinDisk",
		    &(Astrobj::Subcontractor<Astrobj::Python::ThinDisk>));

  Py_InitializeEx(0);

  PyObject *pSys = PyImport_ImportModule("sys");
  PyObject *pPath = PyObject_GetAttrString(pSys, "path");
  PyObject *pDir = PyUnicode_FromString(".");
  Py_XDECREF(pSys);
  PyList_Reverse(pPath);
  PyList_Append(pPath, pDir);
  Py_XDECREF(pDir);
  PyList_Reverse(pPath);
  Py_XDECREF(pPath);

  Py_XDECREF(PyImport_ImportModule("numpy"));
  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Failed importing numpy");
  }
  Gyoto::eat_import_array();

# if PY_VERSION_HEX < 0x03070000
  if (!PyEval_ThreadsInitialized()) {
    PyEval_InitThreads();
    mainPyThread = PyEval_SaveThread();
  }
# endif

  if (PyErr_Occurred()) {
    PyErr_Print();
    GYOTO_ERROR("Failed");
  }
}
