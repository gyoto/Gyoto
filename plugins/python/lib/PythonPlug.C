/*
    Copyright 2011 Thibaut Paumard

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

// include Metric headers
#include "GyotoPython.h"
#include <Python.h>
using namespace Gyoto;

static PyThreadState* mainPyThread=NULL;

extern "C" void __GyotopythonInit() {
  Spectrum::Register("Python",
		     &(Spectrum::Subcontractor<Spectrum::Python>));
  Py_InitializeEx(0);
  if (!PyEval_ThreadsInitialized()) {
    PyEval_InitThreads();
    mainPyThread = PyEval_SaveThread();
  }
}
