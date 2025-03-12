/*
    Copyright 2019, 2025 Thibaut Paumard

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
#include "GyotoProperty.h"
#include "GyotoError.h"

#include <cstdlib>

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL GyotoPython_ARRAY_API
#include <numpy/arrayobject.h>

#include <iostream>
extern "C" int mk_video(int argc, char ** argv) {
  GYOTO_DEBUG << " in mk_video()" << std::endl;

  PyStatus status;
  PyConfig config;
  std::string code =
    "import gyoto.animate\n"
    "gyoto.animate.main()\n";
  int retval;

  GYOTO_DEBUG << "initializing Python" << std::endl ;
  PyConfig_InitIsolatedConfig(&config);
  status = PyConfig_SetBytesArgv(&config, argc, argv);
  if (PyStatus_Exception(status)) goto exception;
  status = Py_InitializeFromConfig(&config);
  if (PyStatus_Exception(status)) goto exception;
  PyConfig_Clear(&config);

  GYOTO_DEBUG << "running Python code" << std::endl << code ;
  retval = PyRun_SimpleString(code.c_str());
  GYOTO_DEBUG << "exiting mk_video" << std::endl;

  return retval;

exception:
  GYOTO_DEBUG << "something went wrong" << std::endl ;
  PyConfig_Clear(&config);
  if (PyStatus_IsError(status)) GYOTO_ERROR(status.err_msg);
  return status.exitcode;
}
