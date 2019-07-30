/*
    Copyright 2019 Thibaut Paumard

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
  wchar_t * wargv[argc];
  size_t sz;
  for (int k=0; k<argc; ++k) {
    wargv[k] = Py_DecodeLocale(argv[k], &sz);
  }
  GYOTO_DEBUG << " setting argv" << std::endl;
  PySys_SetArgv(argc, wargv);
  GYOTO_DEBUG << " done" << std::endl;
  std::string code =
    "import gyoto.animate\n"
    "gyoto.animate.main()\n";
  GYOTO_DEBUG << "trying to run Python code: " << std::endl << code ;
  PyRun_SimpleString(code.c_str());
  GYOTO_DEBUG << "back to mk_video" << std::endl;
  for (int k=0; k<argc; ++k) {
    PyMem_RawFree(wargv[k]);
  }
  return 0;
}
