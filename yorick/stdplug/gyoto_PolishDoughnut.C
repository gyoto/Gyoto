/*
    Copyright (c) 2012-2013 Thibaut Paumard

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

#include "../ygyoto.h"
#include "yapi.h"
#include <GyotoPolishDoughnut.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Astrobj;

// on_eval worker
void ygyoto_PolishDoughnut_eval(SmartPointer<Astrobj::Generic>* OBJ_, int argc) {

  static char const * knames[]={
    "unit",
    "lambda", "centralenthalpyperunitvolume", "centraltemperature", "beta",
    "spectraloversampling",
    "l0", "Wsurface", "Wcentre", "rcusp", "rcentre",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Astrobj, PolishDoughnut, knames,
		     YGYOTO_ASTROBJ_GENERIC_KW_N+12);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE2(lambda);
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(centralEnthalpyPerUnitVolume);
  YGYOTO_WORKER_GETSET_DOUBLE2(centralTemp);
  YGYOTO_WORKER_GETSET_DOUBLE2(beta);
  YGYOTO_WORKER_GETSET_LONG2(spectralOversampling);
  YGYOTO_WORKER_GET_DOUBLE(getL0);
  YGYOTO_WORKER_GET_DOUBLE(getWsurface);
  YGYOTO_WORKER_GET_DOUBLE(getWcentre);
  YGYOTO_WORKER_GET_DOUBLE(getRcusp);
  YGYOTO_WORKER_GET_DOUBLE(getRcentre);

  YGYOTO_WORKER_CALL_GENERIC(Astrobj);
}


extern "C" {
  void Y__gyoto_PolishDoughnut_register_as_Astrobj(){
    ygyoto_Astrobj_register("PolishDoughnut",&ygyoto_PolishDoughnut_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_PolishDoughnut(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Astrobj, Astrobj::Generic, PolishDoughnut, astrobj);
    if ((*OBJ)->kind()!="PolishDoughnut")
      y_error("Expecting Astrobj of kind PolishDoughnut");
    ygyoto_PolishDoughnut_eval(OBJ, argc);
  }

}
