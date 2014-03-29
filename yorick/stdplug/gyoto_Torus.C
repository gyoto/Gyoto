/*
    Copyright 2011, 2013 Thibaut Paumard

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

#include <GyotoTorus.h>
#include <GyotoFactory.h>
#include "../ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

// on_eval worker
void ygyoto_Torus_eval(SmartPointer<Astrobj::Generic>* OBJ_, int argc) {

  static char const * knames[]={
    "unit", "largeradius", "smallradius", "spectrum", "opacity",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Astrobj, Torus, knames, YGYOTO_ASTROBJ_GENERIC_KW_N+5);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(largeRadius);
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(smallRadius);
  YGYOTO_WORKER_GETSET_OBJECT2(spectrum,Spectrum);
  YGYOTO_WORKER_GETSET_OBJECT2(opacity,Spectrum);

  YGYOTO_WORKER_CALL_GENERIC(Astrobj);
}


extern "C" {
  void Y__gyoto_Torus_register_as_Astrobj(){
    ygyoto_Astrobj_register("Torus",&ygyoto_Torus_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_Torus(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Astrobj, Astrobj::Generic, Torus, astrobj);
    if ((*OBJ)->kind().compare("Torus"))
      y_error("Expecting Astrobj of kind Torus");
    ygyoto_Torus_eval(OBJ, argc);
  }

}
