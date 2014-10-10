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

#include <GyotoFixedStar.h>
#include <GyotoFactory.h>
#include "../ygyoto.h"
#include "yapi.h"

using namespace Gyoto;

#include <iostream>
using namespace std;
using namespace Gyoto::Astrobj;

#define OBJ ao

// on_eval worker
void ygyoto_FixedStar_eval(SmartPointer<Astrobj::Generic>* ao_, int argc) {

  static char const * knames[]={
    "unit", "radius", "position", "rotating", "spectrum", "opacity", 
    "deltamaxoverradius", "deltamaxoverdistance",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Astrobj, FixedStar, knames, YGYOTO_ASTROBJ_GENERIC_KW_N+8);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(radius);
  YGYOTO_WORKER_GETSET_VECTOR(Pos, 3);
  YGYOTO_WORKER_GETSET_LONG2(rotating);
  YGYOTO_WORKER_GETSET_OBJECT2(spectrum,Spectrum);
  YGYOTO_WORKER_GETSET_OBJECT2(opacity,Spectrum);
  YGYOTO_WORKER_GETSET_DOUBLE2(deltaMaxOverRadius);
  YGYOTO_WORKER_GETSET_DOUBLE2(deltaMaxOverDistance);

  YGYOTO_WORKER_CALL_GENERIC(Astrobj);
}


extern "C" {
  void Y__gyoto_FixedStar_register_as_Astrobj(){
    ygyoto_Astrobj_register("FixedStar",&ygyoto_FixedStar_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_FixedStar(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Astrobj, Astrobj::Generic, FixedStar, astrobj);
    if ((*ao)->kind().compare("FixedStar"))
      y_error("Expecting Astrobj of kind Star");
    ygyoto_FixedStar_eval(ao, argc);
  }

}
