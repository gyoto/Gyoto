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

#include <GyotoKerrBL.h>
#include <GyotoFactory.h>
#include "ygyoto.h"
#include "yapi.h"

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Metric;

void ygyoto_KerrBL_eval(SmartPointer<Metric::Generic> *OBJ_, int argc) {

  static char const * knames[]={
    "unit",
    "spin", "makecoord",		\
    YGYOTO_METRIC_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Metric, KerrBL, knames, YGYOTO_METRIC_GENERIC_KW_N+3);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE(Spin);

  if (kiargs[++k]>=0) { // makecoord
    if (debug()) cerr << "DEBUG: In ygyoto_KerrBL_eval(): get_coord" << endl;
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    long dims[]={1,8};
    long ntot=1;
    double * coord_=ygeta_d(kiargs[k],&ntot,0);
    if (ntot<7) y_error("YINIT should have >= 7 elements");
    double * cst_=ygeta_d(piargs[0],&ntot,0);
    if (ntot!=4) y_error("CST should have 4 elements");
    yarg_drop(1);
    double * coord=ypush_d(dims);
#   ifdef GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_EXPR(OBJ);
#   endif
    (*OBJ)->MakeCoord(coord_,cst_,coord);
  }
  
  YGYOTO_WORKER_CALL_GENERIC(Metric);
}


extern "C" {
  void Y__gyoto_KerrBL_register_as_Metric(){
    ygyoto_Metric_register("KerrBL",&ygyoto_KerrBL_eval);
  }

  void
  Y_gyoto_KerrBL(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT(Metric, KerrBL);
    if ((*OBJ)->getKind() != "KerrBL")
      y_error("Expecting Metric of kind KerrBL");
    ygyoto_KerrBL_eval(OBJ, argc);
  }

}
