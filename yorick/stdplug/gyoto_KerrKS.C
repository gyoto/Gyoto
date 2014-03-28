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

#include <GyotoKerrKS.h>
#include <GyotoFactory.h>
#include "ygyoto.h"
#include "yapi.h"

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Metric;

#define OBJ gg

void ygyoto_KerrKS_eval(SmartPointer<Metric::Generic> *gg_, int argc) {

  static char const * knames[]={
    "unit", "spin",
    YGYOTO_METRIC_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Metric, KerrKS, knames, YGYOTO_METRIC_GENERIC_KW_N+2);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE2(spin);

  YGYOTO_WORKER_CALL_GENERIC(Metric);
  
}


extern "C" {
  void Y__gyoto_KerrKS_register_as_Metric(){
    ygyoto_Metric_register("KerrKS",&ygyoto_KerrKS_eval);
  }

  void
  Y_gyoto_KerrKS(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Metric, Metric::Generic, KerrKS, metric);
    if ((*OBJ)->kind() != "KerrKS")
      y_error("Expecting Metric of kind KerrKS");
    ygyoto_KerrKS_eval(gg, argc);
  }

}
