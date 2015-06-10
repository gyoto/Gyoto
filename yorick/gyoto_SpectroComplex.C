/*
    Copyright 2013 Thibaut Paumard

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

#include <GyotoComplexSpectrometer.h>
#include <GyotoFactory.h>
#include "ygyoto.h"
#include "yapi.h"

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Spectrometer;

namespace YGyoto {
  ygyoto_Spectrometer_eval_worker_t SpCplxEval;
}

// on_eval worker
void YGyoto::SpCplxEval(SmartPointer<Spectrometer::Generic> *OBJ_, int argc) {
  static char const * knames[]={
    "unit",
    "append", "remove",	"cardinal",		\
    YGYOTO_SPECTROMETER_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Spectrometer, Complex,
		     knames, YGYOTO_SPECTROMETER_GENERIC_KW_N+4);

  YGYOTO_WORKER_SET_UNIT;

  YGYOTO_WORKER_RUN( append(*yget_Spectrometer(iarg)) );
  YGYOTO_WORKER_RUN( remove(ygets_l(iarg)) );

  /* cardinal */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if ((*rvset)++) y_error(rmsg);
    GYOTO_DEBUG << "get cardinal" << endl;
    ypush_long((*OBJ)->getCardinal());
  }
  
  YGYOTO_WORKER_CALL_GENERIC(Spectrometer);
  
  if (!(*paUsed) && (iarg=piargs[0])>=0 && !yarg_nil(iarg)) {
    iarg+=*rvset;
    if ((*rvset)++) y_error(rmsg);
    ++(*paUsed);
    size_t cardinal=(*OBJ) -> getCardinal();
    size_t item = ygets_l(iarg);
    if (item > cardinal) y_error("index overreach array bounds");
    if (item <= 0) {
      item += cardinal;
      if (item <= 0) y_error("index overreach array bounds");
    }
    *ypush_Spectrometer() = (*(*OBJ))[item-1];
  }

  GYOTO_DEBUG << "done\n";
}


extern "C" {
  void Y__gyoto_SpCplx_register_as_Spectrometer(int argc){
    ygyoto_Spectrometer_register(Complex::Kind,&YGyoto::SpCplxEval);
  }

  void
  Y_gyoto_SpectroComplex(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Spectrometer, Spectrometer::Generic,
			     Complex, spectrometer);
    if ((*OBJ)->kindid() != Complex::Kind)
      y_error("Expecting Spectrometer of kind Complex");
    YGyoto::SpCplxEval(OBJ, argc);
  }

}
