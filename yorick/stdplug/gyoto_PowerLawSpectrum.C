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

#include "ygyoto.h"
#include "yapi.h"
#include <iostream>
#include "GyotoFactory.h"
#include "GyotoPowerLawSpectrum.h"

namespace YGyoto {
  namespace Spectrum {
    ygyoto_Spectrum_eval_worker_t PowerLawYEval;
  }
}

using namespace Gyoto;
using namespace YGyoto;
using namespace Gyoto::Spectrum;
using namespace YGyoto::Spectrum;
using namespace std;

void YGyoto::Spectrum::PowerLawYEval(SmartPointer<Generic> * OBJ_, int argc) {

  static char const * knames[]={
    "unit",
    "constant", "exponent",
    YGYOTO_SPECTRUM_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Spectrum, PowerLaw,
		     knames, YGYOTO_SPECTRUM_GENERIC_KW_N+3);

  YGYOTO_WORKER_SET_UNIT;

  YGYOTO_WORKER_GETSET_DOUBLE2(constant);
  YGYOTO_WORKER_GETSET_DOUBLE2(exponent);
  
  YGYOTO_WORKER_CALL_GENERIC(Spectrum);

}

extern "C" {
  void
  Y_gyoto_PowerLawSpectrum(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Spectrum, Gyoto::Spectrum::Generic, PowerLaw, spectrum);
    if ((*OBJ)->kind().compare("PowerLaw"))
      y_error("Expecting Spectrum of kind PowerLaw");
    PowerLawYEval(OBJ, argc);
  }

  void Y__gyoto_PowerLawSpectrum_register_as_Metric(){
    ygyoto_Spectrum_register("PowerLaw",&PowerLawYEval);
  }
}
