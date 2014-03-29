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

#include "ygyoto.h"
#include <yapi.h>
#include <pstdlib.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
using namespace std;
using namespace Gyoto;

YGYOTO_YUSEROBJ(Screen, Screen)

extern "C" {
  void gyoto_Screen_eval(void *obj, int argc) {
    SmartPointer<Screen> *OBJ_ = &((gyoto_Screen*)obj)->smptr;
    static char const * knames[]={
      "unit",
      "metric",
      "time","fov","resolution", "mask", "maskwrite",
      "distance", "dmax", "inclination", "paln", "argument",
      "freqobs",
      "projection", "observerpos",
      "fourvel", "screen1", "screen2", "screen3",
      "spectro",
      "skycoord",  "raycoord",
      "xmlwrite", "clone",
      0
    };
    YGYOTO_WORKER_INIT1(Screen, Screen, knames, 24)

    YGYOTO_WORKER_SET_UNIT;
    YGYOTO_WORKER_GETSET_OBJECT2(metric,Metric);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(Time);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(FieldOfView);
    YGYOTO_WORKER_GETSET_LONG(Resolution);

    /* MASK */
    if ((iarg=kiargs[++k])>=0) { // mask
      GYOTO_DEBUG << "mask=\n";
      iarg+=*rvset;
      if (yarg_nil(iarg)) {
	if ((*rvset)++) y_error(rmsg);
	size_t npix = (*OBJ) -> getResolution();
	long dims[] = {2, npix, npix};
	double * out = ypush_d(dims);
	memcpy(out, (*OBJ)->mask(), npix*npix*sizeof(double));
      } else if (yarg_string(iarg)) {
	(*OBJ)->fitsReadMask(ygets_q(iarg));
      } else {
	long ntot;
	long dims[Y_DIMSIZE];
	double const * const in = ygeta_d(iarg, &ntot, dims);
	if (dims[0]==0 && ntot && *in==0) (*OBJ) -> mask(NULL, 0);
	else if (dims[0]==2 && dims[1]==dims[2]) {
	  (*OBJ)->mask(in, dims[1]);
	} else
	  y_error("MASK must be square image");
      }
    }

    /* MASKWRITE */
    if ((iarg=kiargs[++k])>=0) {
      char * fname=const_cast<char*>("");
      iarg+=*rvset;
      if (!yarg_nil(iarg)) fname=ygets_q(iarg);
      (*OBJ)->fitsWriteMask(ygets_q(iarg));
    }


    YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(distance);
    YGYOTO_WORKER_GETSET_DOUBLE(Dmax);
    YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(inclination);
    YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(PALN);
    YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(argument);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(FreqObs);

    /* PROJECTION */
    if ((iarg=kiargs[++k])>=0) { // Set Projection
      iarg+=*rvset;
      long ntot;
      double *proj=ygeta_d(iarg, &ntot, NULL);
      switch (ntot) {
      case 4:
	(*OBJ)->setProjection(proj[0], proj[1], proj[2], proj[3]);
	break;
      case 3:
	(*OBJ)->setProjection(proj[0], proj[1], proj[2]);
	break;
      }
    }

    YGYOTO_WORKER_GETSET4(ObserverPos);
    YGYOTO_WORKER_GETSET4(FourVel);
    YGYOTO_WORKER_GETSET4(Screen1);
    YGYOTO_WORKER_GETSET4(Screen2);
    YGYOTO_WORKER_GETSET4(Screen3);
    YGYOTO_WORKER_GETSET_OBJECT2(spectrometer,Spectrometer);

    ///// METHODS //////
    /* SKYCOORD METHOD */
    if ((iarg=kiargs[++k])>=0) { // skycoord
      if ((*rvset)++) y_error(rmsg);
      long ntot=1;
      double *pos=ygeta_d(iarg, &ntot, NULL);
      if (ntot<4) y_error("POS argument should have at lest 4 elements");
	
      long dims[] = {1, 3};
      double * skypos=ypush_d(dims);
	
      (*OBJ)->coordToXYZ(pos, skypos);
    }
      
    /* RAYCOORD METHOD */
    if ((iarg=kiargs[++k])>=0) { // raycoord
      if ((*rvset)++) y_error(rmsg);
	
      long ntot=1;
      double *pos=ygeta_d(iarg, &ntot, NULL);
      if (ntot<2) y_error("X_Y argument should have at lest 4 elements");
	
      long dims[] = {1,8};
      yarg_drop(1);
      double * coord=ypush_d(dims);
	
      (*OBJ)->getRayCoord(pos[0], pos[1], coord);
    }
      
    YGYOTO_WORKER_XMLWRITE;
    YGYOTO_WORKER_CLONE(Screen);

  }
}

extern "C" {

  void Y_gyoto_Screen(int argc) {
    YGYOTO_CONSTRUCTOR_INIT2(Screen, Screen, Screen, screen);
    gyoto_Screen_eval(OBJ, argc);
  }


}
