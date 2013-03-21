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
#include <yapi.h>
#include <pstdlib.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
using namespace std;


using namespace Gyoto;

extern "C" {
  typedef struct gyoto_Screen { SmartPointer<Screen> screen; } gyoto_Screen;
  void gyoto_Screen_free(void *obj);
  void gyoto_Screen_print(void *obj);
  void gyoto_Screen_eval(void *obj, int n);
  static y_userobj_t gyoto_Screen_obj =
    {const_cast<char*>("gyoto_Screen"), &gyoto_Screen_free, &gyoto_Screen_print, &gyoto_Screen_eval, 0, 0};

  // SCREEN CLASS
  void gyoto_Screen_free(void *obj) {
    if (((gyoto_Screen*)obj)->screen) {
      ((gyoto_Screen*)obj)->screen=NULL;
    } else printf("null pointer\n");
  }

  void gyoto_Screen_print(void *obj) {
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    try { rest = Factory(((gyoto_Screen*)obj)->screen).format(); }
    YGYOTO_STD_CATCH;
    size_t pos=0, len;
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO screen ",0);
    //    SmartPointer<Screen> gg = ((gyoto_Screen*)obj)->screen;
    //    y_print(gg->getKind().c_str(),0);
#endif
  }

  void gyoto_Screen_eval(void *obj, int argc) {
    int rvset[1]={0}, paUsed[1]={0};
    int k=-1;
    char const * rmsg="Cannot set return value more than once";
    char const * pmsg="Cannot use positional argument more than once";

    SmartPointer<Screen> *OBJ = NULL; //&((gyoto_Screen*)obj)->screen;

    if (!obj) {
      obj = ypush_obj(&gyoto_Screen_obj, sizeof(gyoto_Screen));
      OBJ = &(((gyoto_Screen*)obj)->screen);
      *OBJ = new Screen();
    } else if (argc==1 && yarg_nil(0)) { // If no parameters, return pointer
      ypush_long( (long) ((gyoto_Screen*)obj)->screen() );
      return;
    } else {
      OBJ = &(((gyoto_Screen*)obj)->screen);
      *ypush_Screen()=*OBJ;
    }
      
    static char const * knames[]={
      "unit",
      "metric",
      "time","fov","resolution",
      "distance", "dmax", "inclination", "paln", "argument",
      "freqobs",
      "projection", "observerpos",
      "fourvel", "screen1", "screen2", "screen3",
      "spectro",
      "skycoord",  "raycoord",
      "xmlwrite", "clone",
      0
    };
#define nkw 22
    static long kglobs[nkw+1];
    int kiargs[nkw];
    int piargs[]={-1,-1,-1,-1};
    // push default return value: need to drop before pushing another one
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    char * unit=NULL;
      
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Screen takes at most 4 positional arguments");
      }
    }

    YGYOTO_WORKER_SET_UNIT;
    YGYOTO_WORKER_GETSET_OBJECT(Metric);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(Time);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(FieldOfView);
    YGYOTO_WORKER_GETSET_LONG(Resolution);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(Distance);
    YGYOTO_WORKER_GETSET_DOUBLE(Dmax);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(Inclination);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(PALN);
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(Argument);
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
    YGYOTO_WORKER_GETSET_OBJECT(Spectrometer);

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

// PUBLIC API

SmartPointer<Screen> *yget_Screen(int iarg) {
  return &((gyoto_Screen*)yget_obj(iarg, &gyoto_Screen_obj))->screen;
}
SmartPointer<Screen> *ypush_Screen() {
  gyoto_Screen* obj = (gyoto_Screen*)ypush_obj(&gyoto_Screen_obj, sizeof(gyoto_Screen));
  return &(obj->screen);
}

int yarg_Screen(int iarg) {
  return yget_obj(iarg,0)==gyoto_Screen_obj.type_name;
}


// YAPI FUNCTIONS

extern "C" {

  void Y_gyoto_Screen(int argc) {
    void *obj = NULL;
    if (yarg_Screen(argc)) obj = yget_obj(--argc, &gyoto_Screen_obj) ;
    gyoto_Screen_eval(obj, argc);
  }


}
