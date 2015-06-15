/*
    Copyright 2011-2015 Thibaut Paumard

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
#include "ygyoto_private.h"
#include "yapi.h"
#include "pstdlib.h"

#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif

#include <cstring>
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

static char ygyoto_Astrobj_names[YGYOTO_TYPE_LEN][YGYOTO_MAX_REGISTERED]
={{0}};
static ygyoto_Astrobj_eval_worker_t *ygyoto_Astrobj_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Astrobj_count=0;

YGYOTO_YUSEROBJ(Astrobj, Astrobj::Generic)
YGYOTO_BASE_CONSTRUCTOR1(Astrobj,astrobj)

extern "C" {
  void gyoto_Astrobj_eval(void *obj, int argc) {
    GYOTO_DEBUG << endl;
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Astrobj*)obj)->smptr());
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    SmartPointer<Astrobj::Generic> * OBJ = &(((gyoto_Astrobj*)obj)->smptr);
    const string kind = (*OBJ)->kind();

    while (n<ygyoto_Astrobj_count && kind.compare(ygyoto_Astrobj_names[n])) ++n;

    if (n<ygyoto_Astrobj_count && ygyoto_Astrobj_evals[n]) {
      (*ygyoto_Astrobj_evals[n])(OBJ, argc);
      return;
    }

    // Possibly call higher-level base class worker

    // push default return value: need to drop before pushing another one
    *ypush_Astrobj()=*OBJ;
    int rvset[1]={0}, paUsed[1]={0};
    int iarg=argc, parg=0;
    int piargs[]={-1,-1,-1,-1};

    enum BASE_NAME { GENERIC, THINDISK };
    BASE_NAME base = GENERIC;
    if (dynamic_cast<ThinDisk const * const>((*OBJ)())) base = THINDISK;

    static char const * knames_thindisk[]={"unit", YGYOTO_THINDISK_GENERIC_KW, 0};
    static char const * knames_generic[]={"unit", YGYOTO_ASTROBJ_GENERIC_KW, 0};
    static long kglobs[YGYOTO_ASTROBJ_BASE_MAX_KW_N+2];
    static int kiargs[YGYOTO_ASTROBJ_BASE_MAX_KW_N+1];

    char ** knames=NULL;
    ygyoto_Astrobj_generic_eval_t * worker;
    // ThinDisk
    switch (base) {
    case GENERIC:
      knames=const_cast<char**>(knames_generic);
      worker = &ygyoto_Astrobj_generic_eval;
      break;
    case THINDISK:
      knames = const_cast<char**>(knames_thindisk);
      worker = &ygyoto_ThinDisk_generic_eval;
      break;
    default:
      y_error("BUG: unknown base type");
    }

    yarg_kw_init(knames, kglobs, kiargs);
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Astrobj takes at most 4 positional arguments");
      }
    }

    char * unit=NULL;
    int k =-1;

    /* UNIT */
    YGYOTO_WORKER_SET_UNIT;

    (*worker)(OBJ, kiargs+k+1, piargs, rvset, paUsed, unit);

  }

}


void ygyoto_Astrobj_register(char const*const name, ygyoto_Astrobj_eval_worker_t* on_eval){
  int n;
  if (ygyoto_Astrobj_count==YGYOTO_MAX_REGISTERED)
    y_error("Too many Astrobjs registered");
  for (n=0; n<ygyoto_Astrobj_count; ++n)
    if (!strcmp(ygyoto_Astrobj_names[n], name)) 
      return;

  strcpy(ygyoto_Astrobj_names[ygyoto_Astrobj_count], name);
  ygyoto_Astrobj_evals[ygyoto_Astrobj_count++]=on_eval;
}

void ygyoto_Astrobj_generic_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*OBJ,
				int *kiargs, int *piargs,
				 int *rvset, int *paUsed, char * unit) {
  int k=-1, iarg;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  // avoid ununsed warning
  if (rmsg && pmsg) {}

  if (debug())
    for (int i=0; i<YGYOTO_ASTROBJ_GENERIC_KW_N; ++i)
      cerr << "DEBUG: Astrobj_generic_eval: kiargs[" << i << "]="
	   << kiargs[i] << endl;

  /* METRIC */
  YGYOTO_WORKER_GETSET_OBJECT2(metric,Metric);
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(rMax);
  YGYOTO_WORKER_GETSET_LONG2(opticallyThin);
  YGYOTO_WORKER_XMLWRITE;

  // kind
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    if (debug()) cerr << "kiargs=" << kiargs << endl;
    if ((*rvset)++) y_error(rmsg);
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*OBJ)->kind().c_str());
  }

  /* SETPARAMETER */
  YGYOTO_WORKER_SETPARAMETER;
  YGYOTO_WORKER_CLONE(Astrobj);
  YGYOTO_WORKER_HELP;

  if (debug()) cerr << "DEBUG: out of Astrobj_generic_eval"<< endl;

}
