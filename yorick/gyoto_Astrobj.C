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

extern "C" {
  // ASTROBJ CLASS
  // Opaque Yorick object
  typedef struct gyoto_Astrobj {
    SmartPointer<Astrobj::Generic> astrobj;
    //char type[YGYOTO_TYPE_LEN];
  } gyoto_Astrobj;
  void gyoto_Astrobj_free(void *obj) {
    if (((gyoto_Astrobj*)obj)->astrobj) {
      ((gyoto_Astrobj*)obj)->astrobj = NULL;
    } else printf("null pointer\n");
  }
  void gyoto_Astrobj_print(void *obj) {
#ifdef GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_EXPR(obj);
#endif
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    size_t pos=0, len=0;
    try {rest = Factory(((gyoto_Astrobj*)obj)->astrobj).format();}
    YGYOTO_STD_CATCH;
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO Astrobj object of type ",0);
    y_print(((gyoto_Astrobj*)obj)->astrobj->getKind().c_str(),0);
#endif
  }
  void gyoto_Astrobj_eval(void *obj, int argc) {
    GYOTO_DEBUG << endl;
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Astrobj*)obj)->astrobj());
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    SmartPointer<Astrobj::Generic> * ao = &(((gyoto_Astrobj*)obj)->astrobj);
    const string kind = (*ao)->getKind();

    while (n<ygyoto_Astrobj_count && kind.compare(ygyoto_Astrobj_names[n])) ++n;

    if (n<ygyoto_Astrobj_count && ygyoto_Astrobj_evals[n]) {
      (*ygyoto_Astrobj_evals[n])(ao, argc);
      return;
    }

    // Possibly call hiegher-level base class worker

    // push default return value: need to drop before pushing another one
    *ypush_Astrobj()=*ao;
    int rvset[1]={0}, paUsed[1]={0};
    int iarg=argc, parg=0;
    int piargs[]={-1,-1,-1,-1};

    enum BASE_NAME { GENERIC, THINDISK };
    BASE_NAME base = GENERIC;
    if (dynamic_cast<ThinDisk const * const>((*ao)())) base = THINDISK;

    static char const * knames_thindisk[]={YGYOTO_THINDISK_GENERIC_KW, 0};
    static char const * knames_generic[]={YGYOTO_ASTROBJ_GENERIC_KW, 0};
    static long kglobs[YGYOTO_ASTROBJ_BASE_MAX_KW_N+1];
    static int kiargs[YGYOTO_ASTROBJ_BASE_MAX_KW_N];

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
      y_error("BUG: unkown base type");
    }

    yarg_kw_init(knames, kglobs, kiargs);
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Astrobj takes at most 4 positional arguments");
      }
    }

    (*worker)(ao, kiargs, piargs, rvset, paUsed);

  }
  static y_userobj_t gyoto_Astrobj_obj =
    {const_cast<char*>("gyoto_Astrobj"), &gyoto_Astrobj_free,
     &gyoto_Astrobj_print, &gyoto_Astrobj_eval, 0, 0};

}

SmartPointer<Astrobj::Generic>* yget_Astrobj(int iarg) {
  return &(((gyoto_Astrobj*)yget_obj(iarg, &gyoto_Astrobj_obj))->astrobj);
}

SmartPointer<Astrobj::Generic>* ypush_Astrobj() {
  return &(((gyoto_Astrobj*)ypush_obj(&gyoto_Astrobj_obj, sizeof(gyoto_Astrobj)))->astrobj);
}

int yarg_Astrobj(int iarg) {
  return yget_obj(iarg,0)==gyoto_Astrobj_obj.type_name;
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

void ygyoto_Astrobj_generic_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*ao,
				int *kiargs, int *piargs,
				int *rvset, int *paUsed) {
  int k=-1, iarg;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  if (debug())
    for (int i=0; i<YGYOTO_ASTROBJ_GENERIC_KW_N; ++i)
      cerr << "DEBUG: Astrobj_generic_eval: kiargs[" << i << "]="
	   << kiargs[i] << endl;

  /* METRIC */
  if ((iarg=*(kiargs++))>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "in";
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Metric() = (*ao)->getMetric();
    } else
      (*ao)->setMetric(*yget_Metric(iarg)) ;
    if (debug()) cerr << "out";
  }

  /* RMAX */
  if ((iarg=*(kiargs++))>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error("rmsg");
      ypush_double((*ao)->getRmax());
    } else
      (*ao)->setRmax(ygets_d(iarg));
  }

  /* FLAG_RADTRANSF */
  if ((iarg=*(kiargs++))>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error("rmsg");
      ypush_long((*ao)->getFlag_radtransf());
    } else
      (*ao)->setFlag_radtransf(ygets_l(iarg));
  }

  // Save to file
  if ((iarg=*(kiargs++))>=0) { // xmlwrite
    iarg+=*rvset;
#ifdef GYOTO_USE_XERCES
    char *filename=ygets_q(iarg);
    Factory(*ao).write(filename);
#else
    y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
  }
    
    // kind
  if ((iarg=*(kiargs++))>=0) {
    iarg+=*rvset;
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    if (debug()) cerr << "kiargs=" << kiargs << endl;
    if ((*rvset)++) y_error(rmsg);
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*ao)->getKind().c_str());
  }

  /* SETPARAMETER */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if ((*rvset)++) y_error("rmsg");
    if ((*paUsed)++) y_error("pmsg");
    string name = ygets_q(iarg);
    string content = ygets_q(*piargs);
    (*ao)->setParameter(name, content, "");
  }

  /* CLONE */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    *ypush_Astrobj() = (*ao)->clone();
  }

  if (debug()) cerr << "DEBUG: out of Astrobj_generic_eval"<< endl;

}

extern "C" {
  void Y_gyoto_Astrobj(int argc) 
  {
    int rvset[1]={0}, paUsed[1]={0};
    SmartPointer<Astrobj::Generic> *ao = NULL;
    //    char *obj_type=(char*)yget_obj(argc-1,0);
    //    if (obj_type && //!strcmp(obj_type, "gyoto_Astrobj")) {
    //    if (yget_obj(argc-1,0) && yarg_typeid(argc-1)==Y_OPAQUE) {
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      *ypush_Astrobj() = *ao; // push back astrobj
    } else { // Constructor mode
      ao = ypush_Astrobj();
      *rvset=1;
    }

    static char const * knames[]={
      YGYOTO_ASTROBJ_GENERIC_KW,
      0
    };
    static long kglobs[YGYOTO_ASTROBJ_GENERIC_KW_N+12];
    int kiargs[YGYOTO_ASTROBJ_GENERIC_KW_N+11];
    int piargs[]={-1,-1,-1,-1};
  
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
  
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Astrobj takes at most 4 positional arguments");
      }
    }

    // if rvset==1, constructor mode:
    if (rvset[1]) {
      if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
	try { *ao = Factory(ygets_q(piargs[0])).getAstrobj(); } 
	YGYOTO_STD_CATCH;
	paUsed[1]=1;
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
      } else y_error("Cannot allocate object of virtual class Astrobj");
    }

    ygyoto_Astrobj_generic_eval(ao, kiargs, piargs, rvset, paUsed);

  }

  void
  Y_is_gyoto_Astrobj(int argc)
  {
    ypush_long(yarg_Astrobj(0));
  }

}
