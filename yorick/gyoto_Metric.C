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
#include "ygyoto_idx.h"
#include "yapi.h"
#include "pstdlib.h"
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
#include <cstring>
using namespace std;
using namespace Gyoto;
using namespace YGyoto;

static char ygyoto_Metric_names[YGYOTO_TYPE_LEN][YGYOTO_MAX_REGISTERED]
={{0}};
static ygyoto_Metric_eval_worker_t *ygyoto_Metric_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Metric_count=0;

extern "C" {
  typedef struct gyoto_Metric { SmartPointer<Metric::Generic> metric; char type[YGYOTO_TYPE_LEN];} gyoto_Metric;
  void gyoto_Metric_free(void *obj);
  void gyoto_Metric_print(void *obj);
  void gyoto_Metric_eval(void *obj, int n);
  static y_userobj_t gyoto_Metric_obj =
    {const_cast<char*>("gyoto_Metric"), &gyoto_Metric_free, &gyoto_Metric_print, &gyoto_Metric_eval, 0, 0};

  // METRIC CLASS
  void gyoto_Metric_free(void *obj) {
    if (((gyoto_Metric*)obj)->metric) {
      ((gyoto_Metric*)obj)->metric=NULL;
    } else printf("null pointer\n");
  }

  void gyoto_Metric_print(void *obj) {
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    try { rest = Factory(((gyoto_Metric*)obj)->metric).format(); }
    YGYOTO_STD_CATCH;
    size_t pos=0, len=0;
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO metric of type ",0);
    SmartPointer<Metric::Generic> gg = ((gyoto_Metric*)obj)->metric;
    y_print(gg->getKind().c_str(),0);
#endif
  }

  void gyoto_Metric_eval(void *obj, int argc) {
    SmartPointer<Metric::Generic> gg = ((gyoto_Metric*)obj)->metric;

    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long( (long) gg() );
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    const string kind = gg->getKind();

    while (n<ygyoto_Metric_count && kind.compare(ygyoto_Metric_names[n])) ++n;

    if (n<ygyoto_Metric_count && ygyoto_Metric_evals[n]) {
      (*ygyoto_Metric_evals[n])(&gg, argc);
      return;
    }

    // Fall-back to generic worker
    static char const * knames[]={
      YGYOTO_METRIC_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_METRIC_GENERIC_KW_N+1];
    int kiargs[YGYOTO_METRIC_GENERIC_KW_N];
    int piargs[]={-1,-1,-1,-1};
    // push back metric by default
    *ypush_Metric()=gg;
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Metric takes at most 4 positional arguments");
      }
    }

    int rvset[1]={0}, paUsed[1]={0};
    ygyoto_Metric_generic_eval(&gg, kiargs, piargs, rvset, paUsed);
  }
}

// PUBLIC API

SmartPointer<Metric::Generic> *yget_Metric(int iarg) {
  return &((gyoto_Metric*)yget_obj(iarg, &gyoto_Metric_obj))->metric;
}

SmartPointer<Metric::Generic> *ypush_Metric() {
  gyoto_Metric* obj = (gyoto_Metric*)ypush_obj(&gyoto_Metric_obj, sizeof(gyoto_Metric));
  return &(obj->metric);
}

int yarg_Metric(int iarg) {
  return yget_obj(iarg,0)==gyoto_Metric_obj.type_name;
}

void ygyoto_Metric_register(char const*const name, ygyoto_Metric_eval_worker_t* on_eval){
  int n;
  if (ygyoto_Metric_count==YGYOTO_MAX_REGISTERED)
    y_error("Too many Metrics registered");
  for (n=0; n<ygyoto_Metric_count; ++n)
    if (!strcmp(ygyoto_Metric_names[n], name)) 
      return;

  strcpy(ygyoto_Metric_names[ygyoto_Metric_count], name);
  ygyoto_Metric_evals[ygyoto_Metric_count++]=on_eval;
  //if ((ygyoto_Metric_count) < YGYOTO_METRIC_MAX_REGISTERED)
  //  strcpy(ygyoto_Metric_names[ygyoto_Metric_count], "");
}

void ygyoto_Metric_generic_eval(Gyoto::SmartPointer<Gyoto::Metric::Generic>*gg,
				int *kiargs, int *piargs,
				int *rvset, int *paUsed) {
  int k=-1, iarg=-1;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  long  ntot, dims[Y_DIMSIZE];

  /* METHODS */  
  if ((iarg=kiargs[++k])>=0) { // get_tdot
    if (debug()) cerr << "In get_tdot" << endl;
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    double * pos=ygeta_d(iarg,&ntot,0);
    if (ntot!=4) y_error("POS must have 4 elements");
    double * vel=ygeta_d(piargs[0],&ntot,0);
    if (ntot!=3) y_error("VEL must have 3 elements");
    ypush_double((*gg)->SysPrimeToTdot(pos, vel));
  }

  if ((iarg=kiargs[++k])>=0) { // nullifycoord
    if (debug()) cerr << "In nullifycoord" << endl;
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    double * pos=ygeta_d(iarg,&ntot,0);
    if (ntot!=4) y_error("POS must have 4 elements");
    double * vel=ygeta_d(piargs[0],&ntot,0);
    if (ntot!=3) y_error("VEL must have 3 elements");
    long dims[]= {1, 8};
    double *coord = ypush_d(dims);
    for (int i=0; i<4; ++i) coord[i]=pos[i];
    for (int i=0; i<3; ++i) coord[i+5]=vel[i];
    (*gg)->nullifyCoord(coord);
  }

  // kind
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*gg)->getKind().c_str());
  }

  // Process SET keywords
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) { // mass= : get mass
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*gg)->getMass());
    } else                // mass=m: set mass        
      (*gg)->setMass(ygets_d(iarg)) ;
  }

  // Unit length
  if ((iarg=kiargs[++k])>=0) { // unitLength()
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("UNITLENGTH is readonly");
    ypush_double((*gg)->unitLength());
  }

  // circularvelocity
  if ((iarg=kiargs[++k])>=0) { // unitLength()
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    long ntot=0;
    long dims[Y_DIMSIZE];
    double * coords = ygeta_d(iarg, &ntot, dims);
    if (!dims[0] || dims[1]<4)
      y_error("syntax: circularvelocity=array(double, 4 ...)");
    long dir =1;
    if (piargs[0] >= 0) dir = ygets_l(piargs[0]) >= 0 ? 1 : -1;
    long d1 = dims[1];
    long npoints=ntot/d1;
    dims[1]=4;
    double * vels = ypush_d(dims);
    for (long i=0; i<npoints; ++i, coords+=d1, vels+=4)
      (*gg)->circularVelocity(coords, vels, dir);
  }

  // Save to file
  if ((iarg=kiargs[++k])>=0) { // xmlwrite
    iarg+=*rvset;
    char *filename=ygets_q(iarg);
#ifdef GYOTO_USE_XERCES
    Factory(*gg).write(filename);
#else
    y_error("This GYOTO was compiled without xerces: no xml i/o");
#endif
  }

  /* CLONE */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    *ypush_Metric() = (*gg)->clone();
  }

  if (*rvset || *paUsed || piargs[0]<0 || !yarg_number(piargs[0])) return;
  //else     /* GET G MU NU */
  
  if (debug()) cerr << "DEBUG: gyoto_Metric_generic_eval: getting gmunu\n";
  double  *x     = ygeta_d(piargs[0], &ntot, NULL);
  if (ntot<4) y_error("X must have at least four elements");

  Idx i_idx (piargs[1], 4);
  if (i_idx.isNuller()) return;
  Idx j_idx (piargs[2], 4);
  if (j_idx.isNuller()) return;
  long ni=i_idx.getNElements();
  long nj=j_idx.getNElements();
  long nelem=ni*nj;

  dims[0]=i_idx.getNDims()+j_idx.getNDims();
  size_t offset=0;
  if (i_idx.getNDims()) dims[++offset]=ni;
  if (j_idx.getNDims()) dims[++offset]=nj;
  double * data=ypush_d(dims);

  size_t i, j;
  for ( j=j_idx.first() ; j_idx.valid() ; j=j_idx.next() )
    for ( i=i_idx.first() ; i_idx.valid() ; i=i_idx.next() )
      *(data++) = (*gg)->gmunu(x, i-1, j-1);

}


// YAPI FUNCTIONS

extern "C" {

  void Y_gyoto_Metric(int argc) {
    int rvset[1]={0}, paUsed[1]={0};
    SmartPointer<Metric::Generic> *gg = NULL;
    int builder=0;
    
    if (yarg_Metric(argc-1)) {
      gg = yget_Metric(--argc);
      // Try calling kind-specific worker
      int n=0;
      const string kind = (*gg)->getKind();
      while (n<ygyoto_Metric_count && kind.compare(ygyoto_Metric_names[n])) ++n;
      if (n<ygyoto_Metric_count && ygyoto_Metric_evals[n]) {
	(*ygyoto_Metric_evals[n])(gg, argc);
	return;
      }
      // push back metric
      *ypush_Metric()=*gg;
    } else { // Constructor mode
      gg = ypush_Metric();
      builder=1;
    }

    static char const * knames[]={
      YGYOTO_METRIC_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_METRIC_GENERIC_KW_N+1];
    int kiargs[YGYOTO_METRIC_GENERIC_KW_N];
    int piargs[]={-1,-1,-1,-1};
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Metric takes at most 4 positional arguments");
      }
    }

    // if builder==1, constructor mode:
    if (builder) {
      if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
	try { *gg = Factory(ygets_q(piargs[0])).getMetric(); } 
	YGYOTO_STD_CATCH;
	paUsed[1]=1;
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
      } else y_error("Cannot allocate object of virtual class Metric");
    }

    ygyoto_Metric_generic_eval(gg, kiargs, piargs, rvset, paUsed);
  }

}
