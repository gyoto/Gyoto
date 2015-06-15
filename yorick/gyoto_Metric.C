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
#include "ygyoto_private.h"
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


// Needed by the YGYOTO_WORKER_* macros

YGYOTO_YUSEROBJ(Metric, Metric::Generic)
YGYOTO_BASE_CONSTRUCTOR1(Metric,metric)

static char ygyoto_Metric_names[YGYOTO_TYPE_LEN][YGYOTO_MAX_REGISTERED]
={{0}};
static ygyoto_Metric_eval_worker_t *ygyoto_Metric_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Metric_count=0;

extern "C" {
  void gyoto_Metric_eval(void *obj, int argc) {
    SmartPointer<Metric::Generic> *OBJ_ = &((gyoto_Metric*)obj)->smptr;

    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long( (long) (*OBJ_)() );
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    const string kind = (*OBJ_)->kind();

    while (n<ygyoto_Metric_count && kind.compare(ygyoto_Metric_names[n])) ++n;

    if (n<ygyoto_Metric_count && ygyoto_Metric_evals[n]) {
      (*ygyoto_Metric_evals[n])(OBJ_, argc);
      return;
    }

    // Fall-back to generic worker
    static char const * knames[]={
      "unit",
      YGYOTO_METRIC_GENERIC_KW, 0
    };

    YGYOTO_WORKER_INIT(Metric, Generic, knames, YGYOTO_METRIC_GENERIC_KW_N+1);

    YGYOTO_WORKER_SET_UNIT;

    YGYOTO_WORKER_CALL_GENERIC(Metric);

  }
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

void ygyoto_Metric_generic_eval(SmartPointer<Metric::Generic>*OBJ,
				int *kiargs, int *piargs,
				int *rvset, int *paUsed, char * unit) {
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
    ypush_double((*OBJ)->SysPrimeToTdot(pos, vel));
  }

  if ((iarg=kiargs[++k])>=0) { // nullifycoord
    if (debug()) cerr << "In nullifycoord" << endl;
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    double * pos=ygeta_d(iarg,&ntot,0);
    if (ntot!=4) y_error("POS must have 4 elements");
    double * vel=ygeta_d(piargs[0],&ntot,0);
    if (ntot!=3) y_error("VEL must have 3 elements");
    dims[0]=1; dims[1]=8;
    double *coord = ypush_d(dims);
    for (int i=0; i<4; ++i) coord[i]=pos[i];
    for (int i=0; i<3; ++i) coord[i+5]=vel[i];
    (*OBJ)->nullifyCoord(coord);
  }

  // kind
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*OBJ)->kind().c_str());
  }

  // coordkind
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("COORDKIND is readonly");
    ypush_long((*OBJ)->coordKind());
  }

  YGYOTO_WORKER_SETPARAMETER;

  // ScalarProd
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    ntot=0;
    double * pos = ygeta_d(iarg, &ntot, dims);
    if (!dims[0] || dims[1]<4)
      y_error("scalarprod: pos must be at least 4 elements long");
    double * u1 = ygeta_d(piargs[0], &ntot, dims);
    if (!dims[0] || dims[1]<4)
      y_error("scalarprod: u1 must be at least 4 elements long");
    double * u2 = ygeta_d(piargs[1], &ntot, dims);
    if (!dims[0] || dims[1]<4)
      y_error("scalarprod: u2 must be at least 4 elements long");
    ypush_double((*OBJ)->ScalarProd(pos, u1, u2));
  }

  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(mass);

  YGYOTO_WORKER_GETSET_DOUBLE2(deltaMin);
  YGYOTO_WORKER_GETSET_DOUBLE2(deltaMax);
  YGYOTO_WORKER_GETSET_LONG2(keplerian);

  // Unit length
  if ((iarg=kiargs[++k])>=0) { // unitLength()
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("UNITLENGTH is readonly");
    ypush_double((*OBJ)->unitLength(unit?unit:""));
  }

  // circularvelocity
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    double * coords = ygeta_d(iarg, &ntot, dims);
    if (!dims[0] || dims[1]<4)
      y_error("syntax: circularvelocity=array(double, 4 ...)");
    double dir =1;
    if (piargs[0] >= 0) dir = ygets_d(piargs[0]) >= 0 ? 1 : -1;
    long d1 = dims[1];
    long npoints=ntot/d1;
    dims[1]=4;
    double * vels = ypush_d(dims);
    for (long i=0; i<npoints; ++i, coords+=d1, vels+=4)
      (*OBJ)->circularVelocity(coords, vels, dir);
  }

  // christoffel
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    double * coords = ygeta_d(iarg, &ntot, dims);
    if (!dims[0] || dims[1]<4)
      y_error("syntax: christoffel=array(double, 4)");

    Idx i_idx (piargs[0], 4);
    if (i_idx.isNuller()) return;
    Idx j_idx (piargs[1], 4);
    if (j_idx.isNuller()) return;
    Idx a_idx (piargs[2], 4);
    if (a_idx.isNuller()) return;
    long ni=i_idx.getNElements();
    long nj=j_idx.getNElements();
    long na=a_idx.getNElements();

    dims[0]=i_idx.getNDims()+j_idx.getNDims()+a_idx.getNDims();
    size_t offset=0;
    if (i_idx.getNDims()) dims[++offset]=ni;
    if (j_idx.getNDims()) dims[++offset]=nj;
    if (a_idx.getNDims()) dims[++offset]=na;
    double * data=ypush_d(dims);
    double dst[4][4][4];

    if (dims[0]==3 && dims[1]==4 && dims[2]==4 && dims[3]==4) {
      (*OBJ)->christoffel(dst, coords);
      memcpy(data, &dst[0][0], 4*4*4*sizeof(double));
    } else {
      size_t i, j, a;
      for ( a=a_idx.first() ; a_idx.valid() ; a=a_idx.next() )
	for ( i=i_idx.first() ; i_idx.valid() ; i=i_idx.next() )
	  for ( j=j_idx.first() ; j_idx.valid() ; j=j_idx.next() )
	    *(data++) = (*OBJ)->christoffel(coords, int(a-1), int(i-1), int(j-1));
    }

  } 

  YGYOTO_WORKER_XMLWRITE;
  YGYOTO_WORKER_CLONE(Metric);
  YGYOTO_WORKER_HELP;

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

  dims[0]=i_idx.getNDims()+j_idx.getNDims();
  size_t offset=0;
  if (i_idx.getNDims()) dims[++offset]=ni;
  if (j_idx.getNDims()) dims[++offset]=nj;
  double * data=ypush_d(dims);
  double dst[4][4];

  if (dims[0]==2 && dims[1]==4 && dims[2]==4) {
    (*OBJ)->gmunu(dst, x);
    memcpy(data, &dst[0][0], 4*4*sizeof(double));
    return;
  }

  size_t i, j;
  for ( j=j_idx.first() ; j_idx.valid() ; j=j_idx.next() )
    for ( i=i_idx.first() ; i_idx.valid() ; i=i_idx.next() )
      *(data++) = (*OBJ)->gmunu(x, int(i-1), int(j-1));

}
