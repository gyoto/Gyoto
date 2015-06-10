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
#include "ygyoto_idx.h"
#include "yapi.h"

#include <cstring>
#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace YGyoto;
using namespace Gyoto::Metric;

#define OBJ gg

void ygyoto_KerrKS_eval(SmartPointer<Metric::Generic> *gg_, int argc) {

  static char const * knames[]={
    "unit", "spin", "horizonsecurity", "gmunu_up", "jacobian",
    YGYOTO_METRIC_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Metric, KerrKS, knames, YGYOTO_METRIC_GENERIC_KW_N+5);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE2(spin);
  YGYOTO_WORKER_GETSET_DOUBLE2(horizonSecurity);

  // GMUNU_UP
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    GYOTO_DEBUG << "getting gmunu_up\n";
    long  ntot, dims[Y_DIMSIZE];
    double  *x     = ygeta_d(iarg, &ntot, NULL);
    if (ntot<4) y_error("X must have at least four elements");

    Idx i_idx (piargs[0], 4);
    if (i_idx.isNuller()) return;
    Idx j_idx (piargs[1], 4);
    if (j_idx.isNuller()) return;
    long ni=i_idx.getNElements();
    long nj=j_idx.getNElements();

    dims[0]=i_idx.getNDims()+j_idx.getNDims();
    size_t offset=0;
    if (i_idx.getNDims()) dims[++offset]=ni;
    if (j_idx.getNDims()) dims[++offset]=nj;
    double * data=ypush_d(dims);
    double dst[4][4];

    (*OBJ)->gmunu_up(dst, x);

    if (dims[0]==2 && dims[1]==4 && dims[2]==4) {
      memcpy(data, &dst[0][0], 4*4*sizeof(double));
    } else {
      size_t i, j;
      for ( j=j_idx.first() ; j_idx.valid() ; j=j_idx.next() )
	for ( i=i_idx.first() ; i_idx.valid() ; i=i_idx.next() )
	  *(data++) = dst[i-1][j-1];
    }
  }

  // JACOBIAN
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    GYOTO_DEBUG << "getting jacobian\n";
    long  ntot, dims[Y_DIMSIZE];
    double  *x     = ygeta_d(iarg, &ntot, NULL);
    if (ntot<4) y_error("X must have at least four elements");

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

    (*OBJ)->jacobian(dst, x);

    if (dims[0]==3 && dims[1]==4 && dims[2]==4 && dims[3]==4) {
      memcpy(data, &dst[0][0][0], 4*4*4*sizeof(double));
    } else {
      size_t i, j, a;
      for ( a=a_idx.first() ; a_idx.valid() ; a=a_idx.next() )
	for ( j=j_idx.first() ; j_idx.valid() ; j=j_idx.next() )
	  for ( i=i_idx.first() ; i_idx.valid() ; i=i_idx.next() )
	    *(data++) = dst[a-1][i-1][j-1];
    }
  }

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
