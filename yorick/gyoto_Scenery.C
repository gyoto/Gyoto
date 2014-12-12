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

#include <float.h>
#include <GyotoScenery.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif
#include <GyotoAstrobj.h>
#include <GyotoUtils.h>
#include "yapi.h"
#include "pstdlib.h"
#include "ygyoto.h"
#include "ygyoto_private.h"

#include <iostream>
#include <cstring>
#include "ygyoto_idx.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

using namespace std;
using namespace Gyoto;
using namespace YGyoto;

YGYOTO_YUSEROBJ(Scenery, Scenery)

extern "C" {

  void Y_gyoto_Scenery(int);

  void gyoto_Scenery_eval(void *obj, int argc) {
    SmartPointer<Scenery> *OBJ_ = &(((gyoto_Scenery*)obj)->smptr);

    { // First look for introspected Properties at the beginning of
      // the argument list. Only four syntaxes possible for each
      // property:
      //   obj, Property=value;
      //   obj, Property=value, unit="unit";
      //   obj(Property=);
      //   obj(Property=, unit="unit");

      long kidx;
      Property const * prop;
      std::string pname="", unit="";
      bool unit_found, prop_pushed=false;
      ypush_nil(); ++argc;
      while ( (  argc > 0                 ) &&  
	      ( (kidx=yarg_key(argc-1)) >=0 ) &&
	      ( (prop=(*OBJ_)->property(pname=yfind_name(kidx)))) ) {
	if ( (kidx=yarg_key(argc-3)) >=0 &&
	     !strcmp(yfind_name(kidx),"unit") ) {
	  unit=ygets_q(argc-4);
	  unit_found=true;
	} else {
	  unit_found=false;
	  unit="";
	}
	if (yarg_nil(argc-2)) {
	  if (prop_pushed++) y_error("Can push only one return value");
	  yarg_drop(1);
	  ypush_property(*OBJ_, *prop, -1, pname, unit);
	}
	else yget_property(*OBJ_, *prop, argc-2, pname, unit);
	argc -= unit_found?4:2;
      }
      if (prop_pushed) return;
      yarg_drop(1); --argc;
    }

    double * impactcoords = NULL;
    bool precompute = 0;

    static char const *knames[] = {
      "get_pointer",
      "unit",
      "metric", "screen", "astrobj", "delta", "tmin", "quantities", "adaptive",
      "maxiter", "integrator", "deltamin", "deltamax", "deltamaxoverr",
      "abstol", "reltol", 
      "xmlwrite", "clone", "clonephoton",
      "impactcoords", "nthreads", "nprocesses",
      "mpispawn", "mpiclone",
      0
    };

    YGYOTO_WORKER_INIT1(Scenery, Scenery, knames, 24)

    // Get pointer
    if (yarg_true(kiargs[++k])) {
      if ((*rvset)++) y_error("Only one return value possible");
      ypush_long((long)((*OBJ)()));
    }

    YGYOTO_WORKER_SET_UNIT;
    YGYOTO_WORKER_GETSET_OBJECT2(metric,Metric);
    YGYOTO_WORKER_GETSET_OBJECT2(screen,Screen);
    YGYOTO_WORKER_GETSET_OBJECT2(astrobj,Astrobj);
    YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(delta);
    YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(tMin);

    /* QUANTITIES */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // quantities=      : Getting
	if ((*rvset)++) y_error("Only one return value possible");
	Quantity_t quant = (*OBJ)->getRequestedQuantities();
	long nk = long((*OBJ)->getScalarQuantitiesCount());
	long rquant = nk>1?1:0;
	long dims[2] = { rquant, nk };
	ystring_t *squant = ypush_q(dims);
	long k = 0;
	string requested=(*OBJ)->requestedQuantitiesString();
	char *tk = strtok(const_cast<char*>(requested.c_str()), " \n\t");
	while (tk!=NULL) {
	  if (k>=nk) y_error("BUG: too many tokens in quantity list");
	  squant[k++] = p_strcpy(tk);
	  tk = strtok(NULL, " \n\t");
	}
      }else {               // quantities=["Q1", "Q2"...]: Setting
	long k, nk;
	ystring_t * squant = ygeta_q(iarg, &nk, NULL);
	string quants = squant[0];
	for (k=1; k<nk; ++k) {
	  quants += " ";
	  quants += squant[k];
	}
	(*OBJ) -> requestedQuantitiesString(quants);
      }
    }

    YGYOTO_WORKER_GETSET_LONG2( adaptive );
    YGYOTO_WORKER_GETSET_LONG2( maxiter );
    YGYOTO_WORKER_GETSET_STRING2( integrator );
    YGYOTO_WORKER_GETSET_DOUBLE2( deltaMin );
    YGYOTO_WORKER_GETSET_DOUBLE2( deltaMax );
    YGYOTO_WORKER_GETSET_DOUBLE2( deltaMaxOverR );
    YGYOTO_WORKER_GETSET_DOUBLE2( absTol );
    YGYOTO_WORKER_GETSET_DOUBLE2( relTol );
    YGYOTO_WORKER_XMLWRITE;
    YGYOTO_WORKER_CLONE(Scenery);

    /* CLONEPHOTON */
    if ((iarg=kiargs[++k])>=0) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Photon() = (*OBJ)->clonePhoton();
    }

    /* IMPACTCOORDS */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // impactcoords=  :   Getting
	precompute = 1;
      } else {              // impaccoords=double(16,res,res): Setting
	long ntot;
	long dims[Y_DIMSIZE];
	size_t res=(*OBJ)->screen()->resolution();
	impactcoords = ygeta_d(iarg, &ntot, dims);
	if (dims[0] != 3 || dims[1] != 16 || dims[2] != res || dims[3] != res)
	  y_error("dimsof(impactcoords) != [3,16,res,res]");
      }
    }

    YGYOTO_WORKER_GETSET_LONG2(nThreads);
    YGYOTO_WORKER_GETSET_LONG2(nProcesses);
#ifdef HAVE_MPI
    YGYOTO_WORKER_RUN( mpiSpawn(ygets_l(iarg)) );
    YGYOTO_WORKER_RUN( mpiClone() );
#else
    if ((iarg=kiargs[++k])>=0) GYOTO_WARNING << "No MPI in this GYOTO" << endl;
    if ((iarg=kiargs[++k])>=0) GYOTO_WARNING << "No MPI in this GYOTO" << endl;
#endif


    // Get ray-traced image if there is a supplementary positional argument
    if (
	!*rvset && // has a return value already been set?
	(((argc>0 && argc<=3 && piargs[argc-1]>=0) || (piargs[1]>=0)) // positional argument?
	 || precompute || impactcoords)
	) { 
      size_t res=(*OBJ)->screen()->resolution();
      if ((*rvset)++) y_error("Only one return value possible");
      if ((*paUsed)++) y_error("Only one keyword may use positional arguments");
      GYOTO_DEBUG << "rank: " << yarg_rank(piargs[0]) << endl;

      SmartPointer<Spectrometer::Generic> spr = (*OBJ)->screen()->spectrometer();
      size_t nbnuobs = spr()? spr->nSamples() : 0;

      Idx i_idx (piargs[0], res);
      if (i_idx.isNuller()) return;
      Idx j_idx (piargs[1], res);
      if (j_idx.isNuller()) return;

      bool is_double=false;

      long ndims, dims[Y_DIMSIZE]={0};
      long nelem;

      if (i_idx.isDouble() || j_idx.isDouble()) {
	if (!i_idx.isDouble() || !j_idx.isDouble())
	  throwError("i and j must be of same type (double or long)");
	is_double=true;
	nelem=i_idx.getNElements();
	long const *idims=i_idx.getDims();
	if (j_idx.getNElements() != nelem) {
	  if (nelem==1) {
	    nelem=j_idx.getNElements();
	    idims=j_idx.getDims();
	  }
	  else throwError("alpha and delta must be conformable");
	}
	for (int m=0; m<idims[0]+1; ++m) dims[m]=idims[m];
      } else {
	nelem=i_idx.getNElements()*j_idx.getNElements();
	dims[0]=0;
	if (precompute)       dims[++dims[0]]=16;
	if ( (precompute + i_idx.getDims()[0]+j_idx.getDims()[0]) >= Y_DIMSIZE )
	  throwError("Too many dimensions");
	long const *idims=i_idx.getDims();
	for (int m=1; m<=idims[0];++m) dims[++dims[0]]=idims[m];
	idims=j_idx.getDims();
	for (int m=1; m<=idims[0];++m) dims[++dims[0]]=idims[m];
      }

      long nk = 0; int rquant = 0; ystring_t * squant = NULL;

      if (piargs[2] < 0 || yarg_nil(piargs[2])) {
	Quantity_t quant = (*OBJ)->getRequestedQuantities();
	nk = (*OBJ)->getScalarQuantitiesCount();
	if (quant & GYOTO_QUANTITY_SPECTRUM) nk += 1;
	if (quant & GYOTO_QUANTITY_BINSPECTRUM) nk += 1;
	rquant = nk>1?1:0;
	long dims[2] = { 1, nk };
	squant = ypush_q(dims);
	size_t k = 0;
	char *tk =
	  strtok(const_cast<char*>((*OBJ)->requestedQuantitiesString().c_str()),
		 " \n\t");
	while (tk!=NULL) {
	  if (k>=nk) y_error("BUG: too many tokens in quantity list");
	  squant[k++] = p_strcpy(tk);
	  tk = strtok(NULL, " \n\t");
	}
      } else {
	GYOTO_DEBUG << "quantities provided online"<<endl;
	rquant = yarg_rank(piargs[2]);
	squant = ygeta_q(piargs[2], &nk, NULL);
	GYOTO_DEBUG << "nk="<<nk<<endl;

      }

      size_t k; int has_sp=0, has_bsp=0;
      string tkk="", qunit="", quantity="", intu, spu, bspu;
      size_t first = 0, last = 0;
      if (nbnuobs)
	for (k=0; k<nk; ++k) {
	  GYOTO_DEBUG << "k=" << k<<", nk="<<nk
		      << ", squant[k]="<<squant[k]<<endl;
	  tkk = squant[k];
	  first = tkk.find("[");
	  last = tkk.size() - 1;
	  if (first < last) {
	    qunit = tkk.substr(first+1, last-first-1);
	    quantity = tkk.substr(0, first);
	    squant[k][first]=0;
	  } else {
	    qunit="";
	    quantity=tkk;
	  }
	  GYOTO_DEBUG << "quantity=" << quantity << ", qunit=" << qunit << endl;
#	  ifndef HAVE_UDUNITS
	  if (qunit != "")
	    GYOTO_WARNING << "gyoto_Scenery(): unit \""<< qunit
			  << "\" ignored, try recompiling Gyoto --with-udunits"
			  << endl;
#         endif
	  if (quantity=="Spectrum") {
	    has_sp=1;
	    spu=qunit;
	  } else if (quantity=="BinSpectrum") {
	    has_bsp=1;
	    bspu=qunit;
	  } else if (quantity=="Intensity") {
	    intu=qunit;
	  }
	}
      if (has_sp) nk+=nbnuobs-1;
      if (has_bsp) nk+=nbnuobs-1;
      if (!has_sp && !has_bsp) nbnuobs=0;

      if ( !precompute && ((rquant>=1)||has_sp||has_bsp) ) {
	++dims[0];
	if (dims[0] >=Y_DIMSIZE) throwError("Too many dimensions");
	dims[dims[0]]=nk;
      }

      GYOTO_DEBUG << "precompute=" << precompute << ", nk=" << nk
		  << ", nbnuobs="<<nbnuobs << ", data=ypush_d({"<< dims[0]
		  << ", " << dims[1] << ", " << dims[2] << ", " << dims[3]
		  << "})\n";

      GYOTO_DEBUG_ARRAY(dims, Y_DIMSIZE);

      double * data=ypush_d(dims);

      Astrobj::Properties prop;
      prop.alloc=false;
      SmartPointer<Screen> screen = (*OBJ) -> screen();
#     ifdef HAVE_UDUNITS
      if (data) (*OBJ)->setPropertyConverters(&prop);
      screen->mapPixUnit();
      if (intu != "") prop.intensityConverter(intu);
      if (spu  != "") prop.spectrumConverter(spu);
      if (bspu != "") prop.binSpectrumConverter(bspu);
      screen->unmapPixUnit();
#     endif

      size_t i, j;
      if (precompute) prop.impactcoords=data;
      else {
	for ( k=0; k<nk-nbnuobs+has_sp+has_bsp; ++k ) {
	  GYOTO_DEBUG << "new quantity '" << squant[k] <<"'"<<endl;
	  if (!strcmp(squant[k], "Intensity")) {
	    if (prop.intensity) y_error("can retrieve property only once");
	    prop.intensity=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "EmissionTime")) {
	    if (prop.time) y_error("can retrieve property only once");
	    prop.time=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "MinDistance")) {
	    if (prop.distance) y_error("can retrieve property only once");
	    prop.distance=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "FirstDistMin")) {
	    if (prop.first_dmin) y_error("can retrieve property only once");
	    prop.first_dmin=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "Redshift")) {
	    if (prop.redshift) y_error("can retrieve property only once");
	    prop.redshift=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "Spectrum")) {
	    if (prop.spectrum) y_error("can retrieve property only once");
	    prop.spectrum=data;
	    prop.offset=nelem;
	    data+=nelem*nbnuobs;
	  } else if (!strcmp(squant[k], "BinSpectrum")) {
	    if (prop.binspectrum) y_error("can retrieve property only once");
	    prop.binspectrum=data;
	    prop.offset=nelem;
	    data+=nelem*nbnuobs;
	  } else if (!strcmp(squant[k], "User1")) {
	    if (prop.user1) y_error("can retrieve property only once");
	    prop.user1=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "User2")) {
	    if (prop.user2) y_error("can retrieve property only once");
	    prop.user2=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "User3")) {
	    if (prop.user3) y_error("can retrieve property only once");
	    prop.user3=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "User4")) {
	    if (prop.user4) y_error("can retrieve property only once");
	    prop.user4=data;
	    data+=nelem;
	  } else if (!strcmp(squant[k], "User5")) {
	    if (prop.user5) y_error("can retrieve property only once");
	    prop.user5=data;
	    data+=nelem;
	  } else y_errorq("unknown quantity: %s", squant[k]);
	}
      }

      Screen::Coord2dSet *ijspec=NULL;
      Screen::Coord1dSet *ispec=NULL, *jspec=NULL;

      if (is_double) {
	size_t sz = i_idx.getNElements() >= j_idx.getNElements()?
	  i_idx.getNElements():j_idx.getNElements();
	if (i_idx.getNElements()==1)
	  ispec = new Screen::RepeatAngle(i_idx.getDVal(), sz);
	else 
	  ispec = new Screen::Angles(i_idx.getDoubleBuffer(), sz);
	if (j_idx.getNElements()==1)
	  jspec = new Screen::RepeatAngle(j_idx.getDVal(), sz);
	else 
	  jspec = new Screen::Angles(j_idx.getDoubleBuffer(), sz);
	ijspec = new Screen::Bucket(*ispec, *jspec);
      } else {

	if (i_idx.isRangeOrScalar())
	  ispec = new Screen::Range
	    (i_idx.range_min(), i_idx.range_max(), i_idx.range_dlt());
	else
	  ispec = new Screen::Indices
	    (reinterpret_cast<size_t const*const>(i_idx.getBuffer()),
	     i_idx.getNElements());

	if (j_idx.isRangeOrScalar())
	  jspec = new Screen::Range
	    (j_idx.range_min(), j_idx.range_max(), j_idx.range_dlt());
	else
	  jspec = new Screen::Indices
	    (reinterpret_cast<size_t const*const>(j_idx.getBuffer()),
	     j_idx.getNElements());

	ijspec = new Screen::Grid(*ispec, *jspec, "\rj = ");
	if (verbose() >= GYOTO_QUIET_VERBOSITY)
	  cout << "\nj = 1/" << jspec->size() << flush;
      }

      (*OBJ)->rayTrace(*ijspec, &prop, impactcoords);

      delete ispec;
      delete jspec;
      delete ijspec;
    } // if (conditions for ray-tracing)
  } //  void gyoto_Scenery_eval(void *obj, int argc);

  // Constructor
  void Y_gyoto_Scenery(int argc) {
    YGYOTO_CONSTRUCTOR_INIT1(Scenery, Scenery, Scenery);
    gyoto_Scenery_eval(OBJ, argc);
  }

void Y_gyoto_Scenery_rayTrace(int argc) {
    size_t imin=1, imax=-1, jmin=1, jmax=-1;
    if (argc<1) y_error("gyoto_Scenery_rayTrace takes at least 1 argument");
    gyoto_Scenery * s_obj=(gyoto_Scenery*)yget_obj(argc-1, &gyoto_Scenery_obj);
    Scenery * scenery = s_obj->smptr;
    if (argc>=2 && !yarg_nil(argc-2)) imin=ygets_l(argc-2);
    if (argc>=3 && !yarg_nil(argc-3)) imax=ygets_l(argc-3);
    if (argc>=4 && !yarg_nil(argc-4)) jmin=ygets_l(argc-4);
    if (argc>=5 && !yarg_nil(argc-5)) jmax=ygets_l(argc-5);

    long res;
    try {res=long(scenery->screen()->resolution());}
    YGYOTO_STD_CATCH;

    if (imax>res) imax=res;
    if (jmax>res) jmax=res;

    Screen::Range irange(imin, imax, 1);
    Screen::Range jrange(jmin, jmax, 1);
    Screen::Grid grid(irange, jrange, "\r j = ");

    double * impactcoords = NULL;
    int ipctout = 0;
    if (argc>=6) {
      int iarg = argc-6;
      long ref = yget_ref(iarg);
      long dims[Y_DIMSIZE] = {3, 16, res, res};
      if (ref >= 0 && yarg_nil(iarg)) {
	impactcoords = ypush_d(dims);
	yput_global(ref, 0);
	ipctout = 1;
      } else {
	long ntot = 0;
	impactcoords = ygeta_d(iarg, &ntot, dims);
	if (dims[0]!=3 || dims[1]!=16 || dims[2]!=res || dims[3]!=res)
	  y_error("Wrong dims for impactcoords");
      }
    }

    long nbnuobs=0, nbdata;
    Quantity_t quantities;
    try {
      quantities = scenery -> getRequestedQuantities();
      if (quantities & (GYOTO_QUANTITY_SPECTRUM | GYOTO_QUANTITY_BINSPECTRUM)) {
	SmartPointer<Spectrometer::Generic> spr=scenery->screen()->spectrometer();
	if (!spr) throwError("Spectral quantity requested but "
			     "no spectrometer specified!");
	nbnuobs = long(spr -> nSamples());
      }
      nbdata= long(scenery->getScalarQuantitiesCount());
    } YGYOTO_STD_CATCH;

    long dims[4]={(nbdata+nbnuobs) > 1 ? 3 : 2, res, res, nbdata+nbnuobs};
    double * vect=ypush_d(dims);

    Astrobj::Properties data;

    size_t curquant=0;
    size_t offset=res*res;

    if (ipctout) data.impactcoords = impactcoords;

    if (quantities & GYOTO_QUANTITY_INTENSITY)
      data.intensity=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_EMISSIONTIME)
      data.time=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_MIN_DISTANCE)
      data.distance=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_FIRST_DMIN)
      data.first_dmin=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_REDSHIFT)
      data.redshift=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_USER1)
      data.user1=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_USER2)
      data.user2=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_USER3)
      data.user3=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_USER4)
      data.user4=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_USER5)
      data.user5=vect+offset*(curquant++);
    if (quantities & GYOTO_QUANTITY_SPECTRUM) {
      data.spectrum=vect+offset*(curquant++);
      data.offset=offset;
    }
    if (quantities & GYOTO_QUANTITY_BINSPECTRUM) {
      data.binspectrum=vect+offset*(curquant++);
      data.offset=offset;
    }

    data.intensity=vect;

    if (verbose() >= GYOTO_QUIET_VERBOSITY)
      cout << endl << flush;

    try {scenery -> rayTrace(grid, &data, ipctout?NULL:impactcoords);}
    YGYOTO_STD_CATCH;

  }

}


