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

#include <iostream>
#include <cstring>
#include "ygyoto_idx.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

using namespace std;
using namespace Gyoto;
using namespace YGyoto;

typedef struct ySceneryThreadWorkerArg {
#ifdef HAVE_PTHREAD
  pthread_mutex_t * mutex;
  pthread_t * parent;
#endif
  Idx *i_idx, *j_idx;
  Scenery *sc;
  Photon * ph;
  Astrobj::Properties *data;
  double * impactcoords;
  size_t res;
} ySceneryThreadWorkerArg ;

static void * ySceneryThreadWorker (void *arg) {
  /*
    This is the real ray-tracing loop. It may be called by multiple
    threads in parallel, launched from ::rayTrace
   */

  ySceneryThreadWorkerArg *larg = static_cast<ySceneryThreadWorkerArg*>(arg);

  // Each thread needs its own Photon, clone cached Photon
  // it is assumed to be already initialized with spectrometer et al.
  Photon * ph = larg -> ph;
#ifdef HAVE_PTHREAD
  if (larg->mutex) {
    GYOTO_DEBUG << "locking mutex\n";
    pthread_mutex_lock(larg->mutex);
    GYOTO_DEBUG << "mutex locked\n";
    ph = larg -> ph -> clone();
    GYOTO_DEBUG << "unlocking mutex\n";
    pthread_mutex_unlock(larg->mutex);
    GYOTO_DEBUG << "mutex unlocked\n";
  }
#endif

  // local variables to store our parameters
  size_t i, j;
  Astrobj::Properties data;
  double * impactcoords = NULL;

  size_t count=0;

  while (1) {
    /////// 1- get input and output parameters and update them for next access
    //// i and j are input, data and impactcoords are where to store
    //// output.  we must get them and increase them so that another
    //// thread can get the next values while we integrate.
#ifdef HAVE_PTHREAD
    // lock mutex so we can safely read and update i, j et al.
    if (larg->mutex) pthread_mutex_lock(larg->mutex);
#endif
    // copy i & j 
    i = larg->i_idx->current(); j = larg->j_idx->current();

    // check whether there remains something to compute
    if (!larg->j_idx->valid()
	|| (larg->j_idx->isLast() && !larg->i_idx->valid())) {
      // terminate, but first...
#ifdef HAVE_PTHREAD
      // ...unlock mutex so our siblings can access i & j and terminate too
      if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif
      break;
    }

    // print some info
    if (larg->i_idx->isFirst() &&
	verbose() >= GYOTO_QUIET_VERBOSITY && !impactcoords) {
      cout << "\rRay-tracing scenery: j = " << j << flush ;
    }
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "i = " << i << ", j = " << j << endl;
#   endif

    // update i & j
    larg->i_idx->next();
    if (!larg->i_idx->valid()) {
      larg->j_idx->next(); larg->i_idx->first();
    }

    // copy output pointers and update them
    data = *larg->data; ++(*larg->data);
    if (larg->impactcoords) impactcoords=larg->impactcoords+((j-1)*larg->res+i-1)*16;

#ifdef HAVE_PTHREAD
    // unlock mutex so our siblings can can access i, j et al. and procede
    if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif

    ////// 2- do the actual work.
    (*larg->sc)(i, j, &data, impactcoords, ph);
    ++count;
  }
#ifdef HAVE_PTHREAD
  if (larg->mutex) {
    delete ph;
    pthread_mutex_lock(larg->mutex);
  }
  GYOTO_MSG << "\nThread terminating after integrating " << count << " photons";
  if (larg->mutex) pthread_mutex_unlock(larg->mutex);
# endif
  return NULL;
}

extern "C" {

  void* _ypush_Scenery();
  void* _yget_Scenery(int iarg);
  int _yarg_Scenery(int iarg);

  void Y_gyoto_Scenery(int);

  // SCENERY CLASS
  // Opaque Yorick object
  typedef struct gyoto_Scenery {
    SmartPointer<Scenery> scenery;
    //  char type[YGYOTO_TYPE_LEN];
  } gyoto_Scenery;
  void gyoto_Scenery_free(void *obj) {
    if (((gyoto_Scenery*)obj)->scenery) {
      ((gyoto_Scenery*)obj)->scenery = NULL;
    } else printf("Freeing unattached Scenery object\n");
  }
  void gyoto_Scenery_print(void *obj) {
    GYOTO_DEBUG << endl;
#ifdef GYOTO_USE_XERCES
    if (debug()) {
      cerr << "DEBUG: Printing Gyoto Scenery"<<endl;
      cerr << "DEBUG:          Pointer: ";
      cerr << ((gyoto_Scenery*)obj)->scenery()<<endl;
      cerr << "DEBUG:          Creating Factory" << endl;
    }
    string rest = Factory(((gyoto_Scenery*)obj)->scenery).format(), sub="";
    if (debug())
      cerr << "DEBUG:          Printing" << endl;
    size_t pos=0, len;
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO Scenery object",0);
    //    y_print(((gyoto_Scenery*)obj)->type,0);
#endif
  }
  void gyoto_Scenery_eval(void *obj, int n) {
    int rvset[1]={0}, paUsed[1]={0}, builder=0;
    gyoto_Scenery *s_obj = (gyoto_Scenery*)obj;
    double * impactcoords = NULL;
    bool precompute = 0;

    if (!s_obj) {
      builder=1;
      s_obj = (gyoto_Scenery*)(obj=_ypush_Scenery());
    } else {
      *ypush_Scenery()=s_obj->scenery;
    }

    // Parse keywords
    static char const *knames[] = {
      "get_pointer",
      "unit",
      "metric", "screen", "astrobj", "delta", "tmin", "quantities",
      "xmlwrite", "clone",
      "impactcoords", "nthreads",
      0
    };
    static long kglobs[13];
    int kiargs[12], piargs[]={-1, -1, -1};
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);

    int iarg=n, parg=0;

    while (iarg>0) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<3) piargs[parg++]=iarg--;
	else
	  y_error("gyoto_Scenery::eval takes at most 3 positional arguments");
      }
    }

    if (builder) {
      if (yarg_string(piargs[0])) {
	if (debug())
	  cerr << "DEBUG: gyoto_Scenery() creating new Scenery from file\n";
#ifdef GYOTO_USE_XERCES
	Factory *factory = new Factory(ygets_q(piargs[0]));
	string kind="";
	kind = factory->getKind();
	if (kind.compare("Scenery"))
	  y_error("Only Scenery supported when reading XML");
	s_obj->scenery = factory -> getScenery();
	delete factory;
#else
	y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
      } else {
	if (debug())
	  cerr << "DEBUG: gyoto_Scenery() creating new empty Scenery\n";
	s_obj->scenery=new Scenery();
      }
    }
    SmartPointer<Scenery> sc=s_obj->scenery;
    int k=-1;

    // Get one member
    if (yarg_true(kiargs[++k])) { // get pointer
      if ((*rvset)++) y_error("Only one return value possible");
      yarg_drop(1);
      ypush_long((long)(sc()));
    }

    char * unit=NULL;

    ///////// ACCESSORS //////////
    /* UNIT */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      GYOTO_DEBUG << "get unit" << endl;
      unit = ygets_q(iarg);
    }

    /* METRIC */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // metric=: Getting
      if ((*rvset)++) y_error("Only one return value possible");
      (*ypush_Metric()) = sc->getMetric();
      } else sc->setMetric(*yget_Metric(iarg));
    }

    /* SCREEN */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // screen=: Getting
	if ((*rvset)++) y_error("Only one return value possible");
	*ypush_Screen()=sc->getScreen();
      } else // Setting
	sc->setScreen(*yget_Screen(iarg));
    }

    /* ASTROBJ */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // astrobj=  :   Getting
	if ((*rvset)++) y_error("Only one return value possible");
	*ypush_Astrobj()=sc->getAstrobj();
      } else                // astrobj=ao:   Setting
	sc->setAstrobj(*yget_Astrobj(iarg));
    }

    /* DELTA */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // delta=      : Getting
	if ((*rvset)++) y_error("Only one return value possible");
	ypush_double(sc->getDelta(unit?unit:""));
      } else                // delta=double: Setting
	sc->setDelta(ygets_d(iarg), unit?unit:"");
    }

    /* TMIN */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // tmin=      : Getting
	if ((*rvset)++) y_error("Only one return value possible");
	ypush_double(sc->getTmin(unit?unit:""));
      } else                // tmin=double: Setting
	sc->setTmin(ygets_d(iarg), unit?unit:"");
    }

    /* QUANTITIES */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // quantities=      : Getting
	if ((*rvset)++) y_error("Only one return value possible");
	Quantity_t quant = sc->getRequestedQuantities();
	size_t nk = sc->getScalarQuantitiesCount();
	long rquant = nk>1?1:0;
	long dims[2] = { rquant, nk };
	ystring_t *squant = ypush_q(dims);
	size_t k = 0;
	char *tk =
	  strtok(const_cast<char*>(sc->getRequestedQuantitiesString().c_str()),
		 " \n\t");
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
	sc -> setRequestedQuantities(quants);
      }
    }

    // Save to file
    if ((iarg=kiargs[++k])>=0) { // xmlwrite
      iarg+=*rvset;
#ifdef GYOTO_USE_XERCES
      char *filename=ygets_q(iarg);
      Factory(sc).write(filename);
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
    }

    /* CLONE */
    if ((iarg=kiargs[++k])>=0) {
      if ((*rvset)++) y_error("Only one return value possible");
      *ypush_Scenery() = sc->clone();
    }

    /* IMPACTCOORDS */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // impactcoords=  :   Getting
	precompute = 1;
      } else {              // impaccoords=double(16,res,res): Setting
	long ntot;
	long dims[Y_DIMSIZE];
	size_t res=sc->getScreen()->getResolution();
	impactcoords = ygeta_d(iarg, &ntot, dims);
	if (dims[0] != 3 || dims[1] != 16 || dims[2] != res || dims[3] != res)
	  y_error("dimsof(impactcoords) != [3,16,res,res]");
      }
    }

    /* NTHREADS */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // nthreads=      : Getting
	if ((*rvset)++) y_error("Only one return value possible");
	ypush_long(sc->getNThreads());
      } else                // nthreads=long  : Setting
         sc->setNThreads(ygets_l(iarg));
    }

    // Get ray-traced image if there is a supplementary positional argument
    if (!builder && // don't ray-trace on construction...
	!*rvset && // has a return value already been set?
	(((n<=3 && piargs[n-1]>=0) || (piargs[1]>=0)) // positional argument?
	 || precompute || impactcoords)
	) { 
      size_t res=sc->getScreen()->getResolution();
      if ((*rvset)++) y_error("Only one return value possible");
      if ((*paUsed)++) y_error("Only one keyword may use positional arguments");
      GYOTO_DEBUG << "rank: " << yarg_rank(piargs[0]) << endl;

      SmartPointer<Spectrometer> spr = sc->getScreen()->getSpectrometer();
      size_t nbnuobs = spr()? spr->getNSamples() : 0;

      Idx i_idx (piargs[0], res);
      if (i_idx.isNuller()) return;
      Idx j_idx (piargs[1], res);
      if (j_idx.isNuller()) return;
      long ni=i_idx.getNElements();
      long nj=j_idx.getNElements();
      long nelem=ni*nj;

      long nk = 0; int rquant = 0; ystring_t * squant = NULL;

      if (piargs[2] < 0 || yarg_nil(piargs[2])) {
	Quantity_t quant = sc->getRequestedQuantities();
	nk = sc->getScalarQuantitiesCount();
	if (quant & GYOTO_QUANTITY_SPECTRUM) nk += 1;
	if (quant & GYOTO_QUANTITY_BINSPECTRUM) nk += 1;
	rquant = nk>1?1:0;
	long dims[2] = { 1, nk };
	squant = ypush_q(dims);
	size_t k = 0;
	char *tk =
	  strtok(const_cast<char*>(sc->getRequestedQuantitiesString().c_str()),
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

      long ndims=i_idx.getNDims()+j_idx.getNDims()
	+ (precompute ? 1 : (((rquant>=1)||has_sp||has_bsp)));
      GYOTO_DEBUG << "i_idx.getNDims()=" << i_idx.getNDims()
		  << ", j_idx.getNDims()" << j_idx.getNDims()
		  << ", precompute ? 1 : (((rquant>=1)||has_sp||has_bsp))"
		  << (precompute ? 1 : (((rquant>=1)||has_sp||has_bsp))) <<endl;
      long dims[4]={ndims};
      size_t offset=0;
      if (precompute)       dims[++offset]=16;
      if (i_idx.getNDims()) dims[++offset]=ni;
      if (j_idx.getNDims()) dims[++offset]=nj;
      if ( !precompute && ((rquant>=1)||has_sp||has_bsp) ) dims[++offset]=nk;
      GYOTO_DEBUG << "precompute=" << precompute << ", nk=" << nk
		  << ", nbnuobs="<<nbnuobs << ", data=ypush_d({"<< dims[0]
		  << ", " << dims[1] << ", " << dims[2] << ", " << dims[3]
		  << "})\n";
      double * data=ypush_d(dims);

      Astrobj::Properties prop;
      SmartPointer<Screen> screen = sc -> getScreen();
#     ifdef HAVE_UDUNITS
      if (data) sc->setPropertyConverters(&prop);
      screen->mapPixUnit();
      if (intu != "") prop.setIntensityConverter(intu);
      if (spu  != "") prop.setSpectrumConverter(spu);
      if (bspu != "") prop.setBinSpectrumConverter(bspu);
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
      if (i_idx.isDouble() ||j_idx.isDouble()) {
	prop.init(nbnuobs);
	Photon ph(sc->getMetric(), sc->getAstrobj(), screen,
		  i_idx.getDVal(), j_idx.getDVal());
	ph.hit(&prop);
      } else {
	screen -> computeBaseVectors();
	double coord[8];
	screen -> getRayCoord(size_t(i_idx.first()),
			      size_t(j_idx.first()),
			      coord);
	Photon ph(sc->getMetric(), sc->getAstrobj(), screen, coord);
	ph.setSpectrometer(screen->getSpectrometer());

	ySceneryThreadWorkerArg larg;
	larg.sc=sc;
	larg.ph=&ph;
	larg.data=&prop;
	larg.impactcoords=impactcoords;
	larg.i_idx=&i_idx;
	larg.j_idx=&j_idx;
	larg.res=res;

#       ifdef HAVE_PTHREAD
	larg.mutex  = NULL;
	pthread_mutex_t mumu = PTHREAD_MUTEX_INITIALIZER;
	pthread_t * threads = NULL;
	pthread_t pself = pthread_self();
	larg.parent = &pself;
	size_t nthreads = sc -> getNThreads();
	if (nthreads >= 2) {
	  threads = new pthread_t[nthreads-1];
	  larg.mutex  = &mumu;
	  for (size_t th=0; th < nthreads-1; ++th) {
	    if (pthread_create(threads+th, NULL,
			       ySceneryThreadWorker,
			       static_cast<void*>(&larg)) < 0)
	      y_error("Error creating thread");
	  }
	}
#       endif

	// Call worker on the parent thread
	(*ySceneryThreadWorker)(static_cast<void*>(&larg));

#       ifdef HAVE_PTHREAD
	// Wait for the child threads
	if (nthreads>=2)
	  for (size_t th=0; th < nthreads-1; ++th)
	    pthread_join(threads[th], NULL);
#       endif
	if (impactcoords==NULL) GYOTO_MSG << endl;
      }
    }
  }

  static y_userobj_t gyoto_Scenery_obj = {
    const_cast<char*>("gyoto_Scenery"),
    &gyoto_Scenery_free, &gyoto_Scenery_print, &gyoto_Scenery_eval,
    0, 0
  };

  int _yarg_Scenery(int iarg) {
    return yget_obj(iarg, 0)==gyoto_Scenery_obj.type_name;
  }

  void * _ypush_Scenery() {
    return ypush_obj(&gyoto_Scenery_obj, sizeof(gyoto_Scenery));
  }

  void * _yget_Scenery(int iarg) {
    return yget_obj(iarg, &gyoto_Scenery_obj);
  }

  // Constructor
  void Y_gyoto_Scenery(int n) {
    void* obj = NULL;
    if (_yarg_Scenery(n-1)) obj = _yget_Scenery(--n);
    gyoto_Scenery_eval(obj, n);
  }

  void Y_gyoto_Scenery_rayTrace(int argc) {
    size_t imin=1, imax=-1, jmin=1, jmax=-1;
    if (argc<1) y_error("gyoto_Scenery_rayTrace takes at least 1 argument");
    gyoto_Scenery * s_obj=(gyoto_Scenery*)yget_obj(argc-1, &gyoto_Scenery_obj);
    Scenery * scenery = s_obj->scenery;
    if (argc>=2 && !yarg_nil(argc-2)) imin=ygets_l(argc-2);
    if (argc>=3 && !yarg_nil(argc-3)) imax=ygets_l(argc-3);
    if (argc>=4 && !yarg_nil(argc-4)) jmin=ygets_l(argc-4);
    if (argc>=5 && !yarg_nil(argc-5)) jmax=ygets_l(argc-5);

    size_t res;
    try {res=scenery->getScreen()->getResolution();}
    YGYOTO_STD_CATCH;

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

    size_t nbnuobs=0, nbdata;
    Quantity_t quantities;
    try {
      quantities = scenery -> getRequestedQuantities();
      if (quantities & (GYOTO_QUANTITY_SPECTRUM | GYOTO_QUANTITY_BINSPECTRUM)) {
	SmartPointer<Spectrometer> spr=scenery->getScreen()->getSpectrometer();
	if (!spr) throwError("Spectral quantity requested but "
			     "no spectrometer specified!");
	nbnuobs = spr -> getNSamples();
      }
      nbdata= scenery->getScalarQuantitiesCount();
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

    try {scenery -> rayTrace(imin, imax, jmin, jmax, &data,
			     ipctout?NULL:impactcoords);}
    YGYOTO_STD_CATCH;

  }

}
// PUBLIC API

SmartPointer<Scenery> *yget_Scenery(int iarg) {
  return &((gyoto_Scenery*)yget_obj(iarg, &gyoto_Scenery_obj))->scenery;
}
SmartPointer<Scenery> *ypush_Scenery() {
  gyoto_Scenery* obj = (gyoto_Scenery*)ypush_obj(&gyoto_Scenery_obj, sizeof(gyoto_Scenery));
  return &(obj->scenery);
}

int yarg_Scenery(int iarg) {
  return yget_obj(iarg,0)==gyoto_Scenery_obj.type_name;
}


