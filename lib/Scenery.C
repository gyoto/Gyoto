/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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

#include "GyotoUtils.h"
#include "GyotoScenery.h"
#include "GyotoPhoton.h"
#include "GyotoFactoryMessenger.h"

#include <cmath>
#include <cfloat>
#include <cstring>
#include <cstdlib>

#define DEFAULT_TMIN -DBL_MAX

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include <sys/time.h>    /* for benchmarking */


using namespace Gyoto;
using namespace std;

/*Scenery::Scenery() :
  gg_(NULL), screen_(NULL), obj_(NULL),
  deltatau_(0.01) {}
*/
Scenery::Scenery() :
  gg_(NULL), screen_(NULL), obj_(NULL), delta_(0.01),
  quantities_(0), ph_(), tmin_(DEFAULT_TMIN), nthreads_(0){}

Scenery::Scenery(const Scenery& o) :
  SmartPointee(o),
  gg_(NULL), screen_(NULL), obj_(NULL), delta_(o.delta_),
  quantities_(o.quantities_), ph_(o.ph_), tmin_(o.tmin_), nthreads_(o.nthreads_)
{
  // We have up to 3 _distinct_ clones of the same Metric.
  // Keep only one.
  if (o.gg_()) gg_=o.gg_->clone();
  if (o.screen_()) {
    screen_=o.screen_->clone();
    screen_->setMetric(gg_);
  }
  if (o.obj_()) {
    obj_=o.obj_->clone();
    obj_->setMetric(gg_);
  }
}
Scenery * Scenery::clone() const { return new Scenery(*this); }

/*Scenery::Scenery(SmartPointer<Metric::Generic> met, SmartPointer<Screen> screen, SmartPointer<Astrobj::Generic> obj) :
  gg_(met), screen_(screen), obj_(obj),
  deltatau_(0.01)
{}
*/
Scenery::Scenery(SmartPointer<Metric::Generic> met, SmartPointer<Screen> screen, SmartPointer<Astrobj::Generic> obj) :
  gg_(met), screen_(screen), obj_(obj), delta_(0.01),
  quantities_(0)
{
  if (screen_) screen_->setMetric(gg_);
  if (obj_) obj_->setMetric(gg_);
}

Scenery::~Scenery() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "freeing metric\n";
# endif
  gg_ = NULL;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "freeing screen\n";
# endif
  screen_ = NULL;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "freeing astrobj\n";
# endif
  obj_ = NULL;
 }

SmartPointer<Metric::Generic> Scenery::getMetric() { return gg_; }

void Scenery::setMetric(SmartPointer<Metric::Generic> met) {
  gg_ = met;
  if (!screen_) screen_ = new Screen ();
  screen_ -> setMetric(gg_);
  if (obj_) obj_ -> setMetric(gg_);
}

SmartPointer<Screen> Scenery::getScreen() { return screen_; }

void Scenery::setScreen(SmartPointer<Screen> screen) {
  screen_ = screen;
  if (gg_) screen_ -> setMetric (gg_) ;
}

SmartPointer<Astrobj::Generic> Scenery::getAstrobj() { return obj_; }
void Scenery::setAstrobj(SmartPointer<Astrobj::Generic> obj) {
  obj_ = obj;
  if (gg_) obj_ -> setMetric (gg_) ;
}

double Scenery::getDelta() const { return delta_; }
void Scenery::setDelta(double d) { delta_ = d; }

void  Scenery::setNThreads(size_t n) { nthreads_ = n; }
size_t Scenery::getNThreads() const { return nthreads_; }

typedef struct SceneryThreadWorkerArg {
#ifdef HAVE_PTHREAD
  pthread_mutex_t * mutex;
  pthread_t * parent;
#endif
  size_t i, j, imin, imax, jmin, jmax;
  Scenery *sc;
  Photon * ph;
  Astrobj::Properties *data;
  double * impactcoords;
} SceneryThreadWorkerArg ;

static void * SceneryThreadWorker (void *arg) {
  /*
    This is the real ray-tracing loop. It may be called by multiple
    threads in parallel, launched from ::rayTrace
   */

  SceneryThreadWorkerArg *larg = static_cast<SceneryThreadWorkerArg*>(arg);

  // Each thread needs its own Photon, clone cached Photon
  // it is assumed to be already initialized with spectrometer et al.
  Photon * ph = larg -> ph;
#ifdef HAVE_PTHREAD
  if (larg->mutex) ph = larg -> ph -> clone();
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
    i = larg->i; j = larg->j;
    if (j > larg->jmax || (j==larg->jmax && i>larg->imax)) {
      // terminate, but first...
#ifdef HAVE_PTHREAD
      // ...unlock mutex so our siblings can access i & j and terminate too
      if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif
      break;
    }
    // update i & j
    ++larg->i;
    if (larg->i > larg->imax) {
      ++larg->j; larg->i=larg->imin;
    }

    // copy output pointers and update them
    data = *larg->data; ++(*larg->data);
    if (larg->impactcoords) {
      impactcoords = larg->impactcoords; larg->impactcoords+=16;
    }

#ifdef HAVE_PTHREAD
    // unlock mutex so our siblings can can access i, j et al. and procede
    if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif

    ////// 2- do the actual work.
    if (i==larg->imin && verbose() >= GYOTO_QUIET_VERBOSITY && !impactcoords)
      cout << "\rj = " << j << " / " << larg->jmax << " " << flush;
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "i = " << i << ", j = " << j << endl;
#   endif
    (*larg->sc)(i, j, &data, impactcoords, ph);
    ++count;
  }
#ifdef HAVE_PTHREAD
  if (larg->mutex) delete ph;
#endif
  GYOTO_MSG << "\nThread terminating after integrating " << count << " photons";
  return NULL;
}

void Scenery::rayTrace(size_t imin, size_t imax,
		       size_t jmin, size_t jmax,
		       Astrobj::Properties *data,
		       double * impactcoords) {

  /*
     Ray-trace now is multi-threaded. What it does is
       - some initialization
       - launch nthreads_ - 1  thread working on of SceneryThreadWorker
       - call SceneryThreadWorker itself rather than sleeping
       - wait for the other threads to be terminated
       - some housekeeping
   */

  const size_t npix = screen_->getResolution();
  imax=(imax<=(npix)?imax:(npix));
  jmax=(jmax<=(npix)?jmax:(npix));
  screen_->computeBaseVectors();
         // Necessary for KS integration, computes relation between
         // observer's x,y,z coord and KS X,Y,Z coord. Will be used to
         // compute photon's initial tangent vector.
     // Note : this is a BUG if this is required, should be done automagically.

  /// initialize photon once. It will be cloned.
  SmartPointer<Spectrometer> spr = screen_->getSpectrometer();
  ph_.setSpectrometer(spr);
  ph_.setTmin(tmin_);
  double coord[8];
  screen_ -> getRayCoord(imin,jmin, coord);
  ph_ . setInitialCondition(gg_, obj_, coord);
  // delta is reset in operator()



  SceneryThreadWorkerArg larg;
  larg.sc=this;
  larg.ph=&ph_;
  larg.data=data;
  larg.impactcoords=impactcoords;
  larg.i=imin;
  larg.j=jmin;
  larg.imin=imin;
  larg.imax=imax;
  larg.jmin=jmin;
  larg.jmax=jmax;

  struct timeval tim;
  double start, end;
  gettimeofday(&tim, NULL);  
  start=tim.tv_sec+(tim.tv_usec/1000000.0);  

#ifdef HAVE_PTHREAD
  larg.mutex  = NULL;
  pthread_mutex_t mumu = PTHREAD_MUTEX_INITIALIZER;
  pthread_t * threads = NULL;
  pthread_t pself = pthread_self();
  larg.parent = &pself;
  if (nthreads_ >= 2) {
    threads = new pthread_t[nthreads_-1];
    larg.mutex  = &mumu;
    for (size_t th=0; th < nthreads_-1; ++th) {
      if (pthread_create(threads+th, NULL,
			 SceneryThreadWorker, static_cast<void*>(&larg)) < 0)
	throwError("Error creating thread");
    }
  }
#endif

  // Call worker on the parent thread
  (*SceneryThreadWorker)(static_cast<void*>(&larg));


#ifdef HAVE_PTHREAD
  // Wait for the child threads
  if (nthreads_>=2)
    for (size_t th=0; th < nthreads_-1; ++th)
      pthread_join(threads[th], NULL);
#endif

  gettimeofday(&tim, NULL);  
  end=tim.tv_sec+(tim.tv_usec/1000000.0);  

  GYOTO_MSG << "\nRaytraced "<< (jmax-jmin+1) * (imax-imin+1)
	    << " photons in " << end-start
	    << "s using " << nthreads_ << " thread(s)\n";

}

void Scenery::operator() (
			  size_t i, size_t j,
			  Astrobj::Properties *data, double * impactcoords,
			  Photon *ph
			  ) {
  double coord[8];
  SmartPointer<Spectrometer> spr = screen_->getSpectrometer();
  size_t nbnuobs = spr() ? spr -> getNSamples() : 0;
  SmartPointer<Metric::Generic> gg = NULL;
  SmartPointer<Astrobj::Generic> obj = NULL;
  if (!ph) {
    // if Photon was passed, assume it was initiliazed already. Don't
    // touch its metric and astrobj. Else, update cached photon. Photon
    // is passed in particular when called in a multi-threaded
    // environment: it may really need to work on a given copy of the object.
    ph = &ph_;
    ph -> setSpectrometer(spr);
    ph -> setTmin(tmin_);
    obj=obj_;
    gg=gg_;
  }
  // Always reset delta
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "reset delta" << endl;
# endif
  ph -> setDelta(delta_);
  ph -> setTmin(tmin_);

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "init nbnuobs" << endl;
# endif
  if (data) data -> init(nbnuobs); // Initialize requested quantities to 0. or DBL_MAX

  if (impactcoords) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "impactcoords set" << endl;
#   endif
    if(impactcoords[0] != DBL_MAX) {
      ph -> setInitialCondition(gg, obj, impactcoords+8);
      ph -> resetTransmission();
      obj_ -> processHitQuantities(ph,impactcoords+8,impactcoords,0.,data);
    }
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "impactcoords not set" << endl;
#   endif
    screen_ -> getRayCoord(i,j, coord);
    ph -> setInitialCondition(gg, obj, coord);
    ph -> hit(data);
  }
}

void Scenery::setRequestedQuantities(Gyoto::Quantity_t quant)
{quantities_=quant;}
void Scenery::setRequestedQuantities(std::string squant) {
  quantities_=0;
  char * tk = strtok(const_cast<char*>(squant.c_str()), " \t\n");
  while (tk != NULL) {
    if (!strcmp(tk, "Intensity"))
      quantities_ |= GYOTO_QUANTITY_INTENSITY;
    else if (!strcmp(tk, "EmissionTime"))
      quantities_ |= GYOTO_QUANTITY_EMISSIONTIME;
    else if (!strcmp(tk, "MinDistance"))
      quantities_ |= GYOTO_QUANTITY_MIN_DISTANCE;
    else if (!strcmp(tk, "FirstDistMin"))
      quantities_ |= GYOTO_QUANTITY_FIRST_DMIN;
    else if (!strcmp(tk, "Redshift"))
      quantities_ |= GYOTO_QUANTITY_REDSHIFT;
    else if (!strcmp(tk, "ImpactCoords"))
      quantities_ |= GYOTO_QUANTITY_IMPACTCOORDS;
    else if (!strcmp(tk, "Spectrum"))
      quantities_ |= GYOTO_QUANTITY_SPECTRUM;
    else if (!strcmp(tk, "BinSpectrum"))
      quantities_ |= GYOTO_QUANTITY_BINSPECTRUM;
    else if (!strcmp(tk, "User1"))
      quantities_ |= GYOTO_QUANTITY_USER1;
    else if (!strcmp(tk, "User2"))
      quantities_ |= GYOTO_QUANTITY_USER2;
    else if (!strcmp(tk, "User3"))
      quantities_ |= GYOTO_QUANTITY_USER3;
    else if (!strcmp(tk, "User4"))
      quantities_ |= GYOTO_QUANTITY_USER4;
    else if (!strcmp(tk, "User5"))
      quantities_ |= GYOTO_QUANTITY_USER5;
    else throwError("ScenerySubcontractor(): unkwon quantity"); 
    tk = strtok(NULL, " \t\n");
  }
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "("<<squant<<"): " << "quantities_=" << quantities_ << endl;
# endif
}
Gyoto::Quantity_t Scenery::getRequestedQuantities() const {
  return quantities_?quantities_:(obj_()?obj_->getDefaultQuantities():0);
}

std::string Scenery::getRequestedQuantitiesString() const {
  string squant = "";
  Quantity_t quantities
    = quantities_?quantities_:(obj_()?obj_->getDefaultQuantities():0);
  if (quantities & GYOTO_QUANTITY_INTENSITY   ) squant+="Intensity ";
  if (quantities & GYOTO_QUANTITY_EMISSIONTIME) squant+="EmissionTime ";
  if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) squant+="MinDistance ";
  if (quantities & GYOTO_QUANTITY_FIRST_DMIN  ) squant+="FirstDistMin ";
  if (quantities & GYOTO_QUANTITY_REDSHIFT    ) squant+="Redshift ";
  if (quantities & GYOTO_QUANTITY_IMPACTCOORDS) squant+="ImpactCoords ";
  if (quantities & GYOTO_QUANTITY_SPECTRUM    ) squant+="Spectrum ";
  if (quantities & GYOTO_QUANTITY_BINSPECTRUM ) squant+="BinSpectrum ";
  if (quantities & GYOTO_QUANTITY_USER1       ) squant+="User1 ";
  if (quantities & GYOTO_QUANTITY_USER2       ) squant+="User2 ";
  if (quantities & GYOTO_QUANTITY_USER3       ) squant+="User3 ";
  if (quantities & GYOTO_QUANTITY_USER4       ) squant+="User4 ";
  if (quantities & GYOTO_QUANTITY_USER5       ) squant+="User5 ";
  return squant;
}

size_t Scenery::getScalarQuantitiesCount() const {
  size_t nquant=0;
  Quantity_t quantities
    = quantities_?quantities_:(obj_()?obj_->getDefaultQuantities():0);
  if (quantities & GYOTO_QUANTITY_INTENSITY   ) ++nquant;
  if (quantities & GYOTO_QUANTITY_EMISSIONTIME) ++nquant;
  if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) ++nquant;
  if (quantities & GYOTO_QUANTITY_FIRST_DMIN  ) ++nquant;
  if (quantities & GYOTO_QUANTITY_REDSHIFT    ) ++nquant;
  //  SPECTRUM is not a SCALAR, don't add the following:
  //  if (quantities & GYOTO_QUANTITY_SPECTRUM    ) ++nquant;
  //  if (quantities & GYOTO_QUANTITY_BINSPECTRUM ) ++nquant;
  //  Idem IMPACTCOORDS:
  //  if (quantities & GYOTO_QUANTITY_IMPACTCOORDS) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER1       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER2       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER3       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER4       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER5       ) ++nquant;
  return nquant;
}

double Scenery::getTmin() const { return tmin_; }
void Scenery::setTmin(double tmin) { tmin_ = tmin; }

#ifdef GYOTO_USE_XERCES
void Scenery::fillElement(FactoryMessenger *fmp) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "fmp -> setMetric (gg_) ;" << endl;
# endif
  if (gg_)     fmp -> setMetric (gg_) ;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG <<"fmp -> setScreen (screen_) ;" << endl;
# endif
  if (screen_) fmp -> setScreen (screen_) ;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG <<"fmp -> setAstrobj (obj_) ;" << endl;
# endif
  if (obj_)    fmp -> setAstrobj (obj_) ;

  if (delta_ != GYOTO_DEFAULT_DELTA) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG <<"fmp -> setParameter (\"Delta\", "<<delta_<<") ;" << endl;
#   endif
    fmp -> setParameter ("Delta", delta_);
  }

  if (getRequestedQuantities()) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG <<"fmp -> setParameter (\"Quantities\", \""
		<<getRequestedQuantitiesString()<<"\") ;" << endl;
#   endif
    fmp -> setParameter("Quantities", getRequestedQuantitiesString());
  }

  if (tmin_ != DEFAULT_TMIN) fmp -> setParameter("MinimumTime", tmin_);
  if (nthreads_) fmp -> setParameter("NThreads", nthreads_);
}

SmartPointer<Scenery> Gyoto::ScenerySubcontractor(FactoryMessenger* fmp) {

  string name="", content="";
  double delta = GYOTO_DEFAULT_DELTA ;
  SmartPointer<Metric::Generic> gg = NULL;
  SmartPointer<Screen> scr = NULL;
  SmartPointer<Astrobj::Generic> ao = NULL;
  string squant = "";
  double tmin = DEFAULT_TMIN;
  size_t nthreads = 0;

  gg = fmp->getMetric();
  scr= fmp->getScreen();
  ao = fmp->getAstrobj();


  while (fmp->getNextParameter(&name, &content)) {
    char* tc = const_cast<char*>(content.c_str());
    if (name=="Delta") delta = atof(tc);
    if (name=="Quantities") squant = content;
    if (name=="MinimumTime") tmin = atof(tc);
    if (name=="NThreads") nthreads = atoi(tc);
  }

  SmartPointer<Scenery> sc = new Scenery(gg, scr, ao);
  sc -> setDelta(delta);
  sc -> setTmin(tmin);
  sc -> setNThreads(nthreads);
  if (squant!="") sc -> setRequestedQuantities(squant);
  return sc;
}
#endif
