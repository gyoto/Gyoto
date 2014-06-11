/*
    Copyright 2011, 2013 Thibaut Paumard, Frederic Vincent

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

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include <sys/time.h>    /* for benchmarking */


using namespace Gyoto;
using namespace std;

Scenery::Scenery() :
  screen_(NULL), delta_(GYOTO_DEFAULT_DELTA),
  quantities_(0), ph_(), nthreads_(0){}

Scenery::Scenery(SmartPointer<Metric::Generic> met,
		 SmartPointer<Screen> screen,
		 SmartPointer<Astrobj::Generic> obj) :
  screen_(screen), delta_(GYOTO_DEFAULT_DELTA),
  quantities_(0), ph_(), nthreads_(0)
{
  metric(met);
  if (screen_) screen_->metric(met);
  astrobj(obj);
}

Scenery::Scenery(const Scenery& o) :
  SmartPointee(o),
  screen_(NULL), delta_(o.delta_), 
  quantities_(o.quantities_), ph_(o.ph_), 
  nthreads_(o.nthreads_)
{
  if (o.screen_()) {
    screen_=o.screen_->clone();
    screen_->metric(ph_.metric());
  }
}
Scenery * Scenery::clone() const { return new Scenery(*this); }

Scenery::~Scenery() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "freeing screen\n";
# endif
  screen_ = NULL;
 }

SmartPointer<Metric::Generic> Scenery::metric() const { return ph_.metric(); }

void Scenery::metric(SmartPointer<Metric::Generic> met) {
  ph_.metric(met);
  if (!screen_) screen_ = new Screen ();
  screen_ -> metric(met);
}

SmartPointer<Screen> Scenery::screen() const { return screen_; }

void Scenery::screen(SmartPointer<Screen> screen) {
  screen_ = screen;
  if (metric()) screen_ -> metric (metric()) ;
}

SmartPointer<Astrobj::Generic> Scenery::astrobj() const {return ph_.astrobj();}
void Scenery::astrobj(SmartPointer<Astrobj::Generic> obj) { ph_.astrobj(obj); }

double Scenery::delta() const { return delta_; }
double Scenery::delta(const string &unit) const {
  return Units::FromGeometrical(delta(), unit, metric());
}

void Scenery::delta(double d) { delta_ = d; }
void Scenery::delta(double d, const string &unit) {
  delta(Units::ToGeometrical(d, unit, metric()));
}

void  Scenery::nThreads(size_t n) { nthreads_ = n; }
size_t Scenery::nThreads() const { return nthreads_; }

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
  if (larg->mutex) {
    pthread_mutex_lock(larg->mutex);
    ph = larg -> ph -> clone();
    pthread_mutex_unlock(larg->mutex);
  }
#endif

  // local variables to store our parameters
  size_t i, j;
  size_t eol_offset =
    larg->sc->screen()->resolution() - larg->imax + larg->imin -1;
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

    data = *larg->data; 
    // update i, j and pointer
    ++larg->i;
    ++(*larg->data);
    if (larg->impactcoords) {
      impactcoords = larg->impactcoords; larg->impactcoords+=16;
    }
    if (larg->i > larg->imax) {
      ++larg->j; larg->i=larg->imin;
      (*larg->data) += eol_offset;
      if (larg->impactcoords) larg->impactcoords+=16*eol_offset;
    }

#ifdef HAVE_PTHREAD
    // unlock mutex so our siblings can can access i, j et al. and procede
    if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif

    ////// 2- do the actual work.
    if (i==larg->imin && verbose() >= GYOTO_QUIET_VERBOSITY && !impactcoords) {
#     ifdef HAVE_PTHREAD
      if (larg->mutex) pthread_mutex_lock(larg->mutex);
#     endif
      cout << "\rj = " << j << " / " << larg->jmax << " " << flush;
#     ifdef HAVE_PTHREAD
      if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#     endif
    }
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "i = " << i << ", j = " << j << endl;
#   endif
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

  const size_t npix = screen_->resolution();
  imax=(imax<=(npix)?imax:(npix));
  jmax=(jmax<=(npix)?jmax:(npix));
  screen_->computeBaseVectors();
         // Necessary for KS integration, computes relation between
         // observer's x,y,z coord and KS X,Y,Z coord. Will be used to
         // compute photon's initial tangent vector.
     // Note : this is a BUG if this is required, should be done automagically.

  /// initialize photon once. It will be cloned.
  SmartPointer<Spectrometer::Generic> spr = screen_->spectrometer();
  ph_.spectrometer(spr);
  ph_.freqObs(screen_->freqObs());
  double coord[8];
  screen_ -> getRayCoord(imin,jmin, coord);
  ph_ . setInitCoord(coord, -1);
  // delta is reset in operator()

  if (data) {
    setPropertyConverters(data);
    size_t first_index=(jmin-1)*npix + imin -1;
    (*data) += first_index;
    if (impactcoords) impactcoords += first_index * 16;
  }

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
  start=double(tim.tv_sec)+(double(tim.tv_usec)/1000000.0);  

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
  end=double(tim.tv_sec)+(double(tim.tv_usec)/1000000.0);  

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
  SmartPointer<Spectrometer::Generic> spr = screen_->spectrometer();
  size_t nbnuobs = spr() ? spr -> nSamples() : 0;

  if (data) data -> init(nbnuobs); // Initialize requested quantities to 0. or DBL_MAX
  if (!(*screen_)(i,j)) return; // return if pixel is masked out

  if (!ph) {
    // if Photon was passed, assume it was initiliazed already. Don't
    // touch its metric and astrobj. Else, update cached photon. Photon
    // is passed in particular when called in a multi-threaded
    // environment: it may really need to work on a given copy of the object.
    ph = &ph_;
    ph -> spectrometer(spr);
    ph -> freqObs(screen_->freqObs());
  }
  // Always reset delta
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "reset delta" << endl;
# endif
  ph -> delta(delta_);

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "init nbnuobs" << endl;
# endif

  if (impactcoords) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "impactcoords set" << endl;
#   endif
    if(impactcoords[0] != DBL_MAX) {
      ph -> setInitCoord(impactcoords+8, -1);
      ph -> resetTransmission();
      astrobj() -> processHitQuantities(ph,impactcoords+8,impactcoords,0.,data);
    }
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "impactcoords not set" << endl;
#   endif
    screen_ -> getRayCoord(i,j, coord);
    ph -> setInitCoord(coord, -1);
    ph -> hit(data);
  }
}

SmartPointer<Photon> Scenery::clonePhoton() const {
  return ph_.clone();
}


void Scenery::setRequestedQuantities(Gyoto::Quantity_t quant)
{quantities_=quant;}
void Scenery::setRequestedQuantities(std::string squant) {
  quantities_=0;
  char * tk = strtok(const_cast<char*>(squant.c_str()), " \t\n");
  string tkk="", quant="", unit = "";
  size_t first = 0, last = 0;

# ifdef HAVE_UDUNITS
  if (screen_) screen_ -> mapPixUnit();
# endif

  while (tk != NULL) {
    tkk = tk;
    first = tkk.find("[");
    last = tkk.size() - 1;
    if (first < last) {
      unit = tkk.substr(first+1, last-first-1);
      quant = tkk.substr(0, first);
    } else {
      unit="";
      quant=tkk;
    }
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "quant=" << quant << ", unit=" << unit << "." << endl;
#   endif
    tk=const_cast<char*>(quant.c_str());

    if (!strcmp(tk, "Intensity")) {
      quantities_ |= GYOTO_QUANTITY_INTENSITY;
      intensityConverter(unit);
    } else if (!strcmp(tk, "EmissionTime"))
      quantities_ |= GYOTO_QUANTITY_EMISSIONTIME;
    else if (!strcmp(tk, "MinDistance"))
      quantities_ |= GYOTO_QUANTITY_MIN_DISTANCE;
    else if (!strcmp(tk, "FirstDistMin"))
      quantities_ |= GYOTO_QUANTITY_FIRST_DMIN;
    else if (!strcmp(tk, "Redshift"))
      quantities_ |= GYOTO_QUANTITY_REDSHIFT;
    else if (!strcmp(tk, "ImpactCoords"))
      quantities_ |= GYOTO_QUANTITY_IMPACTCOORDS;
    else if (!strcmp(tk, "Spectrum")) {
      quantities_ |= GYOTO_QUANTITY_SPECTRUM;
      spectrumConverter(unit);
    } else if (!strcmp(tk, "BinSpectrum")) {
      quantities_ |= GYOTO_QUANTITY_BINSPECTRUM;
      binSpectrumConverter(unit);
    } else if (!strcmp(tk, "User1"))
      quantities_ |= GYOTO_QUANTITY_USER1;
    else if (!strcmp(tk, "User2"))
      quantities_ |= GYOTO_QUANTITY_USER2;
    else if (!strcmp(tk, "User3"))
      quantities_ |= GYOTO_QUANTITY_USER3;
    else if (!strcmp(tk, "User4"))
      quantities_ |= GYOTO_QUANTITY_USER4;
    else if (!strcmp(tk, "User5"))
      quantities_ |= GYOTO_QUANTITY_USER5;
    else throwError("ScenerySubcontractor(): unknown quantity"); 
    tk = strtok(NULL, " \t\n");
  }

# ifdef HAVE_UDUNITS
  if (screen_) screen_ -> unmapPixUnit();
# endif

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "("<<squant<<"): " << "quantities_=" << quantities_ << endl;
# endif
}
Gyoto::Quantity_t Scenery::getRequestedQuantities() const {
  return quantities_?quantities_:(astrobj()?astrobj()->getDefaultQuantities():0);
}

void Scenery::intensityConverter(string unit) {
# ifdef HAVE_UDUNITS
  // default is SI
  if (unit=="") unit="J.m-2.s-1.sr-1.Hz-1";
  intensity_converter_ = new Units::Converter("J.m-2.s-1.sr-1.Hz-1", unit);
# else
  if (unit!="")
    GYOTO_WARNING
      << "Unit ignored, please recompile gyoto with --with-udunits"
      << endl;
# endif
}

void Scenery::setPropertyConverters(Astrobj::Properties * data) {
# ifdef HAVE_UDUNITS
  data -> intensityConverter(intensity_converter_);
  data -> spectrumConverter(spectrum_converter_);
  data -> binSpectrumConverter(binspectrum_converter_);
# endif
}

void Scenery::spectrumConverter(string unit) {
# ifdef HAVE_UDUNITS
  // default is SI
  if (unit=="") unit="J.m-2.s-1.sr-1.Hz-1";
  spectrum_converter_ = new Units::Converter("J.m-2.s-1.sr-1.Hz-1", unit);
# else
  if (unit!="")
    GYOTO_WARNING
      << "Unit ignored, please recompile gyoto with --with-udunits"
      << endl;
# endif
}

void Scenery::binSpectrumConverter(string unit) {
# ifdef HAVE_UDUNITS
  // default is SI
  if (unit=="") unit="J.m-2.s-1.sr-1";
  binspectrum_converter_ = new Units::Converter("J.m-2.s-1.sr-1", unit);
# else
  if (unit!="")
    GYOTO_WARNING
      << "Unit ignored, please recompile gyoto with --with-udunits"
      << endl;
# endif
}

std::string Scenery::getRequestedQuantitiesString() const {
  string squant = "";
  Quantity_t quantities
    = quantities_?quantities_:(astrobj()?astrobj()->getDefaultQuantities():0);
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
    = quantities_?quantities_:(astrobj()?astrobj()->getDefaultQuantities():0);
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

double Scenery::tMin() const { return ph_.tMin(); }
double Scenery::tMin(const string &unit) const {
  return ph_.tMin(unit);
}

void Scenery::tMin(double tmin) { ph_.tMin(tmin); }
void Scenery::tMin(double tmin, const string &unit) {
  ph_.tMin(tmin, unit);
}

void Scenery::adaptive(bool mode) { ph_.adaptive(mode); }
bool Scenery::adaptive() const { return ph_.adaptive(); }

void Scenery::secondary(bool sec) { ph_.secondary(sec); }
bool Scenery::secondary() const { return ph_.secondary(); }

void Scenery::maxiter(size_t miter) { ph_.maxiter(miter); }
size_t Scenery::maxiter() const { return ph_.maxiter(); }

#ifdef GYOTO_USE_XERCES
void Scenery::fillElement(FactoryMessenger *fmp) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "fmp -> metric (metric()) ;" << endl;
# endif
  if (metric())     fmp -> metric (metric()) ;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG <<"fmp -> screen (screen_) ;" << endl;
# endif
  if (screen_) fmp -> screen (screen_) ;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG <<"fmp -> astrobj (astrobj()) ;" << endl;
# endif
  if (astrobj())    fmp -> astrobj (astrobj()) ;

  if (delta_ != GYOTO_DEFAULT_DELTA) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG <<"fmp -> setParameter (\"Delta\", "<<delta_<<") ;" << endl;
#   endif
    fmp -> setParameter ("Delta", delta_);
  }

  if (!adaptive()) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG <<"fmp -> setParameter (\"NonAdaptive\") ;" << endl;
#   endif
    fmp -> setParameter ("NonAdaptive");
  }

  if (maxiter_ != GYOTO_DEFAULT_MAXITER)
    fmp -> setParameter("MaxIter", maxiter_);

  if (getRequestedQuantities()) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG <<"fmp -> setParameter (\"Quantities\", \""
		<<getRequestedQuantitiesString()<<"\") ;" << endl;
#   endif
    fmp -> setParameter("Quantities", getRequestedQuantitiesString());
  }

  if (tMin() != -DBL_MAX) fmp -> setParameter("MinimumTime", tMin());
  if (nthreads_) fmp -> setParameter("NThreads", nthreads_);
}

SmartPointer<Scenery> Gyoto::Scenery::Subcontractor(FactoryMessenger* fmp) {

  string name="", content="", unit="";
  SmartPointer<Metric::Generic> gg = NULL;
  SmartPointer<Screen> scr = NULL;
  SmartPointer<Astrobj::Generic> ao = NULL;
  string squant = "";

  gg = fmp->metric();
  scr= fmp->screen();
  ao = fmp->astrobj();

  SmartPointer<Scenery> sc = new Scenery(gg, scr, ao);

  while (fmp->getNextParameter(&name, &content, &unit)) {
    char* tc = const_cast<char*>(content.c_str());
    if (name=="Delta")       sc -> delta(atof(tc), unit);;
    if (name=="Quantities")  sc -> setRequestedQuantities(tc);
    if (name=="MinimumTime") sc -> tMin(atof(tc), unit);
    if (name=="NThreads")    sc -> nThreads(atoi(tc));
    if (name=="MaxIter")     sc -> maxiter(atoi(tc));
    if (name=="Adaptive")    sc -> adaptive(true);
    if (name=="NonAdaptive") sc -> adaptive(false);
    if (name=="PrimaryOnly") sc -> secondary(false);
  }

  return sc;
}
#endif
