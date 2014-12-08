/*
    Copyright 2011, 2013-2014 Thibaut Paumard, Frederic Vincent

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

#ifdef HAVE_MPI
#include "GyotoFactory.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/intercommunicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <string>
#include <boost/serialization/string.hpp>
#include <boost/serialization/array.hpp>
namespace mpi = boost::mpi;
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include <sys/time.h>    /* for benchmarking */


using namespace Gyoto;
using namespace std;

/// Properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Scenery)
GYOTO_WORLDLINE_PROPERTY_END(Scenery, Object::properties)

///

Scenery::Scenery() :
  screen_(NULL), delta_(GYOTO_DEFAULT_DELTA),
  quantities_(0), ph_(), nthreads_(0), nprocesses_(0)
#ifdef HAVE_MPI
  , mpi_team_(NULL)
#endif
{}

Scenery::Scenery(SmartPointer<Metric::Generic> met,
		 SmartPointer<Screen> scr,
		 SmartPointer<Astrobj::Generic> obj) :
  screen_(scr), delta_(GYOTO_DEFAULT_DELTA),
  quantities_(0), ph_(), nthreads_(0), nprocesses_(0)
#ifdef HAVE_MPI
  , mpi_team_(NULL)
#endif
{
  metric(met);
  if (screen_) screen_->metric(met);
  astrobj(obj);
}

Scenery::Scenery(const Scenery& o) :
  SmartPointee(o),
  screen_(NULL), delta_(o.delta_), 
  quantities_(o.quantities_), ph_(o.ph_), 
  nthreads_(o.nthreads_), nprocesses_(0)
#ifdef HAVE_MPI
  , mpi_team_(NULL)
#endif
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
# ifdef HAVE_MPI
  if (!Scenery::am_worker) mpiTerminate();
# endif
 }

SmartPointer<Metric::Generic> Scenery::metric() const { return ph_.metric(); }

void Scenery::metric(SmartPointer<Metric::Generic> met) {
  ph_.metric(met);
  if (!screen_) screen_ = new Screen ();
  screen_ -> metric(met);
}

SmartPointer<Screen> Scenery::screen() const { return screen_; }

void Scenery::screen(SmartPointer<Screen> scr) {
  screen_ = scr;
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

void Scenery::initCoord(std::vector<double> c) { ph_ . initCoord(c); }
std::vector<double> Scenery::initCoord() const { return ph_.initCoord();}

void  Scenery::nThreads(size_t n) { nthreads_ = n; }
size_t Scenery::nThreads() const { return nthreads_; }

void  Scenery::nProcesses(int n) { nprocesses_ = n; }
int Scenery::nProcesses() const { return nprocesses_; }

typedef struct SceneryThreadWorkerArg {
#ifdef HAVE_PTHREAD
  pthread_mutex_t * mutex;
  pthread_t * parent;
#endif
  Screen::Coord2dSet & ij;
  size_t cnt, npix;
  Scenery *sc;
  Photon * ph;
  Astrobj::Properties *data;
  double * impactcoords;
  SceneryThreadWorkerArg(Screen::Coord2dSet & ijin);
  bool is_pixel;
} SceneryThreadWorkerArg ;

SceneryThreadWorkerArg::SceneryThreadWorkerArg(Screen::Coord2dSet & ijin)
  :ij(ijin)
{

}


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
  GYOTO_ARRAY<size_t, 2> ijb;
  GYOTO_ARRAY<double, 2> ad;

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

    // copy i & j or alpha and delta

    if (larg->is_pixel) ijb =  *(larg->ij);
    else ad = larg->ij.angles();

    if (!larg->ij.valid()) {
      // terminate, but first...
#ifdef HAVE_PTHREAD
      // ...unlock mutex so our siblings can access i & j and terminate too
      if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif
      break;
    }

    size_t lcnt = larg->cnt++;
 
    ++(larg->ij);

#ifdef HAVE_PTHREAD
    // unlock mutex so our siblings can can access i, j et al. and procede
    if (larg->mutex) pthread_mutex_unlock(larg->mutex);
#endif

    // Store current cell
    data = *larg->data; 
    size_t cell=lcnt;
    if (larg->is_pixel && data.alloc) cell=(ijb[1]-1)*larg->npix+ijb[0]-1;
    data += cell;
    impactcoords=larg->impactcoords?larg->impactcoords+16*cell:NULL;

    if (larg->is_pixel)
      (*larg->sc)(ijb[0], ijb[1], &data, impactcoords, ph);
    else (*larg->sc)(ad[0], ad[1], &data, ph);
    
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

void Scenery::rayTrace(Screen::Coord2dSet & ij,
		       Astrobj::Properties *data,
		       double * impactcoords) {
  GYOTO_DEBUG_EXPR(am_worker);

#if defined HAVE_MPI
  if (nprocesses_ && !mpi_team_) {
    mpiSpawn(nprocesses_);
    mpiClone();
  }
#endif

  /*
     Ray-trace now is multi-threaded. What it does is
       - some initialization
       - launch nthreads_ - 1  thread working on of SceneryThreadWorker
       - call SceneryThreadWorker itself rather than sleeping
       - wait for the other threads to be terminated
       - some housekeeping
   */

  if (!screen_) {
    if (am_worker) throwError("No screen, have you called mpiClone()?");
    else throwError("Scenery::rayTrace() needs a Screen to work on");
  }
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
  screen_ -> getRayCoord(size_t(1),size_t(1), coord);
  ph_ . setInitCoord(coord, -1);
  // delta is reset in operator()

  if (data) setPropertyConverters(data);

  GYOTO_ARRAY<size_t, 2> ijb;
  GYOTO_ARRAY<double, 2> ad;

  bool alloc=data?data->alloc:false;
  size_t npix=screen_->resolution();

#ifdef HAVE_MPI
  if (mpi_team_) {
    // We are in an MPI content, either the manager or a worker.
    // dispatch over workers and monitor

    if (!am_worker) {
      mpi_tag tag=raytrace;
      mpiTask(tag);
    }

    size_t nbnuobs=0;
    Quantity_t quantities = (am_worker || !data)?GYOTO_QUANTITY_NONE:*data;
    bool has_ipct=am_worker?false:bool(impactcoords);
    bool is_pixel=(ij.kind==Screen::pixel);
    
    mpi::broadcast(*mpi_team_, quantities, 0);
    mpi::broadcast(*mpi_team_, has_ipct, 0);
    mpi::broadcast(*mpi_team_, is_pixel, 0);

    if (quantities & (GYOTO_QUANTITY_SPECTRUM | GYOTO_QUANTITY_BINSPECTRUM)) {
      if (!spr) throwError("Spectral quantity requested but "
			     "no spectrometer specified!");
      nbnuobs = spr -> nSamples();
    }
    size_t nelt= getScalarQuantitiesCount(&quantities);
    if (quantities & GYOTO_QUANTITY_SPECTRUM)     nelt += nbnuobs;
    if (quantities & GYOTO_QUANTITY_BINSPECTRUM)  nelt += nbnuobs;
    if (quantities & GYOTO_QUANTITY_IMPACTCOORDS) nelt += 16;
    double * vect = new double[nelt];
    Astrobj::Properties *locdata = new Astrobj::Properties();
    size_t offset=1;
    size_t curquant=0;

    if (Scenery::am_worker) {
      // set all converters to the trivial one, conversion is
      // performed in the manager.
      intensity_converter_ = NULL;
      spectrum_converter_ = NULL;
      binspectrum_converter_ = NULL;
      setPropertyConverters(locdata);
    }

    if (quantities & GYOTO_QUANTITY_INTENSITY) {
      locdata->intensity=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_EMISSIONTIME) {
      locdata->time=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) {
      locdata->distance=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_FIRST_DMIN) {
      locdata->first_dmin=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_REDSHIFT) {
      locdata->redshift=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_IMPACTCOORDS) {
      locdata->impactcoords=vect+offset*curquant; curquant+=16;
    }
    if (quantities & GYOTO_QUANTITY_USER1) {
      locdata->user1=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_USER2) {
      locdata->user2=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_USER3) {
      locdata->user3=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_USER4) {
      locdata->user4=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_USER5) {
      locdata->user5=vect+offset*(curquant++);
    }
    if (quantities & GYOTO_QUANTITY_SPECTRUM) {
      locdata->spectrum=vect+offset*(curquant++);
      locdata->offset=int(offset);
    }
    if (quantities & GYOTO_QUANTITY_BINSPECTRUM) {
      locdata->binspectrum=vect+offset*(curquant++);
      locdata->offset=int(offset);
    }

    mpi::status s;
    if (!am_worker) { // We are the manager
      int working = mpi_team_->size()-1;
      // First tell the workers to join our task force
      // The corresponding recv is in gyoto-scenery-worker.c

      vector<size_t> cell(working+1);
      size_t cnt=0;

      while (working) {
	// receive one result, need to track back where it belongs and
	// store it there
	int w;

	// Wait for worker to ask for task.
	// tag may be raytrace_done if worker has result to report, 
	// give_task if worker has no data yet.
	s = mpi_team_ -> recv(mpi::any_source, mpi::any_tag, vect, nelt);

	w = s.source();
	
	size_t cs=cell[w]; // remember where to store results

	// give new task or decrease working counter
	if (ij.valid()) {
	  cell[w]=cnt++;
	  if (is_pixel) {
	    ijb=*ij;
	    if (alloc) cell[w]=(ijb[1]-1)*npix+ijb[0]-1;
	    mpi_team_ -> send(w, raytrace, ijb);
	  } else {
	    ad = ij.angles();
	    mpi_team_ -> send(w, raytrace, ad);
	  }
	  if (impactcoords) {
	    mpi_team_ -> send(w, Scenery::impactcoords, impactcoords+cell[w]*16, 16);
	  }

	  ++ij;

	} else {
	  if (is_pixel)
	    mpi_team_ -> send(w, raytrace_done, ijb);
	  else
	    mpi_team_ -> send(w, raytrace_done, ad);
	  --working;
	}

	// Now that the worker is back to work, triage data it has just delivered

	if (s.tag()==Scenery::raytrace_done && data) {
	  // Copy each relevant quantity, performing conversion if needed
	  if (data->intensity)
	    data->intensity[cs]=data->intensity_converter_?
	      (*data->intensity_converter_)(*locdata->intensity):
	      *locdata->intensity;
	  if (data->time) data->time[cs]=*locdata->time;
	  if (data->distance) data->distance[cs]=*locdata->distance;
	  if (data->first_dmin) data->first_dmin[cs]=*locdata->first_dmin;
	  if (data->redshift) data->redshift[cs]=*locdata->redshift;
	  if (data->impactcoords)
	    for (size_t k=0; k<16; ++k)
	      data->impactcoords[cs*16+k]=locdata->impactcoords[k];
	  if (data->user1) data->user1[cs]=*locdata->user1;
	  if (data->user2) data->user2[cs]=*locdata->user2;
	  if (data->user3) data->user3[cs]=*locdata->user3;
	  if (data->user4) data->user4[cs]=*locdata->user4;
	  if (data->user5) data->user5[cs]=*locdata->user5;
	  if (data->spectrum)
	    for (size_t c=0; c<nbnuobs; ++c)
	      data->spectrum[cs+c*data->offset]=
		data->spectrum_converter_?
		(*data->spectrum_converter_)(locdata->spectrum[c]):
		locdata->spectrum[c];
	  if (data->binspectrum)
	    for (size_t c=0; c<nbnuobs; ++c)
	      data->binspectrum[cs+c*data->offset]=
		data->binspectrum_converter_?
		(*data->binspectrum_converter_)(locdata->binspectrum[c]):
		locdata->binspectrum[c];
	}

      }
      if (verbose()) cout << endl;
    } else {
      // We are a worker, do we need to query for impactcoords?
      double ipct[16];
      if (has_ipct) impactcoords=&ipct[0];

      // First send dummy result, using tag "give_tag".
      // Manager will ignore the results and send first coordinates.
      mpi_team_->send(0, give_task, vect, nelt);
      while (true) {
	// Receive new coordinates to work on.
	if (is_pixel)
	  s = mpi_team_->recv(0, mpi::any_tag, ijb);
	else
	  s = mpi_team_->recv(0, mpi::any_tag, ad);
	if (s.tag()==raytrace_done) {
	  break;
	}
	// Receive impactcoords if needed
	if (has_ipct)
	  s = mpi_team_->recv(0, Scenery::impactcoords, impactcoords, 16);
	locdata->init(nbnuobs);
	if (is_pixel)
	  (*this)(ijb[0], ijb[1], locdata, impactcoords, &ph_);
	else
	  (*this)(ad[0], ad[1], locdata, &ph_);
	// send result
	mpi_team_->send(0, raytrace_done, vect, nelt);
      }
    }
    delete locdata;
    delete [] vect;
    return;
  }
#endif

  SceneryThreadWorkerArg larg(ij);
  larg.sc=this;
  larg.ph=&ph_;
  larg.data=data;
  larg.cnt=0;
  larg.npix=npix;
  larg.impactcoords=impactcoords;
  larg.is_pixel= (ij.kind==Screen::pixel);

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

  GYOTO_MSG << "\nRaytraced "<< ij.size()
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

void Scenery::operator() (
			  double a, double d,
			  Astrobj::Properties *data,
			  Photon *ph
			  ) {

  double coord[8];
  SmartPointer<Spectrometer::Generic> spr = screen_->spectrometer();
  size_t nbnuobs = spr() ? spr -> nSamples() : 0;

  if (data) data -> init(nbnuobs);

  if (!ph) {
    ph = &ph_;
    ph -> spectrometer(spr);
    ph -> freqObs(screen_->freqObs());
  }
  // Always reset delta
  ph -> delta(delta_);

  screen_ -> getRayCoord(a, d, coord);
  ph -> setInitCoord(coord, -1);
  ph -> hit(data);

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

size_t Scenery::getScalarQuantitiesCount(Quantity_t *q) const {
  size_t nquant=0;
  Quantity_t quantities;
  if (q) quantities=*q;
  else
    quantities=quantities_?
      quantities_:
      (astrobj()?astrobj()->getDefaultQuantities():GYOTO_QUANTITY_NONE);
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

void Scenery::integrator(std::string type) {ph_.integrator(type);}
std::string Scenery::integrator() const { return ph_.integrator();}
double Scenery::deltaMin() const {return ph_.deltaMin();}
double Scenery::deltaMax() const {return ph_.deltaMax();}
void  Scenery::deltaMin(double h1) {ph_.deltaMin(h1);}
void  Scenery::deltaMax(double h1) {ph_.deltaMax(h1);}
double Scenery::deltaMaxOverR() const { return ph_.deltaMaxOverR();}
void Scenery::deltaMaxOverR(double t) {ph_.deltaMaxOverR(t);}

double Scenery::absTol() const {return ph_.absTol();}
void Scenery::absTol(double t) {ph_.absTol(t);}
double Scenery::relTol() const {return ph_.relTol();}
void Scenery::relTol(double t) {ph_.relTol(t);}

void Scenery::secondary(bool sec) { ph_.secondary(sec); }
bool Scenery::secondary() const { return ph_.secondary(); }

void Scenery::maxiter(size_t miter) { ph_.maxiter(miter); }
size_t Scenery::maxiter() const { return ph_.maxiter(); }

#ifdef GYOTO_USE_XERCES
void Scenery::fillElement(FactoryMessenger *fmp) {
  if (metric())     fmp -> metric (metric()) ;
  if (screen_)      fmp -> screen (screen_) ;
  if (astrobj())    fmp -> astrobj (astrobj()) ;

  if (getRequestedQuantities())
    fmp -> setParameter("Quantities", getRequestedQuantitiesString());

  fmp -> setParameter("MinimumTime", tMin());
  fmp -> setParameter("NThreads", nthreads_);
  fmp -> setParameter("NProcesses", nprocesses_);

  Object::fillElement(fmp);
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

  int mpi=0;

  while (fmp->getNextParameter(&name, &content, &unit)) {
    char* tc = const_cast<char*>(content.c_str());
    if (name=="Quantities")  sc -> setRequestedQuantities(tc);
    else if (name=="MinimumTime") sc -> tMin(atof(tc), unit);
    else if (name=="NThreads")    sc -> nThreads(atoi(tc));
    else if (name=="NProcesses")  sc -> nProcesses(atoi(content.c_str()));
    else if (name=="Metric" || name=="Screen" || name=="Astrobj") ;
    else sc -> setParameter(name, content, unit);
  }
  return sc;
}
#endif

bool Gyoto::Scenery::am_worker=false;

void Gyoto::Scenery::mpiSpawn(int nbchildren) {
  nprocesses_=nbchildren;
#ifdef HAVE_MPI
  int flagi=0, flagt=0;
  if (MPI_Initialized(&flagi)) throwError("Error running MPI_Initialized()");
  if (!flagi) {
    GYOTO_WARNING << "MPI_Init not called yet: "
		  << "not spawning processes" << endl;
    return;
  }
  if (MPI_Finalized(&flagt)) throwError("Error running MPI_Finalized()");
  if (flagt) {
    GYOTO_WARNING << "MPI_Finalize already called: "
		  << "not spawning processes" << endl;
    return;
  }
  if (mpi_team_) {
    if (mpi_team_->size()==nbchildren+1) return;
    mpiTerminate();
  }
  if (nbchildren) {
    char * exec = const_cast<char*>("gyoto-mpi-worker." GYOTO_SOVERS); 

    MPI_Comm children_c;
    MPI_Comm_spawn(exec,
		   MPI_ARGV_NULL, nbchildren,
		   MPI_INFO_NULL, 0, MPI_COMM_SELF, &children_c,
		   MPI_ERRCODES_IGNORE);

    mpi_team_ = new mpi::communicator(mpi::intercommunicator (children_c, mpi::comm_take_ownership).merge(false));
  }
#else
  GYOTO_WARNING << "No MPI in this Gyoto" << endl;
#endif 
}

void Gyoto::Scenery::mpiTerminate() {
#ifdef HAVE_MPI
  if (mpi_team_) {
    mpi_tag tag=terminate;
    mpiTask(tag);
    mpi_team_->barrier();
    delete mpi_team_;
    mpi_team_=NULL;
  }
#endif
}

void Gyoto::Scenery::mpiClone() {
#ifdef HAVE_MPI
  if (!mpi_team_) return;
  std::string xmldata=
    Gyoto::Factory(this).format();
  int errcode;
  mpi_tag tag=read_scenery;
  mpiTask(tag);
  broadcast(*mpi_team_, xmldata, 0);
#endif
}

void Gyoto::Scenery::mpiTask(mpi_tag &tag) {
#ifdef HAVE_MPI
  if (!mpi_team_) return;
  mpi::broadcast(*mpi_team_, tag, 0);
#endif
}

