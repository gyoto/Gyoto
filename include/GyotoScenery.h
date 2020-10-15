/**
 * \file GyotoScenery.h
 * \brief Ray-tracing framework
 *
 *  A Metric, an Astrobj and a Screen.
 */

/*
    Copyright 2011-2016, 2018-2019 Thibaut Paumard

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

#ifndef __GyotoScenery_H_ 
#define __GyotoScenery_H_ 

namespace Gyoto{
  class Scenery;
}

#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>
#include <GyotoObject.h>
#include <GyotoAstrobj.h>
#include <GyotoMetric.h>
#include <GyotoScreen.h>
#include <GyotoPhoton.h>
#include <GyotoConverters.h>

#ifdef HAVE_MPI
#include "GyotoFactory.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif


/**
 * \class Gyoto::Scenery
 * \brief Ray-tracing scene
 *
 * An Scenery contains:
 *    - a Metric: used in Astrobj, Screen and Photon;
 *    - a Screen: sets the field-of-view, the position of the camera,
 *      the observation time, and the Spectrometer;
 *    - an Astrobj: light emitter.
 *
 *
 * In addition, Quantities may be specified (or the default Quantity
 * will be produced: generally Intensity). Not all Astrobj implement
 * all Quantities. The order in which Quantities are listed is not
 * relevant (it is not stored). Possible Quantities:
 *
 * - Intensity: the intensity that reaches the object, integrated over
 *        the line-of-sight;
 * - EmissionTime: date of emission;
 * - MinDistance: minimum distance between the Photon reaching each
 *        pixel and the Astrobj;
 * - FirstDistMin: last closest approach between Photon and Astrobj;
 * - Redshift;
 * - ImpactCoords: 8-coordinates of the object and photon at impact;
 * - Spectrum: I<SUB>&nu;</SUB> computed at various values frequencies,
 *        corresponding to the Screen's Spectrometer.
 * - BinSpectrum:
 *   &int;<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP>I<SUB>&nu;</SUB>d&nu;
 *   computed between various (&nu;<SUB>1</SUB>, &nu;<SUB>2</SUB>
 *   pairs corresponding to the Screen's Spectrometer. This is what a
 *   physical spectrometer measures.
 *
 * In addition, it is possible to ray-trace an image using several
 * cores on a single machine (if Gyoto has been compiled with POSIX
 * threads support). The number of threads can be specified using
 * NThreads entity. Setting NThreads to 0 is equivalent to setting it
 * to 1. Beware that setting NThreads to a number higher than the
 * actual number of cores available on the machine usually leads to a
 * decrease in performance.
 *
 * Finally, Scenery accepts a number of numerical tuning parameters
 * that are passed directly to the underlying photons (actually, the
 * Scenery object holds a Photon instance which stores many
 * parameters, including the Metric and Astrobj):
 * Adaptive/NonAdaptive, Delta, MinimumTime, MaxIter, PrimaryOnly.
 *
 * Thus a fully populated Scenery XML looks like that (the values are
 * examples, they are not necessary the default nor the best or even
 * good values):
 * \code <?xml version="1.0" encoding="UTF-8"
 * standalone="no"?> <Scenery>
 *
 *  <Metric kind = "MetricKind">
 *    <MetricProperties/>
 *  </Metric>
 *
 *  <Screen>
 *    <ScreenProperties/>
 *  </Screen>
 *
 *  <Astrobj kind = "AstrobjKind">
 *    <AstrobjParameters/>
 *  </Astrobj>
 *
 *  <Quantities> Spectrum Intensity ...</Quantities>
 *
 *  <NThreads> 2 </NThreads>  
 *
 *  Next come the numerical tuning parameters:
 *  Integration step, initial in case of adaptive, reset for
 *  for each ray being traced:
 *  <Delta unit="geometrical"> 1 </Delta>
 *
 *  Adaptive or NonAdaptive:
 *  <Adaptive/>
 *
 *  The integrator to use for integrating the photons:
 *  <Integrator>runge_kutta_fehlberg78</Integrator>
 *  The "Legacy" integrator is coded in
 *  Metric::Generic::myrk4_adaptive(), may be re-implemented in othe
 *  metrics, and therefore takes its tuning parameters in the Metric
 *  section. The other integrators (runge_kutta_fehlberg78,
 *  runge_kutta_cash_karp54, runge_kutta_dopri5,
 *  runge_kutta_cash_karp54_classic) accept the following tuning
 *  parameters, directly in the Scenery section:
 *
 *  Absolute and relative tolerance for the adaptive step:
 *  <AbsTol>1e-11</AbsTol>
 *  <RelTol>1e-11</RelTol>
 *  Normally, you should not need to tune the other three. If you need
 *  to, try using a higher order integrator:
 *  maximum integration step:
 *  <DeltaMax> 100 </DeltaMax>
 *  delta_max/R where R is the current distance to the origin:
 *  <DeltaMaxOverR> 0.1 </DeltaMaxOverR>
 *  minimum step:
 *  <DeltaMin>1e-20</DeltaMin>
 *
 *  A few safe-guards to avoid infinite loops:
 *
 *  Maximum number of iterations for each ray:
 *  <Maxiter> 1000000 </Maxiter>
 *
 *  Minimum date a photon may reach backwards in time:
 *  <MinimumTime unit="yr">25e3</MinimumTime>
 *
 *  This one is an experimental, poorly specified feature:
 *  <!--PrimaryOnly/-->
 *
 * </Scenery>
 * \endcode
 */
class Gyoto::Scenery
: public Gyoto::SmartPointee,
  public Gyoto::Object
{
  friend class Gyoto::SmartPointer<Gyoto::Scenery>;
  
  
  // Data : 
  // -----
 protected:

  /**
   * Screen, the camera for this scenery.
   */
  SmartPointer<Screen> screen_;


  /**
   * Default integration step for the photons
   */
  double delta_; // default integration step for the photons

  /// Quantities to compute
  /**
   * Bitwise OR of quantities that will be computed, for instance:
   * \code
   * GYOTO_QUANTITY_INTENSITY | GYOTO_QUANTITY_EMISSIONTIME | ...
   * \endcode
   */
  Gyoto::Quantity_t quantities_;

  /**
   * Used internally to not always reallocate memory when operator()
   * is called and to store all the parameters which affect the
   * integration, except delta_.
   */
  Gyoto::Photon ph_; ///< Template Photon.

  /**
   * When compiled with libpthread, Scenery::rayTrace() may compute
   * several points of the image in parallel threads. This is the
   * number of threads to use.
   */
  size_t nthreads_; ///< Number of parallel threads to use in rayTrace()

  int nprocesses_; ///< Number of parallel processes to use in rayTrace()

# ifdef HAVE_UDUNITS
  /// See Astrobj::Properties::intensity_converter_
  Gyoto::SmartPointer<Gyoto::Units::Converter> intensity_converter_;
  /// See Astrobj::Properties::intensity_converter_
  Gyoto::SmartPointer<Gyoto::Units::Converter> spectrum_converter_;
  /// See Astrobj::Properties::intensity_converter_
  Gyoto::SmartPointer<Gyoto::Units::Converter> binspectrum_converter_;
# endif

 public:
  GYOTO_OBJECT_THREAD_SAFETY;
# ifdef HAVE_MPI
  /// Team of processes for MPI
  /**
   * Rank 0 is the manager, other ranks are workers, instances of the
   * gyoto-mpi-worker executable.
   */
  boost::mpi::communicator * mpi_team_;
# endif
  /// True in instance of gyoto-mpi-worker, otherwise false.
  static bool am_worker;

  /// Spawn gyoto-mpi-worker processes
  /**
   * If nbchildren is -1 set #mpi_team_ to MPI_COMM_WORLD else spawn
   * nbchildren processes and set #nprocesses_ accordingly. If a
   * different number of workers are already running, terminate them
   * first. If nbchildren is 0, just terminate running workers.
   *
   * The approach of Gyoto to MPI is that a manager process (of rank 0
   * within a given MPI communicator) will distribute ray-tracing
   * tasks across worker processes. Several scenarii are supported,
   * including spawning instances of the gyoto-mpi-worker.version
   * executable, where "version" matches the version component in the
   * library name (typically a number, possibly followed by
   * "unreleased").
   *
   * In all cases, the manager process needs to call this function,
   * either with -1 if the worker processes are already running or
   * &gt;1 if workers need to be spawned.
   *
   * \param[in] nbchildren number of processes to spawn.
   */
  void mpiSpawn(int nbchildren);

  /// Terminate worker processes
  void mpiTerminate ();

  /// Send a copy of self to the mpi workers
  /**
   * Always call mpiClone() before ray-tracing if workers are running.
   */
  void mpiClone();

  /// Tags that may be sent to communicate with workers using MPI_Send()
  enum mpi_tag {give_task, read_scenery, terminate,
		raytrace, raytrace_done, ready,
		impactcoords, noimpactcoords};

  /// Send a tag to workers
  void mpiTask(mpi_tag &tag);

  /// Become an MPI worker
  /**
   * Worker processes need to call this function after having called
   * MPI_Init().
   */
  static void mpiWorker();

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_WORLDLINE;
  Scenery(); ///< Set everything to defaults
  Scenery (const Scenery& o); ///< Copy constructor
  Scenery * clone() const; ///< Cloner

  /// Constructor setting Scenery::gg_, Scenery::screen_, and Scenery::obj_ 
  /**
   * To ensure consistency, the Metric will be forcibly attached to
   * the Screen and to the Astrobj (if they are not NULL).
   */
  Scenery(SmartPointer<Metric::Generic>, SmartPointer<Screen>, SmartPointer<Astrobj::Generic>);
  
  ~Scenery();

  // Mutators / assignment
  // ---------------------
 public:
  // Accessors
  // ---------
  SmartPointer<Metric::Generic> metric() const; ///< Get ph_.Worldline::metric_
  /**
   * The provided Metric will also be atached to the Screen and the Astrobj.
   */
  void metric(SmartPointer<Metric::Generic>);  ///< Set Scenery::gg_
  SmartPointer<Screen> screen() const; ///< Get Scenery::screen_

  /**
   * The Metric attached to the Scenery will be attached to the Screen
   */
  void screen(SmartPointer<Screen>);///< Set Scenery::screen_
  SmartPointer<Astrobj::Generic> astrobj() const; ///< Get ph_.obj_
  /**
   * The Metric attached to the Scenery will be attached to the Astrobj
   */
  void astrobj(SmartPointer<Astrobj::Generic>); ///< Set ph_.obj_


  /**
   * \brief Get clone of template Photon
   */
  SmartPointer<Photon> clonePhoton() const; ///< Clone the internal Photon

  /**
   * \brief Get clone of template Photon, intitializing it to pixel
   */
  SmartPointer<Photon> clonePhoton(size_t i, size_t j);

  /**
   * \brief Get clone of template Photon, intitializing it to direction
   */
  SmartPointer<Photon> clonePhoton(double a, double d);

  void updatePhoton(); ///< Update values in cached Photon

  double delta() const ; ///< Get default step in geometrical units
  double delta(const std::string &unit) const ;  ///< Get default step in specified units
  void delta(double); ///< set default step in geometrical units
  void delta(double, const std::string &unit);   ///< set default step in specified units

  void initCoord(std::vector<double> c);
  std::vector<double> initCoord() const;


  /// Set Scenery::quantities_
  /**
   * \param quant Bitwise OR of desired quantities, e.g. \code GYOTO_QUANTITY_SPECTRUM | GYOTO_QUANTITY_MIN_DISTANCE \endcode
   */
  void setRequestedQuantities(Quantity_t quant) ;

  /// Set Scenery::quantities_ from string
  /**
   * \param squant Coma-separated list of quantities, e.g. "Spectrum
   * MinDistance". The order is not relevant.
   */
  void requestedQuantitiesString(std::string const &squant) ;

  /// Get Scenery::quantities_
  Quantity_t getRequestedQuantities() const ;

  /// Get a string representation of Scenery::quantities_
  std::string requestedQuantitiesString() const ;

  /// Get number of requested quantities of scalar nature
  /**
   * This is all quantities except Spectrum, BinSpectrum and ImpactCoords.
   */
  size_t getScalarQuantitiesCount(Quantity_t *q=NULL) const ;

  /// Get number of requested quantities of spectral nature
  /**
   * This is Spectrum, SpectrumStokesQ, SpectrumStokesU, SpectrumStokesV and BinSpectrum.
   */
  size_t getSpectralQuantitiesCount(Quantity_t *q=NULL) const ;

  /// Get ph_.tmin_
  double tMin() const ;
  /// Get ph_.tmin_ in specified unit
  double tMin(const std::string &unit) const ;
  /// Set ph_.tmin_
  void tMin(double);
  /// Set ph_.tmin_ in specified unit
  void tMin(double, const std::string &unit);

  void adaptive (bool mode) ; ///< Set ph_.adaptive_
  bool adaptive () const ; ///< Get ph_.adaptive_

  /// Passed to #ph_
  void integrator(std::string type);
  /// Passed to #ph_
  std::string integrator() const;

  /// Passed to #ph_
  double deltaMin() const;
  /// Passed to #ph_
  void deltaMin(double h1);

  /// Passed to #ph_
  double deltaMax() const;

  /// Passed to #ph_
  void deltaMax(double h1);

  /// Passed to #ph_
  double deltaMaxOverR() const;
  /// Passed to #ph_
  void deltaMaxOverR(double t);

  /// Passed to #ph_
  void absTol(double);
  /// Passed to #ph_
  double absTol()const;
  /// Passed to #ph_
  void relTol(double);
  /// Passed to #ph_
  double relTol()const;

  /// Passed to #ph_
  void maxCrossEqplane(double);
  /// Passed to #ph_
  double maxCrossEqplane()const;
  
  void secondary (bool sec) ; ///< Set ph_.secondary_
  bool secondary () const ; ///< Get ph_.secondary_

  void integ31(bool integ); ///< Set WorldlinIntegState integ_31_
  bool integ31() const ; ///< Get WorldlinIntegState integ_31_

  void parallelTransport (bool pt) ; ///< Set ph_.parallel_transport_
  bool parallelTransport () const ; ///< Get ph_.parallel_transport_

  void maxiter (size_t miter) ; ///< Set ph_.maxiter_
  size_t maxiter () const ; ///< Get ph_.maxiter_

  void nThreads(size_t); ///< Set nthreads_;
  size_t nThreads() const ; ///< Get nthreads_;

  void nProcesses(size_t); ///< Set nprocesses_;
  size_t nProcesses() const ; ///< Get nprocesses_;

  /// Set Scenery::intensity_converter_
  void intensityConverter(std::string unit);
  /// Set Scenery::spectrum_converter_
  void spectrumConverter(std::string unit);
  /// Set Scenery::binspectrum_converter_
  void binSpectrumConverter(std::string unit);

  /// Copy converters to Astrobj::Properties instance
  /**
   * Copy Scenery::intensity_converter_, Scenery::spectrum_converter_
   * and Scenery::binspectrum_converter_ to there alter ego in *prop.
   */
  void setPropertyConverters(Gyoto::Astrobj::Properties *prop);

  // Worker:
 public:
  /// Perform ray-tracing
  /**
   * For each directions specified, launch a Photon back in time to
   * compute the various quantities.
   *
   * At this time, the computed quantities depend on on the pointers
   * in *data which are not NULL.
   *
   * rayTrace() uses
   * - setPropertyConverters() to set the converters in *data;
   * - Astrobj::Properties::init() to initialize each cell in *data;
   * - Astrobj::Properties::operator++() to step through the arrays in *data.
   *
   * data must have been instantiated prior to calling rayTrace and
   * the various pointers in *data must be NULL or point to the first
   * cell in an array of size at least Screen::npix_ squared.
   *
   * If MPI support is built-in, MPI_Init() has been called, and
   * nprocesses_ is &ge;1, then rayTrace() will use several processes,
   * launching them using mpiSpawn() if necessary.
   *
   * Else, if Scenery::nthreads_ is &ge;2 and Gyoto has been compiled with
   * pthreads support, rayTrace() will use Scenery::nthreads_ threads
   * and launch photons in parallel. This works only if the
   * Astrobj::Generic::clone() and Metric::Generic::clone() methods
   * have been properly implemented for the specific astrobj and
   * metric kind, and if they are both thread-safe. At the moment,
   * unfortunately, Lorene metrics are known to not be thread-safe.
   *
   * \param[in] ij Screen::Coord2dSet specification of rays to trace. e.g.:
   *
   * \code
   * Screen::Range irange(imin, imax, di); 
   * Screen::Range jrange(jmin, jmax, dj);
   * Screen::Grid ij(irange, jrange); 
   * \endcode
   *
   * \param[in, out] data Pointer to a preallocated
   * Astrobj::Properties instance which sets which quantities must be
   * computed and where to store the output.
   *
   * \param[in] impactcoords Optional pointer to an array of
   * pre-computed impact coordinates. If impactcoords is provided,
   * rayTracing is skipped and the quantities in *data are fill
   * assuming that the impact coordinates are correct. This only makes
   * sense in optically thick mode, when ray-tracing several sceneries
   * for which the shape of the object is identical but their emission
   * distributions are not. impactcoords can be computed using the
   * ImpactCoords quantity.
   */
  void rayTrace(
#ifdef GYOTO_SWIGIMPORTED
		Coord2dSet & ij,
#else
		Screen::Coord2dSet & ij,
#endif
		Astrobj::Properties *data,
		double * impactcoords=NULL);

  /// Ray-trace a single pixel in Scenery::screen_
  /**
   * Almost identical to rayTrace(), but for a single pixel.
   *
   * If ph is passed, it is assumed to have been properly initialized
   * (with the right metric and astrobj etc.) already. Else, use
   * &Scenery::ph_.
   */
  void operator() (size_t i, size_t j, Astrobj::Properties *data,
		   double * impactcoords = NULL, Photon * ph = NULL);

  /// Ray-trace single direction
  /**
   * Almost identical to rayTrace(), but for a single direction.
   *
   * If ph is passed, it is assumed to have been properly initialized
   * (with the right metric and astrobj etc.) already. Else, use
   * &Scenery::ph_.
   */
  void operator() (double alpha, double delta, Astrobj::Properties *data,
		   Photon * ph = NULL);

#ifdef GYOTO_USE_XERCES
 public:
  // Override fillProperty() to issue InitCoord only if it was set
  void fillProperty(FactoryMessenger *fmp, Property const &p) const ;
  // Override fillElement to fill metric, screen and astrobj first
  void fillElement(FactoryMessenger *fmp) const;
  /// Instanciate Scenery from an XML description.
  static SmartPointer<Scenery> Subcontractor(Gyoto::FactoryMessenger*);

#endif
 
};

#endif
