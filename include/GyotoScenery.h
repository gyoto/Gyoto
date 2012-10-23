/**
 * \file GyotoScenery.h
 * \brief Ray-tracing framework
 *
 *  A Metric, an Astrobj and a screen.
 */

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

#ifndef __GyotoScenery_H_ 
#define __GyotoScenery_H_ 

namespace Gyoto{
  class Scenery;
}

#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>
#include <GyotoAstrobj.h>
#include <GyotoMetric.h>
#include <GyotoScreen.h>
#include <GyotoPhoton.h>
#include <GyotoConverters.h>

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
 * relevant (it is not stored). A value of the integration step for
 * the Photon's trajectory can be specified in Delta. It will be used
 * as the initial step, which is adaptive. Possible Quantities:
 *
 * - Intensity: the intensity that reaches the object, integrated over
 *        the line-of-sight;
 * - EmissionTime: date of emission;
 * - MinDistance: minimum distance between the Photon reaching each
 *        pixel and the Astrobj;
 * - FirstDistMin: last closest approach between Photon and Astrobj;
 * - Redshift;
 * - ImpactCoords: 8-coordinates of the object and photon at impact;
 * - Spectrum: I_{nu} computed at various values frequencies,
 *        corresponding to the Screen's Spectrometer.
 *
 * In addition, it is possible to ray-trace an image using several
 * cores on a single machine (if Gyoto has been compiled with POSIX
 * threads support). The number of threads can be specified using
 * NThreads entity. Setting NThreads to 0 is equivalent to setting it
 * to 1. Beware that setting NThreads to a number higher than the
 * actual number of cores available on the machine usually leads to a
 * decrease in performance.
 *
 * Thus a fully populated Scenery XML looks like that:
 * \code
 * <?xml version="1.0" encoding="UTF-8" standalone="no"?>
 * <Scenery>
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
 *  <Delta> 1. </Delta>
 *
 *  <NThreads> 2 </NThreads>  
 *
 * </Scenery>
 * \endcode
 */
class Gyoto::Scenery : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Scenery>;
  
  
  // Data : 
  // -----
 protected:
  /**
   * The Metric, or stage, for this scenery.
   */
  SmartPointer<Metric::Generic> gg_;

  /**
   * Screen, the camera for this scenery.
   */
  SmartPointer<Screen> screen_;

  /**
   * The astrophysical emitting light in this scenery... the actor.
   */
  SmartPointer<Astrobj::Generic> obj_;

  /**
   * Default integration step for the photons
   */
  double delta_; // default integration step for the photons

  /**
   * The list of quantities that will be computed, for instance:
   * GYOTO_QUANTITY_INTENSITY | GYOTO_QUANTITY_EMISSIONTIME | ...
   */
  Gyoto::Quantity_t quantities_;

  /**
   * Used internally to not always reallocate memory when operator() is called.
   */
  Gyoto::Photon ph_; ///< a Photon.

  /**
   * Computation does not go back before tmin_. Default is -DBL_MAX. tmin_ is
   * always expressed in geometrical units, it is essentially a tuning
   * parameter for the ray-tracing process. tmin should be chosen to
   * always be longer than the distance between the screen and the
   * object.
   */
  double tmin_; ///< Time limit for the integration (geometrical units)

  /**
   * When compiled with libpthread, Scenery::rayTrace() may compute
   * several points of the image in parallel threads. This is the
   * number of threads to use.
   */
  size_t nthreads_; ///< number of parallel threads to use in ::rayTrace

# ifdef HAVE_UDUNITS
  Gyoto::SmartPointer<Gyoto::Units::Converter> intensity_converter_;
  Gyoto::SmartPointer<Gyoto::Units::Converter> spectrum_converter_;
  Gyoto::SmartPointer<Gyoto::Units::Converter> binspectrum_converter_;
# endif

  // Constructors - Destructor
  // -------------------------
 public:
  Scenery(); ///< Set everything to defaults
  Scenery (const Scenery& o); ///< Copy constructor
  Scenery * clone() const; ///< Cloner

  /**
   * To ensure consistency, the Metric will be forcibly attached to
   * the Screen and to the Astrobj.
   */
  Scenery(SmartPointer<Metric::Generic>, SmartPointer<Screen>, SmartPointer<Astrobj::Generic>);
  
  ~Scenery();

  // Mutators / assignment
  // ---------------------
 public:
  // Accessors
  // ---------
  SmartPointer<Metric::Generic> getMetric(); ///< Get Metric
  /**
   * The provided Metric will also be atached to the Screen and the Astrobj.
   */
  void setMetric(SmartPointer<Metric::Generic>);  ///< Set Metric
  SmartPointer<Screen> getScreen(); ///< Get Screen object

  /**
   * The Metric attached to the Scenery will be attached to the Screen
   */
  void setScreen(SmartPointer<Screen>);///< Set screen object
  SmartPointer<Astrobj::Generic> getAstrobj();
  /**
   * The Metric attached to the Scenery will be attached to the Astrobj
   */
  void setAstrobj(SmartPointer<Astrobj::Generic>);
  double getDelta() const ; ///< get default step in geometrical units
  double getDelta(const std::string &unit) const ;  ///< get default step in specified units
  void setDelta(double); ///< set default step in geometrical units
  void setDelta(double, const std::string &unit);   ///< set default step in specified units

  void setRequestedQuantities(Quantity_t) ;
  void setRequestedQuantities(std::string) ;
  Quantity_t getRequestedQuantities() const ;
  std::string getRequestedQuantitiesString() const ;
  size_t getScalarQuantitiesCount() const ;

  double getTmin() const ;///< get tmin_
  double getTmin(const std::string &unit) const ;///< get tmin_
  void setTmin(double); ///< set tmin_;
  void setTmin(double, const std::string &unit); ///< set tmin_;

  void setNThreads(size_t); ///< set nthreads_;
  size_t getNThreads() const ; ///< get nthreads_;

  void setIntensityConverter(std::string unit);
  void setSpectrumConverter(std::string unit);
  void setBinSpectrumConverter(std::string unit);
  void setPropertyConverters(Gyoto::Astrobj::Properties *);

  // Worker:
 public:
  void rayTrace(size_t imin, size_t imax, size_t jmin, size_t jmax,
		Astrobj::Properties* data, double * impactcoords = NULL);

  void operator() (size_t i, size_t j, Astrobj::Properties *data,
		   double * impactcoords = NULL, Photon * ph = NULL);

#ifdef GYOTO_USE_XERCES
 public:
    void fillElement(FactoryMessenger *fmp); ///< called from Factory
#endif
 
};

#ifdef GYOTO_USE_XERCES
namespace Gyoto {
  SmartPointer<Scenery> ScenerySubcontractor(Gyoto::FactoryMessenger*);
}
#endif

#endif
