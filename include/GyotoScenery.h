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
 * - ImpactR, ImpactX, ImpactY, ImpactZ: R, X, Y, and Z of emission point;
 * - Spectrum: I_{nu} computed at various values frequencies,
 *        corresponding to the Screen's Spectrometer.
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
  double getDelta() const ;
  void setDelta(double);

  void setRequestedQuantities(Quantity_t) ;
  void setRequestedQuantities(std::string) ;
  Quantity_t getRequestedQuantities() const ;
  std::string getRequestedQuantitiesString() const ;
  size_t getScalarQuantitiesCount() const ;

  // Worker:
 public:
  void rayTrace(size_t imin, size_t imax, size_t jmin, size_t jmax,
		Astrobj::Properties* data, int save=0);

  void operator() (size_t i, size_t j, Astrobj::Properties *data);

#ifdef GYOTO_USE_XERCES
 public:
    void fillElement(FactoryMessenger *fmp); /// < called from Factory
#endif
 
};

#ifdef GYOTO_USE_XERCES
namespace Gyoto {
  SmartPointer<Scenery> ScenerySubcontractor(Gyoto::FactoryMessenger*);
}
#endif

#endif
