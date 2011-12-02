/**
 * \file GyotoStar.h
 * \brief Mass-less, spherical object following a timelike geodesic
 *
 *  A Gyoto::Star evolves in a Gyoto::Metric following time-like
 *  geodesics and is a Gyoto::Astrobj::Generic suitable for
 *  ray-tracing.
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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


#ifndef __GyotoStar_H_ 
#define __GyotoStar_H_ 

namespace Gyoto{
  namespace Astrobj { class Star; }
}

#include <GyotoMetric.h>
#include <GyotoUniformSphere.h>
#include <GyotoSpectrum.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::Star
 * \brief Mass-less, spherical object following a timelike geodesic
 *
 * Gyoto can compute the Star's orbit in a Gyoto::Metric and perform
 * ray-tracing on this target. The XML description of a Star looks
 * like:
\code
<Astrobj kind = "Star">
  <Metric kind = "KerrBL">
    <Spin> 0. </Spin>
  </Metric>
  <Radius> 2. </Radius>
  <Velocity> 0. 0. 0.037037 </Velocity>
  <Position> 600. 9. 1.5707999999999999741 0 </Position>
  <Spectrum kind="BlackBody">
    <Temperature> 6000 </Temperature>
  </Spectrum>
  <Opacity kind="PowerLaw">
    <Exponent> 0 </Exponent>
    <Constant> 0.1 </Constant>
  </Opacity>
  <OpticallyThin/>
</Astrobj>
\endcode
 * 
 * The Metric element can be of any kind. This Metric sets the
 * coordinate system.
 *
 * The Star is a coordinate sphere of radius Radius in solid motion.
 *
 * Position sets the initial 4-coordinate of the centre of the
 * sphere. Velocity contains its initial 3-velocity (the time
 * derivatives of the 3 space coordinates).
 *
 * Like many Astrobj::Generic impementations, a Star can be
 * OpticallyThin or OpticallyThick.
 *
 * Spectrum and Opacity (if OpticallyThin) are the descriptions of two
 * Gyoto::Spectrum::Generic sub-classes.
 *
 */
class Gyoto::Astrobj::Star :
  public Gyoto::Astrobj::UniformSphere,
  public Gyoto::Worldline {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Star>;
  
  // Data : 
  // -----

  // Constructors - Destructor
  // -------------------------
 public:
 /**
  * Create Star object and set initial condition.
  * \param gg: Gyoto::SmartPointer to the Gyoto::Metric in this part of the Universe
  * \param radius star radius
  * \param pos[4] initial position
  * \param v[3] initial velocity
  */
  Star(SmartPointer<Metric::Generic> gg, double radius,
       double pos[4], double v[3]) ;                        ///< Standard constructor

 /**
  * Create Star object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  Star(); ///< Default constructor
  
  Star(const Star& orig); ///< Copy constructor
  virtual Star * clone() const ;

  virtual ~Star() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual std::string className() const ; ///< "Star"
  virtual std::string className_l() const ; ///< "star"

  virtual void setMetric(SmartPointer<Metric::Generic>);
  virtual SmartPointer<Metric::Generic> getMetric() const;

  /**
   * The mass of a Star is always 1. Stars do not perturb the
   * metric. The only relevant point is that Stars are massive
   * particules, their exact mass is of no importance.
   */
  virtual double getMass() const ; ///< Return 1.

 public:
  virtual double getRmax();
  virtual void unsetRmax();
  //  void setCoordSys(int); ///< Get coordinate system for integration
  //  int  getCoordSys(); ///< Set coordinate system for integration
  void setInitialCondition(double coord[8]); ///< Same as Worldline::setInitialCondition(gg, coord, sys,1)

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ; ///< called from Factory
  static Astrobj::Subcontractor_t Subcontractor;
  static void Init();
#endif

 public:
  virtual void getCartesian(double const * const dates, size_t const n_dates,
		double * const x, double * const y,
		double * const z, double * const xprime=NULL,
		double * const yprime=NULL,  double * const zprime=NULL) ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;

};


#endif
