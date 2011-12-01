/**
 * \file GyotoStar.h
 * \brief Mass-less, spherical object following a timelike geodesic
 *
 *  A Gyoto::UniformSphere evolves in a Gyoto::Metric following time-like
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


#ifndef __GyotoUniformSphere_H_ 
#define __GyotoUniformSphere_H_ 

namespace Gyoto{
  namespace Astrobj { class UniformSphere; }
}

#include <GyotoMetric.h>
#include <GyotoAstrobj.h>
#include <GyotoSpectrum.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::UniformSphere
 * \brief Mass-less, spherical object following a timelike geodesic
 *
 * Gyoto can compute the UniformSphere's orbit in a Gyoto::Metric and perform
 * ray-tracing on this target. The XML description of a UniformSphere looks
 * like:
\code
<Astrobj kind = "UniformSphere">
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
 * The UniformSphere is a coordinate sphere of radius Radius in solid motion.
 *
 * Position sets the initial 4-coordinate of the centre of the
 * sphere. Velocity contains its initial 3-velocity (the time
 * derivatives of the 3 space coordinates).
 *
 * Like many Astrobj::Generic impementations, a UniformSphere can be
 * OpticallyThin or OpticallyThick.
 *
 * Spectrum and Opacity (if OpticallyThin) are the descriptions of two
 * Gyoto::Spectrum::Generic sub-classes.
 *
 */
class Gyoto::Astrobj::UniformSphere :
  public Gyoto::Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::UniformSphere>;
  
  // Data : 
  // -----
 protected:
  double radius_ ; ///< star radius
  SmartPointer<Spectrum::Generic> spectrum_; ///< star emission law
  SmartPointer<Spectrum::Generic> opacity_; ///< if optically thin, opacity law

  // Constructors - Destructor
  // -------------------------
 public:
 /**
  * Create UniformSphere object and set initial condition.
  * \param gg: Gyoto::SmartPointer to the Gyoto::Metric in this part of the Universe
  * \param radius star radius
  */
  UniformSphere(std::string kind,
		SmartPointer<Metric::Generic> gg, double radius) ;
      ///< Standard constructor

 /**
  * Create UniformSphere object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  UniformSphere(std::string kind); ///< Default constructor
  
  UniformSphere(const UniformSphere& orig); ///< Copy constructor

  virtual ~UniformSphere() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual std::string className() const ; ///< "UniformSphere"
  virtual std::string className_l() const ; ///< "uniformsphere"

  virtual void setSpectrum(SmartPointer<Spectrum::Generic>);
  virtual SmartPointer<Spectrum::Generic> getSpectrum() const;
  virtual void setOpacity(SmartPointer<Spectrum::Generic>);
  virtual SmartPointer<Spectrum::Generic> getOpacity() const;

 public:
  double getRadius() const ; ///< Get radius_
  void   setRadius(double); ///< Set radius_

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ; /// < called from Factory
  virtual void setGenericParameters(FactoryMessenger *fmp) ;
  /// < Interpret common XML sections
#endif

  virtual double operator()(double const coord[4]) ;

 protected:
  virtual void getCartesian(double const * const dates, size_t const n_dates,
		double * const x, double * const y,
		double * const z, double * const xprime=NULL,
		double * const yprime=NULL,  double * const zprime=NULL) =0;

  virtual void getVelocity(double const pos[4], double vel[4]) = 0;

  virtual double emission(double nu_em, double dsem,
			  double cp[8], double co[8]=NULL) const;
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const;
  virtual double transmission(double nuem, double dsem, double*) const ;

};


#endif
