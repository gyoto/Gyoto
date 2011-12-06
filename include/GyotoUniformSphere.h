/**
 * \file GyotoUniformSphere.h
 * \brief Optically thick or thin, spherical objects
 *
 *  Gyoto::Astrobj::UniformSphere is an abstract type from which
 *  uniform, spherical objects inherit (in particular, the
 *  Gyoto::Astrobj::Star and Gyoto::Astrobj::FixedStar classes).
 *
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
#include <GyotoStandardAstrobj.h>
#include <GyotoSpectrum.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::UniformSphere
 * \brief Optically thick or thin, spherical objects
 *
 *  Gyoto::Astrobj::UniformSphere is an abstract type from which
 *  uniform, spherical objects inherit (in particular, the
 *  Gyoto::Astrobj::Star and Gyoto::Astrobj::FixedStar classes).
 *
 *  These objects are coordinate-spherical: they comprise all the
 *  points within a given radius from a centre. The distance is the
 *  usual Euclidian distance in a Cartesian coordinate system which is
 *  trivially determined by the coordinate system in which the Metric
 *  is expressed.
 *
 *  The sphere is in solid motion: all the points have the same
 *  4-velocity. The centre of the sphere may move. This motion and the
 *  velocity are provided by the derived classes through the
 *  getCartesian() and getVelocity() methods.
 *
 *  The spheres can be optically thick or optically thin. In the
 *  optically thin case, the opacity law provided as a Gyoto::Spectrum
 *  also sets both the emissivity. Another Gyoto::Spectrum provides
 *  the emission law of the source, which is uniform.
 *
 *  Gyoto::Astrobj::UniformSphere::setGenericParameters() take care of
 *  interpreting the XML elements describing the parameters of the
 *  sphere:
\code
   <Radius> value </Radius>
   <Spectrum kind="..."> parameters for this spectrum kind </Spectrum>
   <Opacity kind="..."> parameters for this spectrum kind </Opacity>
\endcode
 * setGenericParameters() also takes care of calling
 * Generic::setGenericParameter().
 */
class Gyoto::Astrobj::UniformSphere :
  public Gyoto::Astrobj::Standard {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::UniformSphere>;
  
  // Data : 
  // -----
 protected:
  double radius_ ; ///< sphere radius
  SmartPointer<Spectrum::Generic> spectrum_; ///< sphere emission law
  SmartPointer<Spectrum::Generic> opacity_; ///< if optically thin, opacity law

  // Constructors - Destructor
  // -------------------------
 public:
 /**
  * Create UniformSphere object.
  * \param kind: specifi kind (e.g. "Star" or "FixedStar")
  * \param gg: Gyoto::SmartPointer to the Gyoto::Metric in this part of the Universe
  * \param radius: sphere radius
  */
  UniformSphere(std::string kind,
		SmartPointer<Metric::Generic> gg, double radius) ;
      ///< Standard constructor

 /**
  * Create UniformSphere object. Use setMetric(), setRadius(),
  * setSpectrum() and setOpacity() to set the members.
  * 
  * \param kind: specifi kind (e.g. "Star" or "FixedStar")
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
  ///< Set spectrum_
  virtual SmartPointer<Spectrum::Generic> getSpectrum() const;
  ///< Get spectrum_
  virtual void setOpacity(SmartPointer<Spectrum::Generic>);
  ///< Set opacity_
  virtual SmartPointer<Spectrum::Generic> getOpacity() const;
  ///< Get opacity_
  double getRadius() const ; ///< Get radius_
  void   setRadius(double); ///< Set radius_

 public:
#ifdef GYOTO_USE_XERCES
  /**
   * The sub-classes implementations of the
   * Astrobj::Generic::fillElement() method should call
   * Astrobj::UniformSphere::fillElement() to fill the common bits.
   */
  virtual void fillElement(FactoryMessenger *fmp) const ;
  ///< Fill the generic XML bits
  /**
   * The sub-classes subcontractor function (see
   * Astrobj::Subcontractor_t) should call this after creating the
   * object to interpret the common bits (Spectrum, Opacity, Radius):
\code
  SmartPointer<MyObject> object = new MyObject (.....);
  fmp -> reset();
  object -> setGenericParameters(fmp);
\endcode
   */
  virtual void setGenericParameters(FactoryMessenger *fmp) ;
  ///< Interpret common XML sections
#endif

  virtual double operator()(double const coord[4]) ;
  ///< Square distance to the center of the sphere

 protected:
  /**
   * If the coordinate system of the Metric object is spherical, use a
   * trivial conversion.
   */
  virtual void getCartesian(double const * const dates, size_t const n_dates,
		double * const x, double * const y,
		double * const z, double * const xprime=NULL,
		double * const yprime=NULL,  double * const zprime=NULL) =0;
  ///< Yield the Cartesian coordinates of the center of the sphere

  virtual void getVelocity(double const pos[4], double vel[4]) = 0;
  ///< Yield velocity of the center of the sphere.

  virtual double emission(double nu_em, double dsem,
			  double cp[8], double co[8]=NULL) const;
  ///< Emission is determined by spectrum_ and opacity_
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const;
  virtual double transmission(double nuem, double dsem, double*) const ;
  ///< Transmission is determined by opacity_

};


#endif
