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
    Copyright 2011-2014, 2018-2019 Frederic Vincent, Thibaut Paumard

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
 *  Gyoto::Astrobj::UniformSphere::setParameters() take care of
 *  interpreting the XML elements describing the parameters of the
 *  sphere:
 *  \code
 *     <Radius> value </Radius>
 *     <Spectrum kind="..."> parameters for this spectrum kind </Spectrum>
 *     <Opacity kind="..."> parameters for this spectrum kind </Opacity>
 *
 *     The following are numerical parameters mostly usefull when the
 *     sphere is far from the compact object. Larger values speed up
 *     computation but may miss the sphere.
 *     <DeltaMaxOverRadius> 0.1 </DeltaMaxOverRadius>
 *     <DeltaMaxOverDistance> 0.1 </DeltaMaxOverDistance>
 *  \endcode
 * setGenericParameters() also takes care of calling
 * setParameter().
 */
class Gyoto::Astrobj::UniformSphere :
  public Gyoto::Astrobj::Standard {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::UniformSphere>;
  
  // Data : 
  // -----
 protected:
  double radius_ ; ///< sphere radius [geometrical units]
  bool isotropic_; ///< if 1, then emission just returns 1
  SmartPointer<Spectrum::Generic> spectrum_; ///< sphere emission law
  SmartPointer<Spectrum::Generic> opacity_; ///< if optically thin, opacity law

  double dltmor_; ///< see deltaMax(double*)
  double dltmod_; ///< see deltaMax(double*)

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

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
  * Create UniformSphere object. Use metric(), radius(),
  * spectrum() and opacity() to set the members.
  * 
  * \param kind: specify kind (e.g. "Star" or "FixedStar")
  */
  UniformSphere(std::string kind); ///< Default constructor
  
  UniformSphere(const UniformSphere& orig); ///< Copy constructor

  virtual ~UniformSphere() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual std::string className() const ; ///< "UniformSphere"
  virtual std::string className_l() const ; ///< "uniformsphere"

  virtual void spectrum(SmartPointer<Spectrum::Generic>);
  ///< Set spectrum_
  virtual SmartPointer<Spectrum::Generic> spectrum() const;
  ///< Get spectrum_
  virtual void opacity(SmartPointer<Spectrum::Generic>);
  ///< Set opacity_
  virtual SmartPointer<Spectrum::Generic> opacity() const;
  ///< Get opacity_
  double radius() const ; ///< Get radius_ in geometrical units
  virtual void   radius(double); ///< Set radius_ in geometrical units
  double radius(std::string const &) const ; ///< Get radius_ in specified unit
  virtual void   radius(double, std::string const &); ///< Set radius_ in specified unit

  double deltaMaxOverRadius() const ; ///< Get dltmor_
  virtual void   deltaMaxOverRadius(double f); ///< Set dltmor_

  double deltaMaxOverDistance() const ; ///< Get dltmod_
  virtual void   deltaMaxOverDistance(double f); ///< Set dltmod_

  bool isotropic() const;
  void isotropic(bool);
  double alpha() const ;
  void alpha(double);

 public:

  virtual double operator()(double const coord[4]) ;
  ///< Square distance to the center of the sphere

  ///< Ensure integration does not miss the object
  /**
   * \param[in] coord current photon position
   * \return max( #dltmor_*#radius_, #dltmod_*operator()(double coord[]) )
   */
  virtual double deltaMax(double*coord);

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

  using Standard::emission;
  virtual double emission(double nu_em, double dsem,
			  state_t const &cp, double const co[8]=NULL) const;
  ///< Emission is determined by spectrum_ and opacity_
  using Standard::integrateEmission;
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   state_t const &c_ph, double const *c_obj=NULL) const;
  virtual double transmission(double nuem, double dsem, state_t const &, double const *) const ;
  ///< Transmission is determined by opacity_

};


#endif
