/**
 * \file GyotoStar.h
 * \brief Mass-less, spherical object following time geodesics
 *
 *  A Gyoto::Star evolves in a Gyoto::Metric following time-like
 *  geodesics and is a Gyoto::Astrobj suitable for ray-tracing.
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
  class Star;
}

#include <GyotoMetric.h>
#include <GyotoAstrobj.h>
#include <GyotoSpectrum.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Star
 * \brief A mass-less, spherical object
 * 
 * Gyoto can compute the Star's orbit in a Gyoto::Metric and to
 * perform ray-tracing on this target.
 */
class Gyoto::Star : public Gyoto::Astrobj, public Gyoto::Worldline {
  friend class Gyoto::SmartPointer<Gyoto::Star>;
  
  // Data : 
  // -----
 protected:
  double radius_ ; ///< star radius
  SmartPointer<Spectrum::Generic> spectrum_;
  SmartPointer<Spectrum::Generic> opacity_;

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
  Star(SmartPointer<Metric> gg, double radius,
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

  virtual void setMetric(SmartPointer<Metric>);
  virtual SmartPointer<Metric> getMetric() const;
  virtual void setSpectrum(SmartPointer<Spectrum::Generic>);
  virtual SmartPointer<Spectrum::Generic> getSpectrum() const;
  virtual void setOpacity(SmartPointer<Spectrum::Generic>);
  virtual SmartPointer<Spectrum::Generic> getOpacity() const;

  /**
   * The mass of a Star is always 1. Stars do not perturb the
   * metric. The only relevant point is that Stars are massive
   * particules, their exact mass is of no importance.
   */
  virtual double getMass() const ; ///< Return 1.

 public:
  virtual double getRmax();
  virtual void unsetRmax();
  double getRadius() const ; ///< Get radius_
  void   setRadius(double); ///< Set radius_
  //  void setCoordSys(int); ///< Get coordinate system for integration
  //  int  getCoordSys(); ///< Set coordinate system for integration
  void setInitialCondition(double coord[8]); ///< Same as Worldline::setInitialCondition(gg, coord, sys,1)

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ; /// < called from Factory
  static Astrobj::Subcontractor_t Subcontractor;
  static void Init();
#endif

  virtual double operator()(double const coord[4]) ;

  /**
   * Star::Impact() computes the distance of closest approach between
   * the Photon and the Astrobj (if data->distance is non-NULL). The
   * value stored is quadratic, it is not exactly the distance.
   *
   * There is a (false) assumption that the star's trajectory is
   * linear and uniform in Cartesian coordinates between to
   * integration points. Since the integration happens in spherical
   * coordinates, this is not true.
   */
  virtual int Impact_(Photon *ph, size_t index, AstrobjProperties* data=NULL);
 protected:
  virtual void getVelocity(double const pos[4], double vel[4]) ;

  virtual double emission(double nu_em, double dsem,
			  double cp[8], double co[8]=NULL) const;
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const;
  virtual double transmission(double nuem, double dsem, double*) const ;

};


#endif
