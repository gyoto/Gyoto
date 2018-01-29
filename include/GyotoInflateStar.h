/**
 * \file GyotoInflateStar.h
 * \brief Mass-less, spherical object following a timelike geodesic
 *
 *  A Gyoto::InflateStar evolves in a Gyoto::Metric following time-like
 *  geodesics and is a Gyoto::Astrobj::Generic suitable for
 *  ray-tracing.
 */

/*
    Copyright 2011, 2013, 2018 Frederic Vincent, Thibaut Paumard

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


#ifndef __GyotoInflateStar_H_ 
#define __GyotoInflateStar_H_ 

namespace Gyoto{
  namespace Astrobj { class InflateStar; }
}

#include <GyotoMetric.h>
#include <GyotoStar.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::InflateStar
 * \brief An Astrobj::Star with growing size
 *
 */
class Gyoto::Astrobj::InflateStar :
  public Gyoto::Astrobj::Star {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::InflateStar>;
  
  // Data : 
  // -----
 private:
  double timeinflateinit_; ///< coordinate time of starting inflation
  double timeinflatestop_; ///< coordinate time of stopping inflation
  double radiusstop_; ///< maximum radius of star
  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

 /**
  * Create InflateStar object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  InflateStar(); ///< Default constructor
  
  InflateStar(const InflateStar& orig); ///< Copy constructor
  virtual InflateStar * clone() const ;

  virtual ~InflateStar() ;                        ///< Destructor
  
  GYOTO_OBJECT_ACCESSORS_UNIT(timeInflateInit);
  GYOTO_OBJECT_ACCESSORS_UNIT(timeInflateStop);
  GYOTO_OBJECT_ACCESSORS_UNIT(radiusStop);

  using Star::radius;
  virtual double radiusAt(double t) const; ///< Radius at a given time
  virtual double radiusAt(double t,
			  const std::string &t_unit) const; ///< Radius at a given time
  virtual double radiusAt(double t,
			  const std::string &t_unit,
			  const std::string &r_unit) const; ///< Radius at a given time

  // Accessors
  // ---------
 public:
  virtual std::string className() const ; ///< "InflateStar"
  virtual std::string className_l() const ; ///< "inflate_star"

 public:
  
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL);
  virtual double emission(double nu_em, double dsem,
			  state_t const &cp, double const co[8]=NULL) const;
};


#endif
