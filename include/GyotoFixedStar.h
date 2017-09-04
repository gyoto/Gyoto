/**
 * \file GyotoFixedStar.h
 * \brief Fixed (i.e. non-moving) star
 *
 *  The target of ray-traced Gyoto::Photon
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


#ifndef __GyotoFixedStar_H_ 
#define __GyotoFixedStar_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class FixedStar; }
}

#include <GyotoUniformSphere.h>
#include <GyotoMetric.h>

/**
 * \class Gyoto::Astrobj::FixedStar. 
 * \brief Fixed (i.e. non-moving) star (or spherical blob)
 *
 *  The target of ray-traced Gyoto::Photon
 */
class Gyoto::Astrobj::FixedStar : public Astrobj::UniformSphere {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::FixedStar>;

 // Data : 
  // -----

 protected:
  
  double pos_[3];///< x, y, z or r, theta, phi
  bool rotating_; ///< Whether the blob has a non-zero 3-velocity at Metric::circularVelocity(); note that the star anyway is fixed, but it can be useful to give it a non-zero velocity (that will enter the redshift computation) to simulate a snapshot on a rotating star

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  
  /**
   * Everything is undefined, call setCoordSys(), setPos() and
   * radius().
   */
  FixedStar();///< Default constructor

  FixedStar(const FixedStar& orig);///< Copy constructor
  virtual FixedStar* clone() const;

  FixedStar(SmartPointer<Gyoto::Metric::Generic> gg, double StPsn[3], double radius);
                   ///< Standard constructor
  
  virtual ~FixedStar() ;                        ///< Destructor
  
 public:
  // Accessors
  // ---------
 public:
  double const * getPos() const; ///< Get const pointer to pos_
  void getPos(double* dst) const; ///< Get a copy of the pos_ array

  std::vector<double> position() const; ///< Get vector copy of #pos_
  void position(std::vector<double> const&); ///< Set #pos_ from vector

  double rMax();

  bool rotating() const ;
  void rotating(bool) ;

  using Generic::metric;
  virtual void metric(SmartPointer<Metric::Generic> metric) ;
  using UniformSphere::radius;
  virtual void radius(double radius); ///< Set radius
  void setPos(const double[3]); ///< Set pos_ array
  //  void setCoordSys(int); ///< set coordinate system
  
  // Outputs
  // -------
 protected:
  virtual void getCartesian(double const * const dates, size_t const n_dates,
		double * const x, double * const y,
		double * const z, double * const xprime=NULL,
		double * const yprime=NULL,  double * const zprime=NULL) ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;


};


#endif
