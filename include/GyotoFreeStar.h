/**
 * \file FreeStar.h
 * \brief UniformShere which follow a user-defined orbit with a constant speed.
 */

/*
    Copyright 2019 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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


#ifndef __GyotoFreeStar_H_ 
#define __GyotoFreeStar_H_ 

namespace Gyoto{
  namespace Astrobj { class FreeStar; }
}

#include <iostream>
#include <fstream>
#include <iomanip>
#include <GyotoMetric.h>
#include <GyotoUniformSphere.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::FreeStar
 * \brief UniformShere following a trajectory specified in getVelocity (non-geodesic) with a constant velocity
 *
 */

class Gyoto::Astrobj::FreeStar :
  public Gyoto::Astrobj::UniformSphere{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::FreeStar>;
  
  // Data : 
  // -----
 private:
  double* posIni_; // 4-position of the star in spherical coordinates
  double* fourveldt_; // 4-velocity of the star in spherical coordinates (dxi/dt, not dtau)
  bool posSet_;

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

 /**
  * Create FreeStar object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. initCoord()
  * can be used for that.
  */
  FreeStar(); ///< Default constructor
  
  FreeStar(const FreeStar& orig); ///< Copy constructor
  virtual FreeStar * clone() const ;

  virtual ~FreeStar() ;                        ///< Destructor
  
 public:
  virtual std::string className() const ; ///< "FreeStar"
  virtual std::string className_l() const ; ///< "free_star"

 public:
  void initPosition(std::vector<double> const &v);
  std::vector<double> initPosition() const;
  void initVelocity(std::vector<double> const &v);
  std::vector<double> initVelocity() const;
  void initCoord(std::vector<double> const &v);
  std::vector<double> initCoord() const;

  void getCartesian(double const * const dates, size_t const n_dates,
          double * const x, double * const y,
          double * const z, double * const xprime=NULL,
          double * const yprime=NULL,
          double * const zprime=NULL);

  void getVelocity(double const pos[4], double vel[4]);

};


#endif