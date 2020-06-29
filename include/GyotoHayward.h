/**
 * \file GyotoHayward.h
 * \brief Metric of a regular rotating black hole or naked worm-hole
 *
 * This is a regular rotating extension of Hayward's metric.
 *
 * See Lamy et al. (2018), Classical and Quantum Gravity, submitted,
 * https://arxiv.org/abs/1802.01635
 */

/*
    Copyright 2018 Frederic Lamy, Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoHayward_H_
#define __GyotoHayward_H_ 

namespace Gyoto {
  namespace Metric { class Hayward; }
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif


/**
 * \class Gyoto::Metric::Hayward
 * \brief Metric of a regular rotating black hole or naked worm-hole
 *
 * This is a regular rotating extension of Hayward's metric. 
 *
 * The metric reads: ds^2=-(1-2*M(r)*r/Sigma)*dt^2-(4*a*M(r)*sin^2(theta)/Sigma)*dt*dphi+(Sigma/Delta)*dr^2+(Sigma)*dtheta^2+sin^2(theta)*[r^2+a^2+2*a^2*r*M(r)*sin^2(theta)/Sigma]*dphi^2,
 * where: 
 * Sigma=r^2+a^2*cos^2(theta)
 * Delta=r^2-2*M(r)*r+a^2
 * M(r)=M*abs(r)^3/(abs(r)^3+2*m*b^2), m being the mass of the black hole seen at infinity and b a magnetic charge parameter.
 *
 * See Lamy et al. (2018), Classical and Quantum Gravity, submitted,
 * https://arxiv.org/abs/1802.01635
 *
 * See Hayward (2006) for the original nonrotating metric,
 * https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.031103
 *
 * See Bambi & Modesto (2013) for a rotating (but singular) extension.
 * https://www.sciencedirect.com/science/article/pii/S0370269313002505?via%3Dihub
 */
class Gyoto::Metric::Hayward : public Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::Hayward>;
  
  // Data : 
  // -----
 protected:
  double charge_; ///< Magnetic charge parameter
  double spin_ ;  ///< Angular momentum parameter
  double a2_ ; ///< spin_*spin_
  double a3_ ; ///< a2_*spin_
  double a4_ ; ///< a2_*a2_
  double b2_; ///< charge_*charge_
  
  // Constructors - Destructor
  // -------------------------
 public: 
  GYOTO_OBJECT;
  Hayward(); ///< Default constructor
  virtual Hayward * clone () const ;

  // Accessors
  // ---------
  void spin(const double spin); ///< Set spin
  double spin() const ; ///< Returns spin
  void charge(const double charge); ///< Set charge
  double charge() const ; ///< Returns charge

  // Methods needed for PolishDoughnut
  virtual double getSpecificAngularMomentum(double rr) const;
  virtual double getPotential(double const pos[4], double l_cst) const;

  // Keplerian equatorial orbits angular velocity
  virtual void circularVelocity(double const coor[4], double vel[4],
				  double dir) const;

  // Actual space-time API
  void gmunu(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;
  double gmunu(double const x[4], int mu, int nu) const ;
  void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;
  double gmunu_up(double const x[4], int mu, int nu) const ;
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], const double pos[4]) const ;
  
  // Optimized
  double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const ;

};

#endif

