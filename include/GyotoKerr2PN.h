/**
 *  \file GyotoKerr2PN.h
 *  \brief Kerr spacetime at second post-Newtonian order
 *         in Cartesian harmonic coordinates (t,x,y,z).
 */

/*
    Copyright 2026 Karim Abd El Dayem, Frederic Vincent & Thibaut Paumard

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

#ifndef __GyotoKerr2PN_h
#define __GyotoKerr2PN_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class Kerr2PN;
  };
};

class Gyoto::Metric::Kerr2PN
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::Kerr2PN>;

protected:
  double spin_; ///< Geometrized units spin parameter
  double Q_; ///< Geometrized units mass quadrupole parameter
  bool quadrupole_is_kerr_; ///< Boolean: True iff the quadrupole is tied to its Kerr value -spin^2 in geometrized units
  
public:
  GYOTO_OBJECT;
  Kerr2PN();
  Kerr2PN(const Kerr2PN & orig);
  virtual ~Kerr2PN();
  virtual Kerr2PN * clone() const ;

  void spin(const double charge); ///< Sets spin
  double spin() const ; ///< Returns spin
  void quadrupole(const double charge); ///< Sets quadrupole
  double quadrupole() const ; ///< Returns quadrupole
  void kerrQuadrupole(bool t) ; ///< Sets boolean quadrupole_is_kerr_
  bool kerrQuadrupole() const ; ///< Returns boolean quadrupole_is_kerr_

  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const ;
  using Generic::gmunu_up;
  //double gmunu_up(double const x[4], int mu, int nu) const ;
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], double const pos[4]) const ;
  int isStopCondition(double const coord[8]) const;
  //void circularVelocity(double const * coor, double* vel, double dir) const;
  //double getPotential(double const pos[4], double l_cst) const;
  //double getSpecificAngularMomentum(double rr) const;
#endif
};
