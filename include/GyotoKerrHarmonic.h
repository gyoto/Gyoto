
/**
 *  \file GyotoKerrHarmonic.h
 *  \brief Kerr spacetime in Cartesian harmonic coordinates (t,x,y,z),
 *         defined from spherical Kerr-Schild coordinates as derived
 *         by Cook & Scheel (1997).
 */

/*
    Copyright 2026 Frederic Vincent, Thibaut Paumard & Karim Abd El Dayem

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

#ifndef __GyotoKerrHarmonic_h
#define __GyotoKerrHarmonic_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class KerrHarmonic;
  };
};

class Gyoto::Metric::KerrHarmonic
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::KerrHarmonic>;

protected:
  double spin_; ///< Geometrized units spin parameter
  
public:
  GYOTO_OBJECT;
  KerrHarmonic();
  KerrHarmonic(const KerrHarmonic & orig);
  virtual ~KerrHarmonic();
  virtual KerrHarmonic * clone() const ;

  void spin(const double charge); ///< Sets spin
  double spin() const ; ///< Returns spin

  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const ;
  
  using Generic::gmunu_up;
  void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const;

  void jacobian(double ARGOUT_ARRAY3[4][4][4], const double x[4]) const ;
  
  int isStopCondition(double const coord[8]) const;
  void circularVelocity(double const * coor, double* vel, double dir) const;
  //double getPotential(double const pos[4], double l_cst) const;
  //double getSpecificAngularMomentum(double rr) const;
#endif
};
