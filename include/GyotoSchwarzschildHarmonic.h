/**
 *  \file GyotoSchwarzschildHarmonic.h
 *  \brief Schwarzschild spacetime in harmonic coordinates.
 */

/*
    Copyright 2021 Frederic Vincent & Thibaut Paumard

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

#ifndef __GyotoSchwarzschildHarmonic_h
#define __GyotoSchwarzschildHarmonic_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class SchwarzschildHarmonic;
  };
};

class Gyoto::Metric::SchwarzschildHarmonic
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::SchwarzschildHarmonic>;
 public:
  GYOTO_OBJECT;
  SchwarzschildHarmonic();
  SchwarzschildHarmonic(const SchwarzschildHarmonic & orig);
  virtual ~SchwarzschildHarmonic();
  virtual SchwarzschildHarmonic * clone() const ;

  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const ;
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], double const pos[4]) const ;
  int isStopCondition(double const coord[8]) const;
  void circularVelocity(double const * coor, double* vel, double dir) const;
    
#endif
};
