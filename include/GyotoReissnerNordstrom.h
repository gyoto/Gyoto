/**
 *  \file GyotoReissnerNordstrom.h
 *  \brief Reissner-Nordstrom charged spherically sym BH (or naked singularity) spacetime.
 */

/*
    Copyright 2026 Frederic Vincent & Thibaut Paumard

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

#ifndef __GyotoReissnerNordstrom_h
#define __GyotoReissnerNordstrom_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class ReissnerNordstrom;
  };
};

class Gyoto::Metric::ReissnerNordstrom
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::ReissnerNordstrom>;

protected:
  double charge_; ///< Dimensionless charge (if <=1 gives a BH spacetime, if >1 gives a naked singularity spacetime)
  
public:
  GYOTO_OBJECT;
  ReissnerNordstrom();
  ReissnerNordstrom(const ReissnerNordstrom & orig);
  virtual ~ReissnerNordstrom();
  virtual ReissnerNordstrom * clone() const ;

  void charge(const double charge); ///< Set charge
  double charge() const ; ///< Returns charge

  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const ;
  double gmunu_up(double const x[4], int mu, int nu) const ;
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], double const pos[4]) const ;
  int isStopCondition(double const coord[8]) const;
  void circularVelocity(double const * coor, double* vel, double dir) const;
  double getPotential(double const pos[4], double l_cst) const;
  double getSpecificAngularMomentum(double rr) const;
#endif
};
