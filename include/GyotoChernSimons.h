/**
 *  \file GyotoChernSimons.h
 *  \brief Chern-Simons 1st order perturbation to KerrBL metric
 *
 */

/*
    Copyright 2013 Frederic Vincent

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

#ifndef __GyotoChernSimons_h
#define __GyotoChernSimons_h

#include <GyotoKerrBL.h>

namespace Gyoto {
  namespace Metric {
    class ChernSimons;
  };
};

class Gyoto::Metric::ChernSimons
: public Gyoto::Metric::KerrBL {
  friend class Gyoto::SmartPointer<Gyoto::Metric::ChernSimons>;
 protected:
  double dzetaCS_; ///< Chern-Simons coupling constant
 public:
  GYOTO_OBJECT;
  ChernSimons();
  ChernSimons(const ChernSimons &o);
  virtual ~ChernSimons();
  virtual ChernSimons * clone() const ;

  void dzetaCS(double d);
  double dzetaCS() const;

  double gmunu(const double * const x, int mu, int nu) const ;
  double gmunu_up(const double * const x, int mu, int nu) const ;
  int diff(const double y[8], const double cst[5], double res[8]) const ;
  void circularVelocity(double const pos[4], double vel [4],
			double dir=1.) const ;
};
#endif
