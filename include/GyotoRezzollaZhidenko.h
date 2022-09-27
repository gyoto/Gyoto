/**
 *  \file GyotoRezzollaZhidenko.h
 *  \brief Spherically-symmetric parametrized metric of Rezzolla\&Zhidenko 2014
 *         See the paper: PRD, 90, 084009
 *         Only epsilon, a0, a1, a2, a3, b0, b1, b2, b3 are allowed non-zero
 */

/*
    Copyright 2013, 2018 Frederic Vincent & Thibaut Paumard

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

#ifndef __GyotoRezzollaZhidenko_h
#define __GyotoRezzollaZhidenko_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class RezzollaZhidenko;
  };
};

class Gyoto::Metric::RezzollaZhidenko
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::RezzollaZhidenko>;
 private:
  double epsilon_; ///< horizon parameter, rH=2/(1+eps)
  double rms_, rmb_; ///< Provide marginally stable and bound orbits if needed
  double* aparam_; ///< The a-parameter vector [a0,a1,a2,a3] used in RZ14
  double* bparam_; ///< The b-parameter vector [b0,b1,b2,b3] used in RZ14
 public:
  GYOTO_OBJECT;
  RezzollaZhidenko();
  RezzollaZhidenko(const RezzollaZhidenko & orig);
  virtual ~RezzollaZhidenko();
  virtual RezzollaZhidenko * clone() const ;

  // accessors
  GYOTO_OBJECT_ACCESSORS(double, epsilon);
  GYOTO_OBJECT_ACCESSORS(double, rms);
  GYOTO_OBJECT_ACCESSORS(double, rmb);
  void aparam(std::vector<double> const &v);
  std::vector<double> aparam() const;
  void bparam(std::vector<double> const &v);
  std::vector<double> bparam() const;


  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const ;
  double N2(const double rr) const;
  double B2(const double rr) const;
  double Nprime(const double rr) const;
  double Bprime(const double rr) const;
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], double const pos[4]) const ;
  int isStopCondition(double const coord[8]) const;
  virtual double getRmb() const;
  virtual double getRms() const;
  virtual double getPotential(double const pos[4], double l_cst) const;
  virtual double getSpecificAngularMomentum(double rr) const;
  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

    
#endif
};
