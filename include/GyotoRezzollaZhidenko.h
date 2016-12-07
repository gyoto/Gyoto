/**
 *  \file GyotoRezzollaZhidenko.h
 *  \brief Spherically-symmetric parametrized metric of Rezzolla\&Zhidenko 2014
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
  double epsilon_;
  double rms_, rmb_;
  double* aparam_;
  double* bparam_;
 public:
  GYOTO_OBJECT;
  RezzollaZhidenko();
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


  double gmunu(const double * const x, int mu, int nu) const ;
  int diff(const double y[8], const double cst[5], double res[8]) const ;
  void circularVelocity(double const pos[4], double vel [4],
			double dir=1.) const ;
  double N2(const double rr) const;
  double B2(const double rr) const;
  double Nprime(const double rr) const;
  double Bprime(const double rr) const;
  int christoffel(double dst[4][4][4], const double * pos) const ;
  int isStopCondition(double const * const coord) const;
  double getRmb() const;
  double getRms() const;
  virtual double getPotential(double pos[4], double l_cst) const;
  virtual double getSpecificAngularMomentum(double rr) const;

    
#endif
};
