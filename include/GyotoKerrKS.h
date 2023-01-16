/**
 *  \file GyotoKerrKS.h
 *  \brief KerrKS metric
 *
 *  Warning: this metric is seldom used and may be buggy.
 */

/*
    Copyright 2011-2015 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoKerrKS_H_
#define __GyotoKerrKS_H_ 

namespace Gyoto {
  namespace Metric { class KerrKS; }
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

/**
 * \class Gyoto::Metric::KerrKS
 * \brief Metric around a Kerr black-hole in Kerr-Schild coordinates
 *  Warning: this metric is seldom used and may be buggy.
 *
 * By default, uses the generic integrator
 * (Metric::Generic::myrk4()). Use
\code
<SpecificIntegrator/>
\endcode
 * to use the specific integretor which is, as of writting, buggy.
 */
class Gyoto::Metric::KerrKS
: public Metric::Generic
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::KerrKS>;
  
  // Data : 
  // -----

 protected:
  double spin_ ;  ///< Angular momentum parameter
  double a2_;     ///< spin_*spin_
  double rsink_;  ///< numerical horizon
  double drhor_;  ///< horizon security

  // Constructors - Destructor
  // -------------------------
 public: 
  GYOTO_OBJECT;
  KerrKS(); ///< Default constructor
  virtual KerrKS* clone () const;         ///< Copy constructor
  
  // Mutators / assignment
  // ---------------------
 public:
  // default operator= is fine
  void spin(const double spin); ///< Set spin

  // Accessors
  // ---------
 public:
  double spin() const ; ///< Returns spin
  void horizonSecurity(double drhor);
  double horizonSecurity() const;
  
  virtual double gmunu(double const x[4], int alpha, int beta) const ;

  virtual void gmunu(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;


  using Gyoto::Metric::Generic::gmunu_up;
  /**
   *\brief The inverse matrix of gmunu
   */ 
  virtual void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const;

  /**
   * \brief The derivatives of gmunu
   *
   * Used in the test suite
   */
  virtual void jacobian(double ARGOUT_ARRAY3[4][4][4], const double x[4]) const ;

  virtual void gmunu_up_and_jacobian(double ARGOUT_ARRAY2[4][4], double ARGOUT_ARRAY3[4][4][4], const double IN_ARRAY1[4]) const;

  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

  virtual int isStopCondition(double const coord[8]) const;

  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit);
};

#endif
