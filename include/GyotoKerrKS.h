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
  bool   generic_integrator_; ///< which integrator to use

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
  void genericIntegrator(bool);
  bool genericIntegrator() const ;
  
  double gmunu(const double * x,
		       int alpha, int beta) const ;

  void gmunu(double g[4][4], const double * pos) const;

  /**
   *\brief The inverse matrix of gmunu
   */ 
  void gmunu_up(double gup[4][4], const double * pos) const;

  /**
   * \brief The derivatives of gmunu
   *
   * Used in the test suite
   */
  void jacobian(double dst[4][4][4], const double * x) const ;

  using Generic::christoffel;
  int christoffel(double dst[4][4][4], const double * x) const ;
  int christoffel(double dst[4][4][4], const double * pos, double gup[4][4], double jac[4][4][4]) const ;

  void nullifyCoord(double coord[8], double &tdot2) const;
  void nullifyCoord(double coord[8]) const;
  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

 public:

  void MakeCst(const double* coord, double* cst) const;
  ///< In Kerr-Schild coordinates [T,x,y,z,Tdot,xdot,ydot,zdot], computes the four constants of the movement : particule mass, energy, angular momentum and Carter's constant.
 protected:

  /**
   * \brief RK4 integrator
   *
   * Wrapper around myrk4(const double * coord, const double* cst , double h, double* res) const
   *
   *
   */
  int myrk4(Worldline * line, const double coord[8], double h, double res[8]) const;//NB non adaptive integration doesn't work for KS ; this function is not implemented

  /**
   * \brief RK4 integrator
   * \param coord [r,theta,phi,t,pr,ptheta]
   * \param cst   [a,E,L,Q],dy/dtau=F(y,cst)
   * \param h     proper time step.
   * \param res   result
   */
  int myrk4(const double * coord, const double* cst , double h, double* res) const;

  /**
   * \brief ?
   *
   * Is it ever called?
   */
  int myrk4_adaptive(Gyoto::Worldline* line, const double * coord, double lastnorm, double normref, double* coord1, double h0, double& h1, double h1max) const;

  /** F function such as dy/dtau=F(y,cst)
   */
  using Generic::diff;
  int diff(const double* coord, const double* cst, double* res) const;
  virtual int isStopCondition(double const * const coord) const;
  void setParticleProperties(Worldline* line, const double* coord) const;

};

#endif
