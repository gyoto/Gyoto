/**
 *  \file GyotoXmasSKS.h
 *  \brief XmasSKS metric
 *
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoXmasSKS_H_
#define __GyotoXmasSKS_H_ 

namespace Gyoto {
  namespace Metric { class XmasSKS; }
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif


/// Default value for difftol_
#define GYOTO_XmasSKS_DEFAULT_DIFFTOL 1e-2

/**
 * \class Gyoto::Metric::XmasSKS
 * \brief Metric around a Schwarzschild-like BH 
 *        in (cartesian) Kerr-Schild coordinates
 */
class Gyoto::Metric::XmasSKS : public Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::XmasSKS>;
  


  // Data : 
  // -----
 protected:
  // INPUT FROM XML
  double Mdm_over_Mbh_ ;  /// DM/BH mass ratio
  double Rdm_ ;           /// DM radius 
  double gamma_;          /// DM slope




  /// Numerical tuning parameter
  /**
   * Small values yield more accurate integration at the expanse of
   * computing time.
   */
  double difftol_;
  double rsink_;  ///< numerical horizon
  double drhor_;  ///< horizon security
  bool   generic_integrator_; ///< which integrator to use
  


  // Constructors - Destructor
  // -------------------------
 public: 
  GYOTO_OBJECT;
  XmasSKS(); ///< Default constructor
  virtual XmasSKS * clone () const ;




  // Accessors
  // ---------
 public:
//  void spin(const double spin); ///< Set spin
//  double spin() const ; ///< Returns spin

  void Mdm_over_Mbh(const double Mdm_over_Mbh); ///< Set x
  double Mdm_over_Mbh() const ;                ///< Returns x

  void Rdm(const double Rdm);                  ///< Set x
  double Rdm() const ;  

  void gamma(const double gamma);                  ///< Set x
  double gamma() const ;  

  double difftol() const; ///< Get difftol_
  void difftol(double t); ///< Set difftol_

  void horizonSecurity(double drhor);
  double horizonSecurity() const;

  void genericIntegrator(bool);
  bool genericIntegrator() const ;

  void gmunu(double g[4][4], const double * pos) const ;
  double gmunu(const double * const x, int mu, int nu) const ;

  void gmunu_up(double gup[4][4], const double * pos) const ;
  double gmunu_up(const double * const x, int mu, int nu) const ;
 
 
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], const double pos[4]) const ;
  
  virtual void circularVelocity(double const pos[4], double vel [4],
        double dir=1.) const ;

  virtual int isStopCondition(double const * const coord) const;

};

#endif

