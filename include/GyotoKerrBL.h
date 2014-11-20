/**
 *  \file GyotoKerrBL.h
 *  \brief KerrBL metric
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

#ifndef __GyotoKerrBL_H_
#define __GyotoKerrBL_H_ 

namespace Gyoto {
  namespace Metric { class KerrBL; }
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif


/// Default value for difftol_
#define GYOTO_KERRBL_DEFAULT_DIFFTOL 1e-2

/**
 * \class Gyoto::Metric::KerrBL
 * \brief Metric around a Kerr black-hole in Boyer-Lindquist coordinates
 */
class Gyoto::Metric::KerrBL : public Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::KerrBL>;
  
  // Data : 
  // -----
 protected:
  double spin_ ;  ///< Angular momentum parameter
  double a2_ ; ///< spin_*spin_
  double a3_ ; ///< a2_*spin_
  double a4_ ; ///< a2_*a2_

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
  KerrBL(); ///< Default constructor
  virtual KerrBL * clone () const ;

  // Accessors
  // ---------
 public:
  void spin(const double spin); ///< Set spin
  double spin() const ; ///< Returns spin

  double difftol() const; ///< Get difftol_
  void difftol(double t); ///< Set difftol_

  void horizonSecurity(double drhor);
  double horizonSecurity() const;
  void genericIntegrator(bool);
  bool genericIntegrator() const ;

  virtual double getRms() const; 

  virtual double getRmb() const; 

  virtual double getSpecificAngularMomentum(double rr) const;
  
  virtual double getPotential(double pos[4], double l_cst) const;

  void gmunu(double g[4][4], const double * pos) const ;
  double gmunu(const double * const x, int mu, int nu) const ;

  /** 
   * \brief g<SUP>&mu;,&nu;</SUP>
   */
  void gmunu_up(double gup[4][4], const double * pos) const ;
  double gmunu_up(const double * const x, int mu, int nu) const ;
 
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], const double * pos) const ;
  
  double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const ;

  void nullifyCoord(double coord[8], double & tdot2) const;
  void nullifyCoord(double coord[8]) const;

  //  friend std::ostream& operator<<(std::ostream& , const KerrBL& ) ;
  //  std::ostream& print(std::ostream&) const ;
  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

 public:
  virtual void MakeCoord(const double coordin[8], const double cst[5], double coordout[8]) const ;
  ///< Inverse function of MakeMomentumAndCst

   ///< Computes pr, ptheta, E and L from rdot, thetadot, phidot, tdot
  void MakeMomentum(const double coordin[8], const double cst[5], double coordout[8]) const;
  ///< Transforms from Boyer-Lindquist coordinates [t,r,th,phi,tdot,rdot,thdot,phidot] to [t,r,th,phi,pt,pr,pth,pphi] where pt,pr... are generalized momenta.

  virtual void setParameter(std::string, std::string, std::string);
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp); ///< called from Factory
#endif

 protected:

  // outside the API
  /* RK4 : y=[r,theta,phi,t,pr,ptheta], 
     cst=[a,E,L,Q,1/Q],dy/dtau=F(y,cst), h=proper time step. 
     For KerrBL geodesic computation.
   */
  int myrk4(Worldline * line, const double coordin[8], double h, double res[8]) const; //external-use RK4
  
 private:
  int myrk4(const double coor[8], const double cst[5], double h, double res[8]) const;///< Internal-use RK4 proxy
  int myrk4_adaptive(Gyoto::Worldline* line, const double coor[8], double lastnorm, double normref, double coor1[8], double h0, double& h1, double h1max) const; ///< Internal-use adaptive RK4 proxy
  /**
   * \brief Ensure conservation of the constants of motion
   *
   * Tweak thetadot if necessary.
   */
  int CheckCons(const double coor_init[8], const double cst[5], double coor_fin[8]) const;

  /**
   * \brief Normalize 4-velocity
   *
   * To 0 or -1. Changes rdot to allow norm conservation.
   */
  void Normalize4v(double coord[8], const double part_mass) const;

  /** F function such as dy/dtau=F(y,cst)
   */
  using Metric::Generic::diff;
  /** 
   * \brief Used in RK4 proxies.
   */
  virtual int diff(const double y[8], const double cst[5], 
		   double res[8]) const ;
  /** Integrator. Computes the evolution of y (initcond=y(0)).
   */
  virtual void computeCst(const double coord[8], double cst[5]) const;
 public:
  void setParticleProperties(Worldline* line, const double* coord) const;
  virtual int isStopCondition(double const * const coord) const;
  
  void observerTetrad(std::string const obskind,
		      double const pos[4], double fourvel[4],
		      double screen1[4], double screen2[4], 
		      double screen3[4]) const;
};

#endif

