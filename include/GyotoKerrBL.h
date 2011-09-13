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
  class KerrBL;
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif


/**
 * \class Gyoto::KerrBL
 * \brief Gyoto::Metric for KerrBL black hole
 */
class Gyoto::KerrBL : public Metric {
  friend class Gyoto::SmartPointer<Gyoto::KerrBL>;
  
  // Data : 
  // -----
 protected:
  double spin_ ;  ///< Angular momentum parameter
  
  // Constructors - Destructor
  // -------------------------
 public: 
  KerrBL(); ///< Default constructor
  KerrBL(double spin, double mass) ; ///< Constructor with spin and mass specification

  // Default is _not_ fine
  //KerrBL(const KerrBL& ) ;                ///< Copy constructor
  
  
  virtual ~KerrBL() ;                        ///< Destructor
  
  
  // Mutators / assignment
  // ---------------------
 public:
  // default operator= is fine
  void setSpin(const double spin); ///< Set spin
  virtual KerrBL * clone () const ;


  // Accessors
  // ---------
 public:
  double getSpin() const ; ///< Returns spin

  double getRms() const; ///< Returns prograde marginally stable orbit

  double getRmb() const; ///< Returns prograde marginally bound orbit
  
  /** Value of metric coefficient $g_{\alpha\beta}$ at point $(x_{1},x_{2},x_{3})$
      in Boyer-Lindquist coordinates
   */
  double gmunu(const double * const x,
		       int alpha, int beta) const ;

 
  /*
   it's necessary to define christoffel even if it's not used. KerrBL derives from Metric where christoffel is virtual pure. If the function is not defined in KerrBL,  it's considered virtual pure here too. Then KerrBL is considered an abstract class, and it's forbidden to declare any object of type KerrBL....
   See Delannoy C++ p.317-318
   NB : and it's not necessary to declare "virtual" a function in a derived class if it has been declared "virtual" in the basis class.
  */
  double christoffel(const double[8],
		     const int, const int, const int) const;
  
  double ScalarProd(const double* pos,
		    const double* u1, const double* u2) const ;

  void nullifyCoord(double coord[8], double & tdot2) const;
  void nullifyCoord(double coord[8]) const;


  //  friend std::ostream& operator<<(std::ostream& , const KerrBL& ) ;
  //  std::ostream& print(std::ostream&) const ;
  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

 public:
  void MakeCoord(const double coordin[8], const double cst[5], double coordout[8]) const ;
  ///< Inverse function of MakeMomentumAndCst

   ///< Computes pr, ptheta, E and L from rdot, thetadot, phidot, tdot
  void MakeMomentum(const double coordin[8], const double cst[5], double coordout[8]) const;
  ///< Transforms from Boyer-Lindquist coordinates [t,r,th,phi,tdot,rdot,thdot,phidot] to [t,r,th,phi,pt,pr,pth,pphi] where pt,pr... are generalized momenta.
 

#ifdef GYOTO_USE_XERCES
  virtual void fillElement(factoryMessenger *fmp); /// < called from Factory
  static Metric::Subcontractor_t Subcontractor;
  static void Init();
#endif

 protected:

  // outside the API
  /* RK4 : y=[r,theta,phi,t,pr,ptheta], cst=[a,E,L,Q,1/Q],dy/dtau=F(y,cst), h=proper time step. For KerrBL geodesic computation.
   */
  int myrk4(Worldline * line, const double coordin[8], double h, double res[8]) const; //external-use RK4
  int myrk4(const double coor[8], const double cst[5], double h, double res[8]) const;//internal-use RK4
  int myrk4_adaptive(Gyoto::Worldline* line, const double coor[8], double lastnorm, double normref, double coor1[8], double h0, double& h1) const;
  int CheckCons(const double coor_init[8], const double cst[5], double coor_fin[8]) const;
  void Normalize4v(double coord[8], const double part_mass) const;

  /** F function such as dy/dtau=F(y,cst)
   */
  using Metric::diff;
  int diff(const double y[8], const double cst[5], double res[8]) const ;
  /** Integrator. Computes the evolution of y (initcond=y(0)).
   */
  void computeCst(const double coord[8], double cst[5]) const;
  void setParticleProperties(Worldline* line, const double* coord) const;
  
};

#endif

