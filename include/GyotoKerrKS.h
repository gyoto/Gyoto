/**
 *  \file GyotoKerrKS.h
 *  \brief KerrKS metric
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
 */
class Gyoto::Metric::KerrKS : public Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::KerrKS>;
  
  // Data : 
  // -----

 protected:
  double spin_ ;  ///< Angular momentum parameter

  
  // Constructors - Destructor
  // -------------------------
 public: 
  KerrKS(); ///< Default constructor
  KerrKS(double spin, double mass) ; ///< Constructor with spin and mass specification

  // Default is fine
  // KerrKS(const KerrKS& ) ;       
  virtual KerrKS* clone () const;         ///< Copy constructor
  
  virtual ~KerrKS() ;                        ///< Destructor
  
  
  // Mutators / assignment
  // ---------------------
 public:
  // default operator= is fine
  void setSpin(const double spin); ///< Set spin

  // Accessors
  // ---------
 public:
  double getSpin() const ; ///< Returns spin
  
  /** Value of metric coefficient $g_{\alpha\beta}$ at point $(x_{1},x_{2},x_{3})$
      in Kerr-Schild coordinates
   */
  double gmunu(const double * x,
		       int alpha, int beta) const ;

 
  /*
   it's necessary to define christoffel even if it's not used. KerrKS derives from Metric where christoffel is virtual pure. If the function is not defined in KerrKS,  it's considered virtual pure here too. Then KerrKS is considered an abstract class, and it's forbidden to declare any object of type KerrKS....
   See Delannoy C++ p.317-318
   NB : and it's not necessary to declare "virtual" a function in a derived class if it has been declared "virtual" in the basis class.
  */
  double christoffel(const double[8],
		     const int, const int, const int) const;
  

  void nullifyCoord(double coord[8], double &tdot2) const;
  void nullifyCoord(double coord[8]) const;
  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

  //  friend std::ostream& operator<<(std::ostream& , const KerrKS& ) ;
  //  std::ostream& print(std::ostream&) const ;
  virtual void setParameter(std::string, std::string, std::string);
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp); ///< called from Factory
#endif

 public:

  void MakeCst(const double* coord, double* cst) const;
  ///< In Kerr-Schild coordinates [T,x,y,z,Tdot,xdot,ydot,zdot], computes the four constants of the movement : particule mass, energy, angular momentum and Carter's constant.
 protected:

  // outside the API
  /* RK4 : y=[r,theta,phi,t,pr,ptheta], cst=[a,E,L,Q],dy/dtau=F(y,cst), h=proper time step. For KerrKS geodesic computation.
   */
  int myrk4(Worldline * line, const double coord[8], double h, double res[8]) const;//NB non adaptive integration doesn't work for KS ; this function is not implemented
  int myrk4(const double * coord, const double* cst , double h, double* res) const;
  int myrk4_adaptive(Gyoto::Worldline* line, const double * coord, double lastnorm, double normref, double* coord1, double h0, double& h1) const;

  /** F function such as dy/dtau=F(y,cst)
   */
  int diff(const double* coord, const double* cst, double* res) const;
  int diff(const double y[8], double res[8]) const ;
  virtual int isStopCondition(double const * const coord) const;
  void setParticleProperties(Worldline* line, const double* coord) const;

};

#endif
