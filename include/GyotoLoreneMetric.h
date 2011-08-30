/**
 *  \file GyotoLoreneMetric.h
 *  \brief Lorene computed metric
 *
 */

/*
    Copyright 2011 Frederic Vincent

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

#ifndef __GyotoLoreneMetric_H_
#define __GyotoLoreneMetric_H_ 

namespace Gyoto {
  class LoreneMetric;
}

#include <GyotoMetric.h>

/**
 * \class Gyoto::LoreneMetric
 * \brief Gyoto::Metric for Lorene computed metric
 */
class Gyoto::LoreneMetric : public Metric {
  friend class Gyoto::SmartPointer<Gyoto::LoreneMetric>;
  
  // Data : 
  // -----
 protected:
  
  // Constructors - Destructor
  // -------------------------
 public: 
  LoreneMetric(); ///< Default constructor

  // Default is _not_ fine
  LoreneMetric(const LoreneMetric& ) ;                ///< Copy constructor
  
  
  virtual ~LoreneMetric() ;                        ///< Destructor
  
  
  // Mutators / assignment
  // ---------------------
 public:

  // Accessors
  // ---------
 public:
  
  /** Value of metric coefficient $g_{\alpha\beta}$ at point $(x_{1},x_{2},x_{3})$
      in Boyer-Lindquist coordinates
   */
  virtual double gmunu(const double * const x,
		       const int alpha, const int beta,
		       const int sys) const ;

 
  /*
   it's necessary to define christoffel even if it's not used. LoreneMetric derives from Metric where christoffel is virtual pure. If the function is not defined in LoreneMetric,  it's considered virtual pure here too. Then LoreneMetric is considered an abstract class, and it's forbidden to declare any object of type LoreneMetric....
   See Delannoy C++ p.317-318
   NB : and it's not necessary to declare "virtual" a function in a derived class if it has been declared "virtual" in the basis class.
  */
  double christoffel(const double[8],
		     const int, const int, const int, 
		     const int) const;
  
  virtual double ScalarProd(const double* pos,
		    const double* u1, const double* u2,
		    const int sys) const ;

  friend std::ostream& operator<<(std::ostream& , const LoreneMetric& ) ;
  std::ostream& print(std::ostream&) const ;

 public:

 protected:
  /* RK4 : y=[r,theta,phi,t,pr,ptheta], cst=[a,E,L,Q],dy/dtau=F(y,cst), h=proper time step. For Kerr geodesic computation.
   */
  int myrk4_BL(const double y[6], const double* cst, double h, double* res) const;
  int myrk4_BL_adaptive(const double coor[8], const double* cst, const double MassPart, double lastnorm, double normref, int noprogress, double* coor1, double h0, double& h1, int onestep) const;//noprogress is 1 if integration makes no progress

  /** F function such as dy/dtau=F(y,cst)
   */
  int diff_BL(const double* y, const double* cst, double* res) const ;
  
};


#endif
