/**
 * \file GyotoRotStar3_1.h
 * \brief Numerical metric around a rotating star in 3+1 formalism
 * 
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


#ifndef __GyotoRotStar3_1_H_
#define __GyotoRotStar3_1_H_ 

#include <iostream>
#include <fstream>

class Star_rot;

namespace Gyoto {
  namespace Metric { class RotStar3_1; }
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>
#include <GyotoSmartPointer.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

/**
 * \class Gyoto::Metric::RotStar3_1
 * \brief Numerical metric around a rotating star in 3+1 formalism
 */
class Gyoto::Metric::RotStar3_1 : public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::RotStar3_1>;

 private:
  char* filename_;
  Star_rot * star_;
  int integ_kind_;//1 if RotStar3_1::myrk4, 0 if Metric::myrk4
 
 public:

  RotStar3_1(const char * lorene_res, const int integ_kind); ///< Constructor
  virtual ~RotStar3_1() ;        ///< Destructor
  virtual RotStar3_1* clone() const ;
           ///< Cloner (uses RotStar3_1(file, integ_kind))

  char const * const getFileName() const;

  int getIntegKind() const ;

  using Metric::Generic::myrk4;
  int myrk4(const double coord[6], double h, double res[6]) const;

  //NB: there is no myrk4(const double coord[8], double h, double res[8]), which makes no pb because this same function is defined in class Metric where it is virtual but not pure. This myrk4(const double coord[6], double h, double res[6]) is a purely internal function, only called by RotStar3_1::myrk4_adaptive(const double coor[6],...)
  
  int myrk4_adaptive(Gyoto::Worldline* line, const double coord[8], double lastnorm, double normref, double coordnew[8], double h0, double& h1) const;

  int myrk4_adaptive(const double coor[6], double lastnorm, double normref, double coornew[6], double cst[2], double& tdot_used, double h0, double& h1, double& hused) const;
  /** F function such as dy/dtau=F(y,cst)
   */
  int diff(const double coord[8], double res[8]) const ;
  int diff(const double y[6], double res[6], int) const ;
  // NB: last int is just to distinguish with the other diff ; it doesn't seem to be possible to redefine diff with [8]->[6]

  void Normalize4v(const double coordin[8], double coordout[8], const double cst[2], double& tdot_used) const;

  double gmunu(const double * x, int mu, int nu) const ;

  double christoffel(const double coord[8], const int alpha, const int mu, 
		     const int nu) const ;

  double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const ;

#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp); ///< called from Factory
  static Metric::Subcontractor_t Subcontractor;
  static void Init();
#endif

};

#endif
