/*
    Copyright 2012 Thibaut Paumard

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

#include "GyotoFunctors.h"
#include "GyotoError.h"
#include <cmath>
#include <cfloat>
#include <iostream>

#define dr_tol 1e-5
#define f_tol 1e-9
#define sign(a) (a>=0?1:-1)

Gyoto::Functor::Double_constDoubleArray::~Double_constDoubleArray() {}

Gyoto::Functor::Double_Double_const::~Double_Double_const() {}

double Gyoto::Functor::Double_Double_const::ridders(double r1, double r2) const
{
  // Ridders' root-finding method applied to operator()
  double val1 = operator()(r1);
  double val2 = operator()(r2);
  double r3, r4, val3, val4;

  do {
    r3 = 0.5*(r1+r2);
    val3 = operator()(r3);

    r4 = r3+(r3-r1) * sign(val1-val2)*val3/sqrt(val3*val3-val1*val2);
    val4 = operator()(r4);

    if (sign(val3)!=sign(val4)) {
      // use r3 and r4 for next iteration
      r1=r3; val1=val3;
      r2=r4; val2=val4;
    } else if (sign(val1) != sign(val4)) {
      // use r1 and r4 for next iteration
      r2=r4; val2=val4;
    } else {
      // use r4 and r2 for next iteration
      r1=r4; val1=val4;
    }

  } while (fabs(r1-r2) >= dr_tol && fabs(val4) > f_tol) ;

  return r4;
}

double Gyoto::Functor::Double_Double_const::secant(double r1, double r2)
{
  // Secant root-finding method applied to operator()
  double val1 = operator()(r1);
  double val2 = operator()(r2);
  double r3, val3=DBL_MAX;
  unsigned int maxiter=20, iter;

  for (iter = 0; iter < maxiter; ++iter) {
    r3=r1-val1*(r2-r1)/(val2-val1);
    val3 = operator()(r3);
    
    if (fabs(r2-r3) <= fabs(r1-r3)) {
      r1=r3; val1=val3;
    } else {
      r2=r3; val2=val3;
    }
    
    if ( (fabs(r1-r2) < dr_tol) || (fabs(val3) < f_tol) ) break;
    if (val2==val1) {
      status=1;
      return r3;
    }
  }
  
  if (iter==maxiter) {
    status=2;
    return r3;
  }

  status=0;
  return r3;
}
