/**
 * \file GyotoCompexMetric.h
 * \brief Combine metrics
 * 
 */

/*
    Copyright 2020 Thibaut Paumard

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


#ifndef __GyotoComplexMetric_H_
#define __GyotoComplexMetric_H_

#include <GyotoMetric.h>
#include <GyotoWIP.h>

namespace Gyoto {
  namespace Metric { class Complex; }
}

/**
 * \class Gyoto::Metric::Complex
 * \brief Combine several metrics
 * 
 * Adds linearly the contribution of several metrics. All suv-metrics
 * must use a common coordinate system.
 *
 */

class Gyoto::Metric::Complex
  : public Gyoto::Metric::Generic, Gyoto::WIP
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::Complex>;

  // Data : 
  // -----
 protected:

  /**
   * \brief Number of objects
   */
  size_t cardinal_;

  /**
   * \brief Array of Astrobj::Generic
   */
  Gyoto::SmartPointer<Gyoto::Metric::Generic> * elements_;

 public:
  GYOTO_OBJECT_THREAD_SAFETY;
  Complex(); ///< Default constructor.
  Complex(const Complex& ) ; ///< Copy constructor.
  virtual Complex* clone() const; ///< "Virtual" copy constructor
  /**
   *  Frees every SmartPointer<Metric::Generic> before freed the array itself.
   */
  virtual ~Complex() ; ///< Destructor
  void append(Gyoto::SmartPointer<Gyoto::Metric::Generic> element);
  ///< Add element at the end of the array.
  void remove(size_t i);
  ///< Remove i-th element from the array.
  size_t getCardinal() const;
  ///< Get the number of elements in the array.
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
  virtual void setParameters(FactoryMessenger *fmp);
#endif
  /**
   * This should work as expected:
   * \code
   *   SmartPointer<Metric::Complex> cplx;
   *   SmartPointer<Metric::TypeA> objA;
   *   SmartPointer<Metric::TypeB> objB;
   *   cplx -> append(objA);
   *   cplx[0] = objB;
   * \endcode
   */
  Gyoto::SmartPointer<Gyoto::Metric::Generic>& operator[](size_t i) ;
  ///< Retrieve i-th element.
  Gyoto::SmartPointer<Gyoto::Metric::Generic> const&operator[](size_t i) const;
  ///< Retrieve a const version of the i-th element.

 public:
  using Generic::gmunu;
  void gmunu(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;
  void jacobian(double ARGOUT_ARRAY3[4][4][4], const double IN_ARRAY1[4]) const;
  int isStopCondition(double const coord[8]) const;

};

#endif
