/**
 * \file GyotoStandardAstrobj.h
 * \brief Astronomical objects defined bya a potential/distance
 *
 *  Many geometrically thick objects can be defined by the value of a
 *  function of the 4 coordinates, and their emission can often be
 *  defined in terms of an emission law and of a transmission law.
 *
 *  This is a base class for this standard case which simplifies a lot
 *  writting new Astrobjs.
 */

/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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


#ifndef __GyotoStandardAstrobj_H_ 
#define __GyotoStandardAstrobj_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "GyotoAstrobj.h"
#include "GyotoFunctors.h"

namespace Gyoto{
  namespace Astrobj {
    class Standard;
  }
}

/**
 * \class Gyoto::Astrobj::Standard
 * \brief Astronomical objects defined bya a potential/distance
 *
 *  Many geometrically thick objects can be defined by the value of a
 *  function of the 4 coordinates, and their emission can often be
 *  defined in terms of an emission law and of a transmission law.
 *
 *  This is a base class for this standard case which simplifies a lot
 *  writting new Astrobjs.
 *
 *  It is either to implement a sub-class of Astrobj::Standard than a
 *  sub-class of Astrobj::Generic. In particular, there is no need to
 *  implement the Generic::Impact() function. Instead, one needs to
 *  implement a few much simpler functions and most of the complex
 *  ray-tracing algorithms and heuristics is implemented in
 *  Standard::Impact(). It is recommended to read first the
 *  introduction in the Gyoto::Astrobj namespace documentation.
 *
 *  The geometrical shape of a Gyoto::Astrobj::Standard object is
 *  yielded by a function of the 4 position vector. This function is
 *  implemented as operator()(). The velocity field of the fluid is
 *  implemented in the getVelocity() method. The emission(),
 *  integrateEmission() and transmission() methods implement the
 *  radiative transfer primitives for this object. Finally, you may
 *  choose to reimplement processHitQuantities() and Impact(), but
 *  this should not be necessary (that is the all point of the
 *  Standard class).
 *
 * Like any other Astrobj::Generic sub-classes, an Astrobj::Standard
 * subclass should register an Astrobj::Subcontractor_t function using
 * the Astrobj::Register() function. See also \ref
 * writing_plugins_page .
 */
class Gyoto::Astrobj::Standard :
  public Gyoto::Astrobj::Generic,
  public Gyoto::Functor::Double_constDoubleArray
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Standard>;


  // Data : 
  // -----
 protected:
  double critical_value_; ///< See operator()(double const coord[4])
  double safety_value_; ///< See operator()(double const coord[4])
  double delta_inobj_; ///< Constant value of the integration step inside object, in units of the compact object's mass M

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

  /**
   *  kind_ =  "Default", rmax_ = 0., rmax_set_ = 0.
   */
  Standard(); ///< Default constructor.

  /**
   *  kind_ =  "Default", rmax_ = radmax, rmax_set_ = 1.
   */
  Standard(double radmax); ///< Set rmax in constructor.

  /**
   *  kind_ =  kind, rmax_ = 0., rmax_set_ = 0.
   */
  Standard(std::string kind); ///< Set kind in constructor.

  /**
   * Make a deep copy of an Astrobj::Standard instance
   */
  Standard(const Standard& ) ; ///< Copy constructor.

  virtual ~Standard() ; ///< Destructor: does nothing.

  // Accessors
  // ---------
 public:
  virtual void safetyValue(double val) ; ///< Set Standard::safety_value_
  virtual double safetyValue() const ; ///< Get Standard::safety_value_

  /**
   *  Get the constant integration step inside the astrobj
   *
   *  \return delta_inobj_ in geometrical units
   */
  double deltaInObj() const; ///< Get Generic::delta_inobj_
  void   deltaInObj(double val); ///< Set Generic::delta_inobj_

  // Outputs
  // -------
 public:
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL)  ;

  /**
   * \brief Function defining the object interior
   *
   * A potential, distance, or whatever function such that
   * operator()(double const coord[4]) < Standard::critical_value_ if
   * and only if coord is inside the object. This function is used by
   * the default implmenetation of Impact(). If Impact() is
   * overloaded, it is not necessary to overload operator()(double
   * coord[4]). The default implementation throws an error.
   */
  virtual double operator()(double const coord[4]) = 0;
  
  /**
   * \brief Fluid velocity field.
   *
   * Fill vel with the 4-vector velocity of the fluid at 4-position pos.
   *
   * \param[in] pos 4-position at which to compute velocity;
   * \param[out] vel 4-velocity at pos.
   */
  virtual void getVelocity(double const pos[4], double vel[4]) = 0 ;

  /**
   * \brief Maximum &delta; inside object
   *
   * Gives the requested integration step &delta;<SUB>t</SUB> (in
   * coordinate time t) between two neighbooring points along a
   * portion of geodesic inside an astrobj; the current implementation
   * only considers a constant delta, equal to Standard::deltaInobj()
   *
   * \param coord input coordinate at which &delta;<SUB>t</SUB> is given
   */
  virtual double giveDelta(double coord[8]);



};


#endif
