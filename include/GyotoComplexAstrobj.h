/**
 * \file GyotoComplexAstrobj.h
 * \brief Combine astronomical objects
 *
 *  An astrobj made of several astrobjs
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


#ifndef __GyotoComplexAstrobj_H_ 
#define __GyotoComplexAstrobj_H_ 

#include <GyotoAstrobj.h>

namespace Gyoto{
  namespace Astrobj {
    class Complex;
  }
}

/**
 * \class Gyoto::Astrobj::Complex
 * \brief Complex astronomical object
 *
 *  A Gyoto::Astrobj::Generic whic contain several
 *  Gyoto::Astrobj::Generic instances.
 *
 *  When implementing a new object, you must:
 *    - make sure the object can be loaded from XML by providing a
 *      subcontractor;
 *    - make sure this subcontractor is registerred in the initialization
 *      routine of your plug-in;
 *    - make sure Impact() works (see below).
 *
 *  In addition, you should make sure that your object plays nicely in
 *  the Yorick plug-in, which means:
 *    - implement the copy constructor and the clone() method;
 *    - implement the fillElement method, used for printing and saving to
 *      XML.
 *
 *  There are basically two ways of making Impact() work: either by
 *  providing your own Impact() function, or by implementing a bunch
 *  of lower level, simpler functions which are used by the generic
 *  Astrobj::Impact(). Those lower functions are not pure virtual, but
 *  the default throws a runtime error:
 *    - operator()() yields a distance or potential defining the interior
 *      of the object;
 *    - getVelocity() yields the velocity field of the fluid ;
 *    - processHitQuantities() fills the Spectrum etc. quantities in the
 *      data parameter of Impact().
 *
 *  processHitQuantities() itself is like Impact() in that you have
 *  two choices: either reimplement it or implement a second lot of
 *  small, low-level functions:
 *    - emission();
 *    - integrateEmission();
 *    - transmission().
 *
 */
class Gyoto::Astrobj::Complex : public Gyoto::Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Complex>;
  
  // Data : 
  // -----
 protected:

  /**
   * Number of objects
   */
  size_t cardinal_;

  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> * elements_;

  double step_max_;

 public:
  /**
   *  kind_ =  "Default", rmax_ = 0., rmax_set_ = 0.
   */
  Complex(); ///< Default constructor.

  Complex(const Complex& ) ; ///< Copy constructor.
  virtual Complex* clone() const; ///< "Virtual" copy constructor
  
  virtual ~Complex() ; ///< Destructor: does nothing.

  // Mutators
  // --------
 public:
  void append(Gyoto::SmartPointer<Gyoto::Astrobj::Generic> element);
  void remove(size_t i);
  size_t getCardinal() const;
  void setMetric(SmartPointer<Metric::Generic> gg);

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
                                             /// < called from Factory
  static Astrobj::Subcontractor_t Subcontractor;
  static void Init();

#endif
  
  // Outputs
  // -------
 public:
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL)  ;
  ///< does a photon at these coordinates impact the object?


  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> operator[](size_t i) ;
  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> const operator[](size_t i) const ;
  
 protected:
};

#endif
