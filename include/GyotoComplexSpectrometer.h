/**
 * \file GyotoComplexSpectrometer.h
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


#ifndef __GyotoComplexSpectrometer_H_ 
#define __GyotoComplexSpectrometer_H_ 

#include <GyotoSpectrometer.h>

namespace Gyoto{
  namespace Spectrometer {
    class Complex;
  }
}

/**
 * \class Gyoto::Spectrometer::Complex
 * \brief Complex astronomical object
 *
 *  A Gyoto::Spectrometer::Generic whic contain several
 *  Gyoto::Spectrometer::Generic instances. It is essentially a
 *  SmartPointer<Spectrometer::Generic> array, which some methods
 *  arround. Indeed, the operator[](size_t i) method is implemented to retrived
 *  the i-th element.
 *
 * In an XML description, the
 *  &lt;Spectrometer&gt; section must be unique, its kind is
 *  "Complex". Each sub-astrobj then appears as a
 *  &lt;SubSpectrometer&gt; subsection:
\code
  <Spectrometer kind = "Complex">
    <SubSpectrometer kind = "ThinInfiniteDiskBL"/>
    <SubSpectrometer kind = "Star">
      <Radius> 2. </Radius>
      <Velocity> 0. 0. 0.037037 </Velocity>
      <Position> 600. 9. 1.5707999999999999741 0 </Position>
      <Spectrum kind="PowerLaw">
	<Exponent> 0 </Exponent>
	<Constant> 0.001 </Constant>
      </Spectrum>
      <Opacity kind="PowerLaw">
	<Exponent> 0 </Exponent>
	<Constant> 0.01 </Constant>
      </Opacity>
      <OpticallyThin/>
    </SubSpectrometer>
  </Spectrometer>
\endcode
 *
 */
class Gyoto::Spectrometer::Complex
: public Gyoto::Spectrometer::Generic,
  public Gyoto::Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer::Complex>;
  
  // Data : 
  // -----
 protected:

  /**
   * Number of objects
   */
  size_t cardinal_;

  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> * elements_;


 public:
  Complex(); ///< Default constructor.
  Complex(const Complex& ) ; ///< Copy constructor.
  virtual Complex* clone() const; ///< "Virtual" copy constructor

  /**
   *  Frees every SmartPointer<Spectrometer::Generic> before freed the array itself.
   */
  virtual ~Complex() ; ///< Destructor

  // Mutators
  // --------
 public:
  /**
   * If the Spectrometer::Complex itself does not have a metric already
   * assigned, it takes it from the new element. Else, it sets the
   * metric in the new element to its own. This ensures that all
   * elements use the same metric (this heuristic is not entirely
   * fool-proof, it's safer to set the metric directly in the
   * Spectrometer::Complex).
   */
  void append(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> element);
  ///< Add element at the end of the array.
  void remove(size_t i);
  ///< Remove i-th element from the array.
  size_t getCardinal() const;
  ///< Get the number of elements in the array.

  virtual void tell(Gyoto::Hook::Teller *msg); ///< Called by Teller

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
  virtual void setParameters(FactoryMessenger *fmp);
#endif
  
  // Outputs
  // -------
 public:


  /**
   * This should work as expected:
\code
  SmartPointer<Spectrometer::Complex> cplx;
  SmartPointer<Spectrometer::TypeA> objA;
  SmartPointer<Spectrometer::TypeB> objB;
  cplx -> append(objA);
  cplx[0] = objB;
\endcode
   */
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> operator[](size_t i) ;
  ///< Retrieve i-th element.
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> const operator[](size_t i) const;
  ///< Retrieve a const version of the i-th element.
  
  static SpectroKind_t const Kind;


 protected:
};

#endif
