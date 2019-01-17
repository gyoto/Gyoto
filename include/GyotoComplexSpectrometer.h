/**
 * \file GyotoComplexSpectrometer.h
 * \brief Combine spectrometer objects
 *
 *  A spectrometer made of several spectrometers
 */

/*
    Copyright 2013-2014, 2016, 2019 Thibaut Paumard

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
 * \brief Complex spectrometer object
 *
 *  A Gyoto::Spectrometer::Generic whic contain several
 *  Gyoto::Spectrometer::Generic instances. It is essentially a
 *  SmartPointer<Spectrometer::Generic> array, which some methods
 *  arround. Indeed, the operator[](size_t i) method is implemented to
 *  retrieve the i-th element.
 *
 * In an XML description, the &lt;Spectrometer&gt; section is unique,
 *  its kind is "Complex". Each sub-spectrometer then appears as a
 *  &lt;SubSpectrometer&gt; subsection. For instance, to compute 10
 *  channels ovr the K infrared band plus 10 channels in the high
 *  energy domain:
 * \code
 * <Spectrometer kind = "Complex">
 *   <SubSpectrometer kind = "wave" nsamples=10 unit="µm">
 *     2.0 2.4
 *   </SubSpectrometer>
 *   <SubSpectrometer kind = "freqlog" nsamples=10 unit="eV">
 *     14 22
 *   </SubSpectrometer>
 * </Spectrometer>
 * \endcode
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
   * \brief Number of subspectrometers
   *
   * Size of elements_.
   */
  size_t cardinal_;

  /**
   * \brief Actual array of SmartPointer<Spectrometer::Generic> objects
   */
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> * elements_;


 public:
  GYOTO_OBJECT_THREAD_SAFETY;
  Complex(); ///< Default constructor.
  Complex(const Complex& ) ; ///< Copy constructor.

  /**
   * \brief Clone an instance
   * 
   * Use this to get a deep copy of an instance;
   * \code
   * SmartPointer<Generic> myclone = orig->clone();
   * \endcode
   *
   * Most implementations will use the copy constructor:
   * \code
   * Generic* Uniform::clone() const { return new Uniform(*this); }
   * \endcode
   */
  virtual Complex* clone() const;

  /**
   *  Frees every SmartPointer<Spectrometer::Generic> before freed the
   *  array itself.
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

  virtual void tell(Gyoto::Hook::Teller *msg);

 public:
#ifdef GYOTO_USE_XERCES
  /**
   * \brief Fill in the XML entity
   *
   * Loops on elements_[i]->fillElement();
   */ 
  virtual void fillElement(FactoryMessenger *fmp) const ;

  /**
   * \brief Main loop in the (templated) subcontractor
   *
   * In the case of Spectrometer::Complex, the setParameter() API is
   * not sufficient: setParameters() needs access to the
   * FactoryMessenger to instantiate children for the
   * SubSpectrometers.
   */
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
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> & operator[](size_t i) ;
  ///< Retrieve i-th element.
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> const & operator[](size_t i) const;
  ///< Retrieve a const version of the i-th element.
  
  /**
   * \brief "Complex"
   *
   * Use this static member attribute to check whether a Spectrometer
   * object spectro is of kind Complex:
   * \code
   * if (spectro->kind() == Complex::Kind) ... ;
   * \endcode
   *
   */
  static kind_t const Kind;


 protected:
};

#endif
