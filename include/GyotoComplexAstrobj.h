/**
 * \file GyotoComplexAstrobj.h
 * \brief Combine astronomical objects
 *
 *  An astrobj made of several astrobjs
 */

/*
    Copyright 2011, 2013-2014, 2016 Thibaut Paumard

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
 * A Gyoto::Astrobj::Generic whic contain several
 * Gyoto::Astrobj::Generic instances. It is essentially a
 * SmartPointer<Astrobj::Generic> array, which some methods
 * around. Indeed, the operator[](size_t i) method is implemented to
 * retrieve the i-th element.
 *
 * In an XML description, the &lt;Astrobj&gt; section must be unique,
 * its kind is "Complex". Each sub-astrobj then appears as a
 * &lt;SubAstrobj&gt; subsection:
 * \code
 *   <Astrobj kind = "Complex">
 *     <SubAstrobj kind = "ThinInfiniteDiskBL"/>
 *     <SubAstrobj kind = "Star">
 *       <Radius> 2. </Radius>
 *       <Velocity> 0. 0. 0.037037 </Velocity>
 *       <Position> 600. 9. 1.5707999999999999741 0 </Position>
 *       <Spectrum kind="PowerLaw">
 * 	<Exponent> 0 </Exponent>
 * 	<Constant> 0.001 </Constant>
 *       </Spectrum>
 *       <Opacity kind="PowerLaw">
 * 	<Exponent> 0 </Exponent>
 * 	<Constant> 0.01 </Constant>
 *       </Opacity>
 *       <OpticallyThin/>
 *     </SubAstrobj>
 *   </Astrobj>
 * \endcode
 *
 */
class Gyoto::Astrobj::Complex : public Gyoto::Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Complex>;
  
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
  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> * elements_;

  /**
   * Currently not settable, always equal to GYOTO_DEFAULT_DELTA
   */
  double step_max_; ///< Maximum &delta; step inside the Astrobj

 public:
  GYOTO_OBJECT_THREAD_SAFETY;
  Complex(); ///< Default constructor.
  Complex(const Complex& ) ; ///< Copy constructor.
  virtual Complex* clone() const; ///< "Virtual" copy constructor

  virtual double deltaMax(double coord[8]);

  /**
   * rMax loops over the elementary astrobjs rmax_ and returns the biggest
   */
  virtual double rMax();

  /**
   *  Frees every SmartPointer<Astrobj::Generic> before freed the array itself.
   */
  virtual ~Complex() ; ///< Destructor

  // Mutators
  // --------
 public:
  /**
   * If the Astrobj::Complex itself does not have a metric already
   * assigned, it takes it from the new element. Else, it sets the
   * metric in the new element to its own. This ensures that all
   * elements use the same metric (this heuristic is not entirely
   * fool-proof, it's safer to set the metric directly in the
   * Astrobj::Complex).
   */
  void append(Gyoto::SmartPointer<Gyoto::Astrobj::Generic> element);
  ///< Add element at the end of the array.
  void remove(size_t i);
  ///< Remove i-th element from the array.
  size_t getCardinal() const;
  ///< Get the number of elements in the array.
  using Generic::metric;
  void metric(SmartPointer<Metric::Generic> gg);
  ///< Set metric in each element.

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
  virtual void setParameters(FactoryMessenger *fmp);
#endif
  
  // Outputs
  // -------
 public:

  /**
   * Astrobj::Complex::Impact(Gyoto::Photon* ph, size_t index,
   * Astrobj::Properties *data) calls the specific implementation of
   * Astrobj::Generic::Impact() for each of its
   * elements twice: the first time, data is set to NULL so that
   * Astrobj::Complex::Impact() only knows whether each object is hit
   * by the Photon. If no object is hit, return. If a single object is
   * hit, call Impact() again only for this object, passing data this
   * time. If several objects are hit, the Photon's trajectory is
   * refined so that the step is at most step_max_ and the Impact()
   * methods for each of the hit objects are called again for each
   * step, passing data. It is therefore important that the
   * transmission of the Photon is not touched by Impact() when
   * data==NULL.
   * 
   */
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL)  ;
  ///< Call Impact() for each of the elements.


  /**
   * This should work as expected:
   * \code
   *   SmartPointer<Astrobj::Complex> cplx;
   *   SmartPointer<Astrobj::TypeA> objA;
   *   SmartPointer<Astrobj::TypeB> objB;
   *   cplx -> append(objA);
   *   cplx[0] = objB;
   * \endcode
   */
  Gyoto::SmartPointer<Gyoto::Astrobj::Generic>& operator[](size_t i) ;
  ///< Retrieve i-th element.
  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> const&operator[](size_t i) const;
  ///< Retrieve a const version of the i-th element.
  
};

#endif
