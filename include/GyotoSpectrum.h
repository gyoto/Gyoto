/**
 * \file GyotoSpectrum.h
 * \brief Spectrum of a simple object (e.g. Star)
 *
 *  Light emitted by an astronomical object
 */

/*
    Copyright 2011-2012 Thibaut Paumard

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

#ifndef __GyotoSpectrum_H_ 
#define __GyotoSpectrum_H_ 

#include "GyotoRegister.h"

namespace Gyoto{
  class FactoryMessenger;
  namespace Spectrum {
    class Generic;
#if defined GYOTO_USE_XERCES
    /**
     * This is a more specific version of the
     * SmartPointee::Subcontractor_t type. A Spectrum::Subcontrator_t
     * is called by the Gyoto::Factory to build an instance of the
     * kind of spectrum specified in an XML file (see
     * Register()). The Factory and Subcontractor_t function
     * communicate through a Gyoto::FactoryMessenger. A template is
     * provided so that you may not have to code anything.
     */
    typedef Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
      Subcontractor_t(Gyoto::FactoryMessenger* fmp);

    /**
     * Instead of reimplementing the wheel, your subcontractor can simply be
     * Gyoto::Astrobj::Subcontractor<MyKind>
     */
    template<typename T> SmartPointer<Spectrum::Generic> Subcontractor
      (FactoryMessenger* fmp) {
      SmartPointer<T> sp = new T();
      sp -> setParameters(fmp);
      return sp;
    }
    ///< A template for Subcontractor_t functions

    void Register(std::string, Gyoto::Spectrum::Subcontractor_t*);
    Gyoto::Spectrum::Subcontractor_t*
      getSubcontractor(std::string, int errmode=0);
    extern Register::Entry* Register_;
    void initRegister();
#endif
  }
}

#include <GyotoSmartPointer.h>
#include <string>
/**
 * \class Gyoto::Spectrum::Generic
 * \brief Spectrum emitted by an Astrobj
 *
 *  Light emitted by e.g. a Star
 *
 */
class Gyoto::Spectrum::Generic : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::Generic>;
 protected:
  std::string kind_; ///< e.g. constants, blackbody...

 public:
  Generic(const std::string kind); ///< Set kind in constructor
  //  Spectrum::Generic(const Spectrum::Generic &); ///< Copy constructor. Default is fine.
  virtual Generic * clone() const; ///< Cloner

  virtual ~Generic() ; ///< Destructor: does nothing.

  const std::string getKind() const; ///< Get spectrum kind

  virtual double operator()(double nu) const =0;
          ///< I_nu = mySpectrum(nu), nu in Hz. Assumes optically thick regime.
  /**
   * Generic implementation assumes emissivity = opacity.
   *
   * \param nu frequency in Hz
   * \param opacity such that opacity*ds=optical thickness.
   * \param ds in geometrical units
   */
  virtual double operator()(double nu, double opacity, double ds) const;
          ///< I_nu in optically thin regime.

  /**
   * \brief Integrate optically thick I_nu
   *
   * See operator()(double nu) const
   *
   * \param nu1, nu2 boundaries for the integration
   * \result I, the integral of I_nu between nu1 and nu2
   */
  virtual double integrate(double nu1, double nu2) ;

  /**
   * \brief Integrate optically thin I_nu
   *
   * See operator()(double nu, double opacity, double ds) const
   *
   * \param nu1, nu2 boundaries for the integration
   * \param opacity the frequency-dependent opacity law given as a
   *        pointer to a Gyoto::Spectrum::Generic sub-class instance
   * \param ds the element length for spatial integration
   * \result I, the integral of I_nu between nu1 and nu2
   */
  virtual double integrate(double nu1, double nu2,
			   const Spectrum::Generic * opacity, double ds) ;

#ifdef GYOTO_USE_XERCES
  /**
   * Spectrum implementations should impement fillElement to save their
   * parameters to XML and call the generic implementation to save
   * generic parts.
   */

  virtual void fillElement(FactoryMessenger *fmp) const ;
                                             ///< called from Factory
  virtual void setParameter(std::string name,
			    std::string content,
			    std::string unit) ;
  ///< To be called by fillElement()

  /**
   * The Subcontractor_t function for each Spectrum kind should look
   * somewhat like this:
   * \code
   * SmartPointer<Astrobj::Generic>
   * Gyoto::Spectrum::MyKind::Subcontractor(FactoryMessenger* fmp) {
   *   SmartPointer<MyKind> ao = new MyKind();
   *   ao -> setParameters(fmp);
   *   return ao;
   * }
   * \endcode
   *
   * Each object kind should implement setParameter(string name,
   * string content) to interpret the individual XML
   * elements. setParameters() can be overloaded in case the specific
   * Astrobj class needs low level access to the FactoryMessenger. See
   * Astrobj::UniformSphere::setParameters().
   */
  virtual void setParameters(FactoryMessenger *fmp);
  ///< Main loop in Subcontractor_t function

#endif
};


#endif
