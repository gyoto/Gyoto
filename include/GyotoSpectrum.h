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

#include <GyotoRegister.h>

namespace Gyoto{
  class FactoryMessenger;
  namespace Spectrum {
    class Generic;
#if defined GYOTO_USE_XERCES
    typedef Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
      Subcontractor_t(Gyoto::FactoryMessenger* fmp);
    void Register(std::string, Gyoto::Spectrum::Subcontractor_t*);
    Gyoto::Spectrum::Subcontractor_t* getSubcontractor(std::string);
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
   * \param opacity, such that opacity*ds=optical thickness.
   * \param ds, in geometrical units
   */
  virtual double operator()(double nu, double opacity, double ds) const;
          ///< I_nu in optically thin regime.

  virtual double integrate(double nu1, double nu2) ;
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
  virtual void setParameter(std::string name, std::string content) ;
  ///< To be called by fillElement()

#endif
};


#endif
