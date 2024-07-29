/**
 * \file GyotoSpectrum.h
 * \brief Spectrum of a simple object (e.g. Star)
 *
 *  Light emitted by an astronomical object
 */

/*
    Copyright 2011-2016, 2022, 2024 Thibaut Paumard

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
#include "GyotoObject.h"

namespace Gyoto{
  namespace Register { class Entry; }
  class FactoryMessenger;

  /// Spectrum of a simple object (e.g. a Gyoto::Astrobj::Star)
  namespace Spectrum {
    class Generic;

    /// A function to build instances of a specific Spectrum::Generic sub-class
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
      Subcontractor_t(Gyoto::FactoryMessenger* fmp, std::vector<std::string> const &);

    /// Subcontractor template
    /**
     * Instead of reimplementing the wheel, your subcontractor can simply be
     * Gyoto::Spectrum::Subcontractor<MyKind>
     *
     * \tparam T Sub-class of Spectrum::Generic 
     */
    template<typename T> SmartPointer<Spectrum::Generic> Subcontractor
      (FactoryMessenger* fmp, std::vector<std::string> const & plugins) {
      SmartPointer<T> sp = new T();
      sp -> plugins(plugins) ;
#ifdef GYOTO_USE_XERCES
      if (fmp) sp -> setParameters(fmp);
#endif
      return sp;
    }

    /// Make a Spectrum kind known to the Factory
    /**
     * Register a new Spectrum::Generic sub-class so that the
     * Gyoto::Factory knows it.
     *
     * \param kind The kind name which identifies this object type in
     * an XML file, as in &lt;Spectrum kind="name"&gt;
     *
     * \param scp A pointer to the subcontractor, which will
     * communicate with the Gyoto::Factory to build an instance of
     * the class from its XML description
     */
    void Register(std::string kind, Gyoto::Spectrum::Subcontractor_t* scp);

    /// Query the Spectrum register Spectrum::Register_
    /**
     * Query the Spectrum register to get the Spectrum::Subcontractor_t
     * correspondig to a given kind name. This function is normally
     * called only from the Factory. If \p plugin is not empty, all
     * the plug-ins listed there will be first loaded using
     * Gyoto::requirePlugin(), then a subcontractor matching both \p
     * name and and one element of \p plugin will be searched for. If
     * \p plugin is an empty vector, then the first subcontractor
     * matching \p name will be returned, and the name of the plug-in
     * it belongs to will be returned in \p plugin upon output.
     *
     * This function is defined using the #GYOTO_GETSUBCONTRACTOR
     * macro.
     *
     * \param[in] name e.g. "PowerLaw"
     * \param[inout] plugin  vector of strings listing plug-ins to look
     *        for \p name in.
     * \param[in] errmode int=0. If errmode==0, failure to find a
     *        registered Spectrum by that name is an error. Else, simply
     *        return NULL pointer in that case.
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Spectrum::Subcontractor_t* getSubcontractor(std::string name,
						       std::vector<std::string> &plugins,
						       int errmode=0);

    /// The Spectrum register
    /**
     * Use the Spectrum::initRegister() once in your program to
     * initiliaze it, the Spectrum::Register() function to fill it, and
     * the Spectrum::getSubcontractor() function to query it.
     */
    extern Register::Entry* Register_;

    /// Empty the Spectrum register Spectrum::Register_
    /**
     *  This must be called once. It is called by
     *  Gyoto::Register::init().
     */
    void initRegister();
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
class Gyoto::Spectrum::Generic
: public Gyoto::SmartPointee,
  public Gyoto::Object
{
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::Generic>;
 protected:

 public:
  GYOTO_OBJECT;

  Generic();
  Generic(const std::string kind); ///< Set kind in constructor
  Generic(const Spectrum::Generic &);
  virtual Generic * clone() const; ///< Cloner

  virtual ~Generic() ; ///< Destructor: does nothing.

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

};


#endif
