/**
 *  \file GyotoSpectrometer.h
 *  \brief Spectroscopic capabilities of a Screen
 *
 *  Describes the spectroscopic capabilites of a Screen.
 *
 */
/*
    Copyright 2011-2016 Thibaut Paumard

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

#ifndef __GyotoSpectrometer_H_ 
#define __GyotoSpectrometer_H_ 

#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>
#include <GyotoObject.h>
#include <GyotoRegister.h>
#include <GyotoHooks.h>
#include <string>

/**
 * \namespace Gyoto::Spectrometer
 * \brief Access to spectrometers
 * 
 * Objects which describe spectrometers (including one-channel
 * devices, a.k.a cameras) must inherit from the
 * Gyoto::Spectrometer::Generic class.
 *
 * To be usable, a Spectrometer::Generic sub-class should register a
 * Spectrometer::Subcontractor_t function using the Spectrometer::Register()
 * function. See also \ref writing_plugins_page .
 */

namespace Gyoto{
  namespace Register { class Entry; }
  class FactoryMessenger;
  namespace Spectrometer {
    /**
     * \class Gyoto::Spectrometer::Generic
     *
     * \brief Base class for spectrometers
     *
     * Example: class Gyoto::Spectrometer::Uniform
     *
     * See Gyoto::Spectrometer for an introduction.
     *
     * Generic inherits from Gyoto::SmartPointee so that it is
     * possible to create a SmartPointer to a Spectrometer.
     *
     * It also inherits from Gyoto::Hook::Teller. This allows a
     * consistent implementation of Spectrometer::Complex (in
     * particular). Any method which mutates a Spectrometer should
     * call tellListeners().
     */
    class Generic;

  /**
   * \brief Type for Spectrometer kind
   * 
   * Spectrometer kindid is a unique numerical identifier for that kind,
   * produced as the address to a static C string variable holding the
   * kind's name. Most of the time, the address is the only
   * significant part as this is more reliable and allows for direct
   * numerical comparison instead of slower string comparison. The
   * value of the string variable can be used for printing and as
   * anchor name for Register(), although the anchor name could be
   * different.
   */
    typedef char const * kind_t;

    /**
     * This is a more specific version of the
     * SmartPointee::Subcontractor_t type. A Spectrometer::Subcontrator_t
     * is called by the Gyoto::Factory to build an instance of the
     * kind of spectrometer specified in an XML file (see
     * Register()). The Factory and Subcontractor_t function
     * communicate through a Gyoto::FactoryMessenger.
     */
    typedef SmartPointer<Gyoto::Spectrometer::Generic>
      Subcontractor_t(Gyoto::FactoryMessenger*, std::vector<std::string> const &);
    ///< A function to build instances of a specific Astrobj::Generic sub-class

    /**
     * \brief Query the Spectrometer register
     *
     * Get the Spectrometer::Subcontractor_t correspondig to a given
     * kind name. This function is normally called only from the
     * Gyoto::Factory. If plugin is specified, only a
     * subcontractor matching both name and plugin will be returned,
     * loading the plug-in if necessary. If plugin is the empty
     * string, then the first subcontractor matching name will be
     * returned, and the name of the plug-in it belongs to will be
     * returned in plugin upon output.
     *
     * \param[in] name Name of the subclass to build, e.g. "Complex"
     *            or "wave".
     * \param[inout] plugin e.g. "stdplug".
     * \param errmode[in] If name is not registered,
     *             getSubcontractor() return NULL errmode==1, throws a
     *             Gyoto::Error if errmode==0.
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Spectrometer::Subcontractor_t* getSubcontractor(std::string name,
							   std::vector<std::string> &plugins,
							   int errmode = 0);

    /**
     * \brief A template for Subcontractor_t functions
     *
     * Instead of reimplementing the wheel, your subcontractor can
     * simply be Gyoto::Spectrometer::Subcontractor<MyKind>. It must
     * however implement setParameters().
     *
     * \tparam T A Spectrometer::Generic sub-class.
     */
    template<typename T> SmartPointer<Spectrometer::Generic> Subcontractor
      (FactoryMessenger* fmp, std::vector<std::string> const &plugins) {
      SmartPointer<T> spectro = new T();
      spectro -> plugins(plugins);
#ifdef GYOTO_USE_XERCES
      if (fmp) spectro -> setParameters(fmp);
#endif
      return spectro;
    }

    /**
     * \brief The Spectrometer register
     *
     * Use the Spectrometer::initRegister() once in your program to
     * initiliaze it, the Spectrometer::Register() function to fill it, and
     * the Spectrometer::getSubcontractor() function to query it.
     */
    extern Gyoto::Register::Entry * Register_;

     /**
      * \brief Initialize the Spectrometer register
      *  This must be called once. It initializes Register_ and
      *  registers the standard kinds (Uniform and Complex).
      */
    void initRegister(); 

    /**
     * \brief Register a new Spectrometer kind
     *
     * Register a new Spectrometer::Generic sub-class so that the
     * Gyoto::Factory knows it.
     *
     * \param name The kind name which identifies this object type in
     * an XML file, as in &lt;Spectrometer kind="name"&gt;. For
     * clarity, this should be the same as the value of kindid_ for this
     * object, but it is not mandatory.
     *
     * \param scp A pointer to the subcontractor, which will
     * communicate with the Gyoto::Factory to build an instance of
     * the class from its XML description. If all parameters can be
     * set using setParameter(), this can be:
     * \code
     * &(Gyoto::Spectrometer::Subcontractor<Myind>)
     * \endcode
     */
    void Register(std::string name, Gyoto::Spectrometer::Subcontractor_t* scp);

  }
}

class Gyoto::Spectrometer::Generic
: public Gyoto::SmartPointee,
  public Gyoto::Object,
  public Gyoto::Hook::Teller
{
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>;
 protected:
  /**
   * \brief Spectrometer kind ID
   *
   * The content is not copied. kindid_ should be set (as a parameter to
   * the Generic() constructor or using kindid()) to the address of a
   * static variable holding the name. This allows checking the kind
   * using pointer comparison rather than string comparison.
   */
  kind_t kindid_;
 public:
  size_t nsamples_; ///< Number of spectral elements
  size_t nboundaries_; ///< Size of the boundaries_ array

  /**
   * \brief Frequency  (in Hz) at the boundaries of the spectral channels
   *
   * Array of size nboundaries_
   *
   * Spectral channel i extends from
   * \code
   * boundaries_[chanind_[2*i]]
   * \endcode
   * to
   * \code
   * boundaries_[chanind_[2*i+1]]
   * \endcode.
   * Channels may or may not be contiguous or ordered.
   */
  double* boundaries_;

  /**
   * \brief Indices in boundaries_
   *
   * Array of size 2*nsamples_
   */
  size_t* chanind_;

  /**
   * \brief Effective frequency (in Hz) of each spectral channel
   *
   * Array of size nsamples_
   */
  double* midpoints_;

  /**
   * \brief Width of each channel
   *
   * Array of size nsamples_ (in Hz)
   */
  double* widths_;

 public:
  GYOTO_OBJECT;

  /**
   * \brief Default constructor
   *
   * Sets each member to 0. 
   */
  Generic();

  /**
   * \brief Constructor setting kindid
   *
   * Sets the other members to 0. This is usually the right
   * constructor to use:
   * \code
   * Complex::Complex : Generic(Complex::Kindid) {}
   * \endcode
   *
   * Always set kindid to the address of a static variable, not to a temporary.
   * Usually your class should have a static member for that purpose:
   * \code
   * class MyKind : public Spectrometer::Generic
   * {
   *   static kind_t Kindid;
   * };
   * kind_t MyKind::Kindid = "MyKind";
   * \endcode
   *
   */
  Generic(kind_t kindid);

  /**
   * \brief Copy constructor
   *
   * Takes care of (deep) copying all the members known to the base class.
   */
  Generic(const Generic& ) ;

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
  virtual Generic * clone() const =0;

  /**
   * \brief Destructor
   *
   * Takes care of deleting the arrays (if the pointers are not NULL).
   */
  virtual ~Generic();

  /**
   * \brief Get kindid_
   *
   * You can check whether the Spectrometer sp is of a given kind
   * MyKind with something like:
   *
   * \code
   * if (sp->kind()) == MyKind::Kindid;
   * \endcode
   *
   * See Uniform::WaveKind, Uniform::WaveLogKind, Uniform::FreqKind,
   * Uniform::FreqLogKind and Complex::Kind.
   *
   */
  virtual kind_t kindid() const ;

  /**
   * \brief Set Generic::kindid_
   *
   * This should rarely be used as the Generic::kindid_ attribute usually is set
   * in the constructor and doesn't change after that.
   *
   * Always set to the address of a static variable, not to a temporary.
   * Usually your class should have a static member for that purpose:
   * \code
   * class MyKind : public Spectrometer::Generic
   * {
   *   static kind_t Kindid;
   * };
   * kind_t MyKind::Kind = "MyKind";
   * ...
   * SmartPointer<MyKind> sp();
   * sp->kindid(MyKind::Kindid)
   * \endcode
   * 
   */
  virtual void  kindid(kind_t) ;

  virtual size_t nSamples() const ; ///< Get Generic::nsamples_.
  virtual size_t getNBoundaries() const ; ///< Get Generic::nboundaries_
  virtual double const * getMidpoints() const  ; ///< Get Generic::midpoints_.
  /**
   * \brief Copy Generic::midpoints_, converting to unit
   * \param data an array of Generic::nsamples_ doubles to fill with result
   * \param unit a string 
   */
  virtual void getMidpoints( double data[], std::string unit);
  /**
   * \brief Copy Generic::boundaries_, converting to unit
   * \param data an array of Generic::nboundaries_ doubles to fill with result
   * \param unit a string 
   */
  virtual void getChannelBoundaries( double data[], std::string unit);
  virtual double const * getChannelBoundaries() const ; ///< Get Generic::boundaries_.
  virtual size_t const * getChannelIndices() const ; ///< Get Generic::chanind_.
  virtual double const * getWidths() const ; ///< Get Generic::widths_.
  /**
   * \brief Copy Generic::widths_, converting to unit
   *
   * Think carefully before using: widths are often used to convert
   * spectral flux density to flux. If flux density is per Herz, you
   * don't need to convert widths.
   *
   * \param data an array of Generic::nboundaries_ doubles to fill with result
   * \param unit a string 
   */
  virtual void getWidths( double data[], std::string unit);

};


#endif
