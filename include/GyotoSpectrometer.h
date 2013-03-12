/**
 *  \file GyotoSpectrometer.h
 *  \brief Spectroscopic capabilities of a Screen
 *
 *  Describes the spectroscopic capabilites of a Screen.
 *
 */
/*
    Copyright 2011 Thibaut Paumard

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
#include <GyotoRegister.h>
#include <string>

namespace Gyoto{
  namespace Spectrometer {
    class Generic;
    class Uniform;
    /**
     * This is a more specific version of the
     * SmartPointee::Subcontractor_t type. A Spectrometer::Subcontrator_t
     * is called by the Gyoto::Factory to build an instance of the
     * kind of spectrometer specified in an XML file (see
     * Register()). The Factory and Subcontractor_t function
     * communicate through a Gyoto::FactoryMessenger.
     */
    typedef SmartPointer<Gyoto::Spectrometer::Generic>
      Subcontractor_t(Gyoto::FactoryMessenger*);
    ///< A function to build instances of a specific Astrobj::Generic sub-class

    /**
     * Query the Spectrometer register to get the Spectrometer::Subcontractor_t
     * correspondig to a given kind name. This function is normally
     * called only from the Factory.
     *
     * \param name e.g. "Star"
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Spectrometer::Subcontractor_t* getSubcontractor(std::string name,
						      int errmode = 1);
    ///< Query the Spectrometer register

#if defined GYOTO_USE_XERCES
    /**
     * Use the Spectrometer::initRegister() once in your program to
     * initiliaze it, the Spectrometer::Register() function to fill it, and
     * the Spectrometer::getSubcontractor() function to query it.
     */
    extern Gyoto::Register::Entry * Register_;
    ///< The Spectrometer register

     /**
      *  This must be called once.
      */
    void initRegister(); 
    ///< Empty the Spectrometer register

    /**
     * Register a new Spectrometer::Generic sub-class so that the
     * Gyoto::Factory knows it.
     *
     * \param name The kind name which identifies this object type in
     * an XML file, as in &lt;Spectrometer kind="name"&gt;
     *
     * \param scp A pointer to the subcontractor, which will
     * communicate whith the Gyoto::Factory to build an instance of
     * the class from its XML description
     */
    void Register(std::string name, Gyoto::Spectrometer::Subcontractor_t* scp);
    ///< Make an Spectrometer kind known to the Factory
#endif

  }
}

class Gyoto::Spectrometer::Generic : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>;
 protected:
  char const * kind_;
 public:
  Generic();
  Generic(SpectroKind_t kind);
  Generic(const Generic& ) ;                ///< Copy constructor
  virtual ~Generic();
  virtual char const * getKind() const;
  virtual void  setKind(char const *) ;
  virtual size_t getNSamples() const =0;
  virtual size_t getNBoundaries() const =0;
  virtual double const * getMidpoints() const =0 ;
  virtual double const * getChannelBoundaries() const =0;
  virtual size_t const * getChannelIndices() const =0;
  virtual double const * getWidths() const =0;
  virtual Generic * clone() const =0;
#ifdef GYOTO_USE_XERCES

  //virtual void setParameters(Gyoto::FactoryMessenger *fmp) ;

  /**
   * Metrics implementations should impement fillElement to save their
   * parameters to XML and call the Metric::fillElement(fmp) for the
   * shared properties
   */

  virtual void fillElement(FactoryMessenger *fmp) =0; ///< called from Factory
  //void processGenericParameters(Gyoto::FactoryMessenger *fmp) ;
#endif
};

class Gyoto::Spectrometer::Uniform : public Gyoto::Spectrometer::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer::Uniform>;
 protected:
  size_t nsamples_; ///< number of spectral elements
  double band_[2]; ///< boundaries of the spectro 
  /**
   * Spectral channel i extends from
\code
boundaries_[chanind_[2*i]]
\endcode
to
\code
boundaries_[chanind_[2*i+1]]
\endcode.
Channels may or may not be contiguous or ordered.
   */
  double* boundaries_; ///< Frequency at the boundaries of the spectral channels
  size_t* chanind_; ///< Indices in boundaries_
  double* midpoints_;
  double* widths_;

  void reset_(); ///< Computes boundaries_, midpoints_ and widths_

 public:
  Uniform() ; ///< Default constructor
  Uniform(size_t nsamples, double band_min, double band_max,
	       SpectroKind_t kind); ///< Constructor setting everything
  Uniform(const Uniform& ) ;                ///< Copy constructor
  Generic * clone() const; ///< Cloner
  virtual ~Uniform() ; ///< Destructor

  void setKind(SpectroKind_t);
  void setKind(std::string);
  void setNSamples(size_t n);
  void setBand(double nu[2]);

  /**
   * \brief Set the spectral band boundaries in specified unit
   *
   * If kind is not specified, member kind_ is used. Else kind_ is updated.
   *
   * unit is actually the unit for 10^nu for freqlog and wavelog. Defaults:
   *  - kind==freq: nu in Hz
   *  - kind==freqlog: 10^nu in Hz
   *  - kind==wave: nu in meters
   *  - kind==wavelog: 10^nu in meters
   * 
   */
  void setBand(double nu[2], std::string unit, std::string kind="");

  std::string getKindStr() const;
  size_t getNSamples() const ;
  size_t getNBoundaries() const ;
  double const * getBand() const ;

  double const * getMidpoints() const ;
  double const * getChannelBoundaries() const ;
  size_t const * getChannelIndices() const ;
  double const * getWidths() const ;

#ifdef GYOTO_USE_XERCES
 public:
    void fillElement(FactoryMessenger *fmp); ///< called from Factory
    static Spectrometer::Subcontractor_t Subcontractor;
#endif

    static SpectroKind_t const WaveKind;
    static SpectroKind_t const WaveLogKind;
    static SpectroKind_t const FreqKind;
    static SpectroKind_t const FreqLogKind;


};


#endif
