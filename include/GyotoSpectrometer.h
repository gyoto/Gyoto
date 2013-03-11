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
  }
}

class Gyoto::Spectrometer::Generic : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>;
 public:
  Generic();
  Generic(const Generic& ) ;                ///< Copy constructor
  virtual ~Generic();
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
  SpectroKind_t kind_; ///< none, freqlog, freq, wavelog or wave.
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

  SpectroKind_t getKind() const ;
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
#endif

};


#ifdef GYOTO_USE_XERCES
namespace Gyoto {
  SmartPointer<Spectrometer::Generic> SpectrometerSubcontractor(FactoryMessenger* fmp);
}
#endif


#endif
