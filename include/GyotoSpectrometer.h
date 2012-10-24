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
  class Spectrometer;
}

class Gyoto::Spectrometer : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer>;
 protected:
  SpectroKind_t kind_; ///< none, freqlog, freq, wavelog or wave.
  size_t nsamples_; ///< number of spectral elements
  double band_[2]; ///< boundaries of the spectro 
  double* boundaries_;
  double* midpoints_;
  double* widths_;

  void reset_(); ///< Computes boundaries_, midpoints_ and widths_

 public:
  Spectrometer() ; ///< Default constructor
  Spectrometer(size_t nsamples, double band_min, double band_max,
	       SpectroKind_t kind); ///< Constructor setting everything
  Spectrometer(const Spectrometer& ) ;                ///< Copy constructor
  Spectrometer * clone() const; ///< Cloner
  virtual ~Spectrometer() ; ///< Destructor

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
  double const * getBand() const ;

  double const * getMidpoints() const ;
  double const * getChannels() const ;
  double const * getWidths() const ;

#ifdef GYOTO_USE_XERCES
 public:
    void fillElement(FactoryMessenger *fmp); ///< called from Factory
#endif

};


#ifdef GYOTO_USE_XERCES
namespace Gyoto {
  SmartPointer<Spectrometer> SpectrometerSubcontractor(FactoryMessenger* fmp);
}
#endif


#endif
