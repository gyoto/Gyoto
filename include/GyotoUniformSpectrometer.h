/**
 *  \file GyotoUniformSpectrometer.h
 *  \brief Uniformly spaced spectrometers
 *
 *  Spectral channels are contiguous and uniformly spaced in either
 *  wavelength, frequency or log10 of
 *  either. Gyoto::Spectrometer::Uniform is registered four times in
 *  the factory: as kind="wave", "wavelog", "freq" and
 *  "freqlog". Example XML entity:
 *  \code
 *  <Spectrometer kind="wave" unit="µm" nsamples=10>
 *   2.0 2.4
 *  </Spectrometer>
 *  \endcode
 *  
 *  The content of the entity yields the band pass expressed in "unit"
 *  or in log10(unit).
 *
 */
/*
    Copyright 2011-2013 Thibaut Paumard

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

#ifndef __GyotoUniformSpectrometer_H_ 
#define __GyotoUniformSpectrometer_H_ 

#include <GyotoSpectrometer.h>

namespace Gyoto{
  namespace Spectrometer {
    class Uniform;
  }
}

/**
 *  \class Gyoto::Spectrometer::Uniform
 *  \brief Uniformly spaced spectrometers
 *
 *  Spectral channels are contiguous and uniformly spaced in either
 *  wavelength, frequency or log10 of
 *  either. Gyoto::Spectrometer::Uniform is registered four times in
 *  the factory: as kind="wave", "wavelog", "freq" and
 *  "freqlog". Example XML entity:
 *  \code
 *  <Spectrometer kind="wave" unit="µm" nsamples=10>
 *   2.0 2.4
 *  </Spectrometer>
 *  \endcode
 *  
 *  The content of the entity yields the band pass expressed in "unit"
 *  or in log10(unit).
 */
class Gyoto::Spectrometer::Uniform : public Gyoto::Spectrometer::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrometer::Uniform>;
 protected:
  /**
   * \brief boundaries of the spectro.
   *
   * Depending on the kind, band_ is stored in Hz, log10(Hz), m or log10(m).
   */
  double band_[2];

  void reset_(); ///< Computes boundaries_, midpoints_ and widths_

 public:
  GYOTO_OBJECT;
  Uniform() ; ///< Default constructor
  Uniform(size_t nsamples, double band_min, double band_max,
	       kind_t kind); ///< Constructor setting everything
  Uniform(const Uniform& ) ;                ///< Copy constructor
  Generic * clone() const; ///< Cloner
  virtual ~Uniform() ; ///< Destructor

  using Generic::kindid;
  void kindid(kind_t);

  /**
   * \brief Set Generic::kind_ from a std::string
   *
   * Generic::kind_ will actually be set to one of Uniform::WaveKind,
   * Uniform::WaveLogKind, Uniform::FreqKind or Uniform::FreqLogKind.
   *
   * \param name std::string, one of "wave", "wavelog", "freq" or
   * "freqlog"
   */
  void kind(std::string const &name);
  std::string kind() const;
 
  using Generic::nSamples;
 /**
   * \brief Set Generic::nsamples_
   */
  void nSamples(size_t n);
 
  void band(std::vector<double> const &nu);
  void band(std::vector<double> const &nu, std::string const &unit);
  std::vector<double>band() const;
  std::vector<double>band(std::string const &unit) const;

 /**
   * \brief Set Uniform::band_
   *
   * \param nu 2-element vector, in Hz, m, log10(Hz) or log10(m)
   * depending on Generic::kind_
   */
  void band(double nu[2]);

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
  void band(double nu[2], std::string const &unit, std::string const &kind);
  void band(double nu[], std::string const &unit);

  double const * getBand() const ; ///< Get Uniform::band_

 public:
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;
#ifdef GYOTO_USE_XERCES
  virtual void setParameters(FactoryMessenger *fmp);
#endif

  /**
   * \brief "wave"
   *
   * Use this static member attribute to check whether a Spectrometer
   * object spectro is of this kind:
   * \code
   * if (spectro->kind() == Uniform::WaveKind) ... ;
   * \endcode
   *
   */
    static kind_t const WaveKind; ///< "wave"

  /**
   * \brief "wavelog"
   *
   * Use this static member attribute to check whether a Spectrometer
   * object spectro is of this kind:
   * \code
   * if (spectro->kind() == Uniform::WaveLogKind) ... ;
   * \endcode
   *
   */
    static kind_t const WaveLogKind;

  /**
   * \brief "freq"
   *
   * Use this static member attribute to check whether a Spectrometer
   * object spectro is of this kind:
   * \code
   * if (spectro->kind() == Uniform::FreqKind) ... ;
   * \endcode
   *
   */
    static kind_t const FreqKind;

  /**
   * \brief "freqlog"
   *
   * Use this static member attribute to check whether a Spectrometer
   * object spectro is of this kind:
   * \code
   * if (spectro->kind() == Uniform::FreqLogKind) ... ;
   * \endcode
   *
   */
    static kind_t const FreqLogKind;


};


#endif
