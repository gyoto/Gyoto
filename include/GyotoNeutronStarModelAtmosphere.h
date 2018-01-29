/**
 * \file GyotoNeutronStarModelAtmosphere.h
 * \brief Neutron star emitting at its surface an analytic
 *  emission, typically blackbody
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
  Copyright (c) 2017, 2018 Frederic Vincent, Thibaut Paumard
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


#ifndef __GyotoNeutronStarModelAtmosphere_H_
#define __GyotoNeutronStarModelAtmosphere_H_

#include <GyotoStandardAstrobj.h>
#include <GyotoNumericalMetricLorene.h>
#include <GyotoNeutronStar.h>


namespace Gyoto{
  namespace Astrobj { class NeutronStarModelAtmosphere; }
}

/**
 * \class Gyoto::Astrobj::NeutronStarModelAtmosphere
 * \brief Neutron star emitting at its surface an 
 *  emission provided by a FITS table
 *
 */
class Gyoto::Astrobj::NeutronStarModelAtmosphere : public Astrobj::NeutronStar {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::NeutronStarModelAtmosphere>;

 private:
  std::string filename_; ///< Optional FITS file name containing the arrays
  /**
   * An array of dimensionality double[nsg_][ni_][nnu_]. In FITS
   * format, the first dimension is nu, the second cosi (direction cosine), 
   * and the third surface gravity.
   */
  double * emission_; ///< I<SUB>&nu;</SUB>(&nu;, surfgrav, cosi;)

  double * surfgrav_; ///< Surface gravity vector
  double * cosi_; ///< Direction cosine vector
  double * freq_; ///< Frequencies vector

  size_t nnu_; ///< Number of frequencies provided in NeutronStarModelAtmosphere::emission_
  size_t ni_; ///< Number of direction cosine
  size_t nsg_; ///< Number of surfgrav values

  bool average_over_angle_; ///< true to average over emission angle

 protected:

 public:
  GYOTO_OBJECT;
  // fillProperty is overridden to remove leading "!" from FITS filename
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;
  NeutronStarModelAtmosphere(); ///< Standard constructor
  NeutronStarModelAtmosphere(const NeutronStarModelAtmosphere& o); ///< Copy constructor
  virtual NeutronStarModelAtmosphere * clone() const ; ///< Cloner
  virtual ~NeutronStarModelAtmosphere() ; ///< Destructor

 public:

    void file(std::string const &f);
  std::string file() const ;
  void averageOverAngle(bool t);
  bool averageOverAngle()const;

  #ifdef GYOTO_USE_CFITSIO
  /// Read parameters and arrays from FITS file
  virtual void fitsRead(std::string filename_);

  /// Write parameters and arrays to FITS file
  virtual void fitsWrite(std::string filename_);
#endif

  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setEmission(double * pattern);
  
  void surfgrav(double * pattern);
  
  /**
   * NeutronStarModelAtmosphere::emission_ is freed if not NULL, 
   * reallocated, and
   * pattern is copied into emission_.
   *
   * Finally, NeutronStarModelAtmosphere::nnu_, 
   * NeutronStarModelAtmosphere::ni_, and
   * NeutronStarModelAtmosphere::nsg_ are set according to naxes.
   *
   * \param pattern Array to copy as emission_. May be NULL in which
   * case emission_ is simply deallocated and set to NULL.
   *
   * \param naxes { nnu_, ni_, nsg_ }.
   */
  virtual void copyIntensity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);

  virtual double const * getIntensity() const;///< Get NeutronStarModelAtmosphere::emission_
  virtual void getIntensityNaxes( size_t naxes[3] ) const ; ///< Get NeutronStarModelAtmosphere::nnu_, NeutronStarModelAtmosphere::ni_, and NeutronStarModelAtmosphere::nsg_

  virtual void copyGridSurfgrav(double const * const pattern = NULL,
			      size_t nsg = 0 );
  virtual double const * getGridSurfgrav() const; ///< Get NeutronStarModelAtmosphere::surfgrav_
  virtual void copyGridCosi(double const * const pattern = NULL,
			      size_t ni = 0 );
  virtual double const * getGridCosi() const; ///< Get NeutronStarModelAtmosphere::cosi_
  virtual void copyGridFreq(double const * const pattern = NULL,
			      size_t ni = 0 );
  virtual double const * getGridFreq() const; ///< Get NeutronStarModelAtmosphere::freq_

 protected:
  void getIndices(size_t i[3], double const co[4], double cosi, double nu=0.) const ;
  ///< Get emission_ cell corresponding to position co[4]

 public:
  
  virtual double emission(double nu_em, double dsem,
			  state_t const &_ph, double const _obj[8]=NULL) const;

};

#endif
