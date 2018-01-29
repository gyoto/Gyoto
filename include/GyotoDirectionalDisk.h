/**
 * \file GyotoDirectionalDisk.h
 * \brief Geometrically thin disk read from FITS file
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane.  The flux emitted
 *   at radius r, making an angle i with respect to the local normal,
 *   at frequency nu is given in a FITS file.
 *  
 *   This astrobj is typically used to compute reflected spectra
 *   in the lamp post model. 
 *
 *   For the time being the metric is imposed to be KerrBL, but should
 *   easily generalized.
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2014-2015, 2018 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoDirectionalDisk_H_ 
#define __GyotoDirectionalDisk_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class DirectionalDisk; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>

/**
 * \class Gyoto::Astrobj::DirectionalDisk
 * \brief Geometrically thin disk read from FITS file
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane.  The flux emitted
 *   at radius r, making an angle i with respect to the local normal,
 *   at frequency nu is given in a FITS file.
 *
 */
class Gyoto::Astrobj::DirectionalDisk : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::DirectionalDisk>;
 private:
  std::string filename_; ///< Optional FITS file name containing the arrays
  /**
   * An array of dimensionality double[nr_][ni_][nnu_]. In FITS
   * format, the first dimension is nu, the second cosi (direction cosine), 
   * and the third r. There is no phi dependence.
   */
  double * emission_; ///< I<SUB>&nu;</SUB>(&nu;, r, cosi;)

  double * radius_; ///< Radius vector
  double * cosi_; ///< Direction cosine vector
  double * freq_; ///< Frequencies vector

  double lampaltitude_; ///< Lamp altitude (z coordinate) in M units

  size_t nnu_; ///< Number of frequencies provided in DirectionalDisk::emission_
  size_t ni_; ///< Number of direction cosine
  size_t nr_; ///< Number of radius values

  double minfreq_computed_; ///< Minimum frequency computed by ATM21
  double maxfreq_computed_; ///< Maximum frequency computed by ATM21

  double minfreq_lampframe_; ///< Minimum frequency emitted by the lamp
  double maxfreq_lampframe_; ///< Maximum frequency emitted by the lamp

  bool average_over_angle_; ///< true to average over emission angle

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  // fillProperty is overridden to remove leading "!" from FITS filename
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;

  using Generic::metric;
  void metric(SmartPointer<Metric::Generic> gg);

  DirectionalDisk(); ///< Standard constructor
  
  DirectionalDisk(const DirectionalDisk& ) ;///< Copy constructor
  virtual DirectionalDisk* clone () const; ///< Cloner
  
  virtual ~DirectionalDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:

  void file(std::string const &f);
  std::string file() const ;
  void averageOverAngle(bool t);
  bool averageOverAngle()const;
  void lampaltitude(double zz);
  double lampaltitude() const ;
  void lampcutoffsinev(std::vector<double> const &v) ;
  std::vector<double> lampcutoffsinev() const ;

#ifdef GYOTO_USE_CFITSIO
  /// Read parameters and arrays from FITS file
  virtual void fitsRead(std::string filename_);

  /// Write parameters and arrays to FITS file
  virtual void fitsWrite(std::string filename_);
#endif

  /// Set DirectionalDisk::emission_
  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setEmission(double * pattern);

  void radius(double * pattern);

  /**
   * DirectionalDisk::emission_ is freed if not NULL, reallocated, and
   * pattern is copied into emission_.
   *
   * Finally, DirectionalDisk::nnu_, DirectionalDisk::ni_, and
   * DirectionalDisk::nr_ are set according to naxes.
   *
   * \param pattern Array to copy as emission_. May be NULL in which
   * case emission_ is simply deallocated and set to NULL.
   *
   * \param naxes { nnu_, ni_, nr_ }.
   */
  virtual void copyIntensity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);

  virtual double const * getIntensity() const;///< Get DirectionalDisk::emission_
  virtual void getIntensityNaxes( size_t naxes[3] ) const ; ///< Get DirectionalDisk::nnu_, DirectionalDisk::ni_, and DirectionalDisk::nr_


  virtual void copyGridRadius(double const * const pattern = NULL,
			      size_t nr = 0 );
  virtual double const * getGridRadius() const; ///< Get DirectionalDisk::radius_
  virtual void copyGridCosi(double const * const pattern = NULL,
			      size_t ni = 0 );
  virtual double const * getGridCosi() const; ///< Get DirectionalDisk::cosi_
  virtual void copyGridFreq(double const * const pattern = NULL,
			      size_t ni = 0 );
  virtual double const * getGridFreq() const; ///< Get DirectionalDisk::freq_

 protected:
  void getIndices(size_t i[3], double const co[4], double cosi, double nu=0.) const ;
  ///< Get emission_ cell corresponding to position co[4]

 public:
  using ThinDisk::emission;
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;

};

#endif
