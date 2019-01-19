/**
 * \file GyotoPatternDisk.h
 * \brief A geometrically thin, optically thick disk
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2011-2015, 2018 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoPatternDisk_H_ 
#define __GyotoPatternDisk_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class PatternDisk; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>

/**
 * \class Gyoto::Astrobj::PatternDisk
 * \brief Geometrically thin disk read from FITS file
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane, extending from r=r_ISCO to r=infinity.  The flux emitted
 *   at radius r and longitude phi at frequency nu is given in a FITS
 *   file.
 *
 */
class Gyoto::Astrobj::PatternDisk : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PatternDisk>;
 private:
  std::string filename_; ///< Optional FITS file name containing the arrays
  /**
   * An array of dimensionality double[nr_][nphi_][nnu_]. In FITS
   * format, the first dimension is nu, the second phi, and the third
   * r.
   */
  double * emission_; ///< I<SUB>&nu;</SUB>(&nu;, r, &phi;)

  double * opacity_; ///< Same dimenstions as emission, or NULL

  /**
   * An array of dimensionality double[nr_][nphi_][2]. In FITS format,
   * the second dimension is phi, and the third r. The first plane in
   * the first FITS dimention is d&phi;/dt, the second dr/dt.
   */
  double * velocity_; ///< velocity(r, &phi;)

  /**
   * In case of adaptive grid.
   */
  double * radius_; ///< Radius vector

  /**
   * XML element: &lt;Omega&gt;.
   * FITS keyword: HIERARCH GYOTO PatternDisk Omega
   */
  double Omega_;  ///< Pattern angular velocity

  /**
   * XML element: &lt;T0&gt;.
   * FITS keyword: HIERARCH GYOTO PatternDisk t0
   */
  double t0_;     ///< Date for which i=0 corresponds to phi=0

  double dnu_; ///< Frequency scale of PatternDisk::emission_ in Hz
  double nu0_; ///< Lowest frequency provided in PatternDisk::emission_ in Hz
  size_t nnu_; ///< Number of frequencies provided in PatternDisk::emission_

  double dphi_; ///< &delta;&phi; between two grid columns
  double phimin_;///< Minimum &phi; in grid
  size_t nphi_; ///< Grid size in the &phi; direction
  double phimax_; ///< Maximum &phi; in grid

  /**
   * XML elment: &lt;RepeatPhi&gt;.
   * FITS keyword: HIERARCH GYOTO PatternDisk RepeatPhi
   */
  size_t repeat_phi_; ///< Number of times the pattern should be repeated to cover [0, 2&Pi;]

  double dr_; ///< Radius step
  size_t nr_; ///< Number of rows in the patternGrid size in the r direction


  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  // fillProperty is overridden to remove leading "!" from FITS filename
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;

  PatternDisk(); ///< Standard constructor
  
  PatternDisk(const PatternDisk& ) ;///< Copy constructor
  virtual PatternDisk* clone () const; ///< Cloner
  
  virtual ~PatternDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  using ThinDisk::innerRadius;
  virtual void   innerRadius(double);
  using ThinDisk::outerRadius;
  virtual void   outerRadius(double);

  /**
   * Unit: radians per geometrical unit time.
   */
  virtual void   patternVelocity(double); ///< Set PatternDisk::Omega_
  virtual double patternVelocity() const; ///< Get PatternDisk::Omega_

  virtual void file(std::string const &f);
  virtual std::string file() const ;

#ifdef GYOTO_USE_CFITSIO
  /// Read parameters and arrays from FITS file
  virtual void fitsRead(std::string filename_);

  /// Write parameters and arrays to FITS file
  virtual void fitsWrite(std::string filename_);
#endif

  /// Set PatternDisk::emission_
  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setEmission(double * pattern);

  /// Set PatternDisk::velocity__
  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setVelocity(double * pattern);

  /// Set PatternDisk::radius_
  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void radius(double * pattern);

  /// Set PatternDisk::emission_
  /**
   * PatternDisk::emission_ is freed if not NULL, reallocated, and
   * pattern is copied into emission_.
   *
   * If PatternDisk::opacity_, PatternDisk::velocity_ or
   * PatternDisk::radius_ have been set previously with mismatching
   * sizes, they are deallocated too.
   *
   * Finally, PatternDisk::nnu_, PatternDisk::nphi_, and
   * PatternDisk::nr_ are set according to naxes.
   *
   * \param pattern Array to copy as emission_. May be NULL in which
   * case emission_ is simply deallocated and set to NULL.
   *
   * \param naxes { nnu_, nphi_, nr_ }.
   */
  virtual void copyIntensity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);

  virtual double const * getIntensity() const;///< Get PatternDisk::emission_
  virtual void getIntensityNaxes( size_t naxes[3] ) const ; ///< Get PatternDisk::nnu_, PatternDisk::nphi_, and PatternDisk::nr_

  /**
   * \brief Set PatternDisk::opacity_
   *
   * PatternDisk::opacity_ is first freed if not NULL and set to NULL.
   *
   * If pattern is not NULL, PatternDisk::emission_ must have been set
   * previously with matching dimensions. PatternDisk::opacity_ is
   * then reallocated, and pattern is copied into opacity_.
   *
   * \param pattern Array to copy as opacity_. May be NULL in which
   * case opacity_ is simply deallocated and set to NULL.
   *
   * \param naxes { nnu_, nphi_, nr_ }.
   */
  virtual void copyOpacity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);
  virtual double const * opacity() const; ///< Get PatternDisk::opacity_

  /// Set PatternDisk::velocity_
  /**
   * PatternDisk::velocity_ is first freed if not NULL and set to NULL.
   *
   * If pattern is not NULL, PatternDisk::emission_ must have been set
   * previously with matching dimensions. PatternDisk::velocity_ is
   * then reallocated, and pattern is copied into velocity_.
   *
   * \param pattern Array to copy as velocity_. May be NULL in which
   * case velocity_ is simply deallocated and set to NULL.
   *
   * \param naxes { nphi_, nr_ }.
   */
  virtual void copyVelocity(double const * const pattern = NULL,
			      size_t const naxes[2] = NULL);
  virtual double const * getVelocity() const;///< Get PatternDisk::velocity_

  /// Set PatternDisk::radius_
  /**
   * PatternDisk::radius_ is first freed if not NULL and set to NULL.
   *
   * If pattern is not NULL, PatternDisk::emission_ must have been set
   * previously with matching dimensions. PatternDisk::radius_ is
   * then reallocated, and pattern is copied into radius_.
   *
   * \param pattern Array to copy as radius_. May be NULL in which
   * case radius_ is simply deallocated and set to NULL.
   *
   * \param nr size of radius array.
   */
  virtual void copyGridRadius(double const * const pattern = NULL,
			      size_t nr = 0 );
  virtual double const * getGridRadius() const; ///< Get PatternDisk::radius_

  virtual void repeatPhi(size_t n); ///< Set PatternDisk::repeat_phi_
  virtual size_t repeatPhi() const; ///< Get PatternDisk::repeat_phi_

  virtual void nu0(double freq); ///< Set PatternDisk::nu0_
  virtual double nu0() const; ///< Get PatternDisk::nu0_

  virtual void dnu(double dfreq); ///< Set PatternDisk::dnu_
  virtual double dnu() const; ///< Get PatternDisk::dnu_

  void phimin(double phimin); ///< Set PatternDisk::phimin_
  double phimin() const; ///< Get PatternDisk::phimin_

  void phimax(double phimax); ///< Set PatternDisk::phimax_
  double phimax() const; ///< Get PatternDisk::phimax_

 protected:
  void getIndices(size_t i[3], double const co[4], double nu=0.) const ;
  ///< Get emission_ cell corresponding to position co[4]

 public:
  using ThinDisk::emission;
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;
  virtual double transmission(double nu_em, double dsem, state_t const &, double const coord[8]) const;

  virtual void getVelocity(double const pos[4], double vel[4])  ;

};

#endif
