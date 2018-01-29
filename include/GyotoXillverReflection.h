/**
 * \file GyotoXillverReflection.h
 * \brief The illumination table specifies how the thin disk is
 * illuminated while the reflection table deduces from that
 * the reflected spectrum as computed by Javier Garcia's XILLVER code.
 * The metric is imposed to be Kerr for simplicity.
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


#ifndef __GyotoXillverReflection_H_
#define __GyotoXillverReflection_H_

#include <GyotoThinDisk.h>


namespace Gyoto{
  namespace Astrobj { class XillverReflection; }
}

/**
 * \class Gyoto::Astrobj::XillverReflection
 * \brief The illumination table specifies how the thin disk is
 * illuminated while the reflection table deduces from that
 * the reflected spectrum as computed by Javier Garcia's XILLVER code.
 *
 */
class Gyoto::Astrobj::XillverReflection
: public Astrobj::ThinDisk,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::XillverReflection>;

 private:
  std::string filenameIllum_; ///< FITS file containing the illumination pattern
  std::string filenameRefl_; ///< FITS file containing the reflection pattern 
  /**
   * An array of dimensionality double[nxi_][ni_][nnu_]. In FITS
   * format, the first dimension is nu, the second is incl (the emission angle),
   * and the third is log(ionization parameter).
   */
  double * reflection_; 

  double * logxi_; ///< log of ionization param
  double * incl_; ///< emission angle
  double * freq_; ///< frequencies vector
  size_t nnu_; ///< Number of frequencies
  size_t ni_; ///< Number of emission angles
  size_t nxi_; ///< Number of log(ionization param)

  double * illumination_;

  double * radius_; ///< radii at which illumination is known
  double * phi_; ///< azimuthal angle at which illumination is known
  size_t nr_; ///< numbar of radii
  size_t nphi_; ///< numbar of phi

  double aa_; ///< Spin of Kerr BH
  double lampradius_; ///< Coordinate radius at which the lamp is in Keplerian rotation
  double timelampphizero_; ///< Time at which lamp is at phi=0

  bool average_over_angle_; ///< true to average over emission angle

 protected:

 public:
  GYOTO_OBJECT;
  // fillProperty is overridden to remove leading "!" from FITS filename
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;
  XillverReflection(); ///< Standard constructor
  XillverReflection(const XillverReflection& o); ///< Copy constructor
  virtual XillverReflection * clone() const ; ///< Cloner
  virtual ~XillverReflection() ; ///< Destructor

 public:

  void timelampphizero(double tt);
  double timelampphizero() const;
  void lampradius(double rr);
  double lampradius() const;
  void fileillumination(std::string const &f);
  std::string fileillumination() const ;
  void filereflection(std::string const &f);
  std::string filereflection() const ;
  
  void averageOverAngle(bool t);
  bool averageOverAngle()const;

  #ifdef GYOTO_USE_CFITSIO
  /// Read parameters and arrays from FITS file
  virtual void fitsReadIllum(std::string filename);
  /// Write parameters and arrays to FITS file
  virtual void fitsWriteIllum(std::string filename);

  /// Read parameters and arrays from FITS file
  virtual void fitsReadRefl(std::string filename);
  /// Write parameters and arrays to FITS file
  virtual void fitsWriteRefl(std::string filename);
#endif

  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setIllumination(double * pattern);
  void setReflection(double * pattern);
  
  
  /**
   * XillverReflection::emission_ is freed if not NULL, 
   * reallocated, and
   * pattern is copied into emission_.
   *
   * Finally, XillverReflection::nnu_, 
   * XillverReflection::ni_, and
   * XillverReflection::nsg_ are set according to naxes.
   *
   * \param pattern Array to copy as emission_. May be NULL in which
   * case emission_ is simply deallocated and set to NULL.
   *
   * \param naxes { nnu_, ni_, nsg_ }.
   */
  virtual void copyIllumination(double const * const pattern = NULL,
				size_t const naxes[2] = NULL);
  virtual double const * getIllumination() const;
  virtual void getIlluminationNaxes( size_t naxes[2] ) const ; ///< Get XillverReflection::nr_, XillverReflection::nphi_

  virtual void copyReflection(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);
  virtual double const * getReflection() const;
  virtual void getReflectionNaxes( size_t naxes[3] ) const ; ///< Get XillverReflection::nnu_, XillverReflection::ni_, XillverReflection::nxi_

  virtual void copyGridReflLogxi(double const * const pattern = NULL,
				 size_t nxi = 0 );
  virtual double const * getGridReflLogxi() const; ///< Get XillverReflection::logxi_
  virtual void copyGridReflIncl(double const * const pattern = NULL,
				size_t ni = 0 );
  virtual double const * getGridReflIncl() const; ///< Get XillverReflection::incl_
  virtual void copyGridReflFreq(double const * const pattern = NULL,
				size_t nnu = 0 );
  virtual double const * getGridReflFreq() const; ///< Get XillverReflection::freq_
  virtual void copyGridIllumRadius(double const * const pattern = NULL,
				   size_t nr = 0 );
  virtual double const * getGridIllumRadius() const; ///< Get XillverReflection::radius_
  virtual void copyGridIllumPhi(double const * const pattern = NULL,
				size_t nphi = 0 );
  virtual double const * getGridIllumPhi() const; ///< Get XillverReflection::phi_
  
 protected:
  
  void getIndicesRefl(size_t i[3], double const co[4], double logxi, double incl,
		      double nu=0.) const ;
  ///< Get reflection_ cell corresponding to position co[4]
  void getIndicesIllum(size_t i[3], double const co[4]) const ;
  ///< Get illumination_ cell corresponding to position co[4]

 public:
  
  virtual double emission(double nu_em, double dsem,
			  state_t const &_ph, double const _obj[8]=NULL) const;

  virtual void updateSpin() ;
  virtual void tell(Gyoto::Hook::Teller *msg);
  using Generic::metric;
  virtual void metric(SmartPointer<Metric::Generic>);

};

#endif
