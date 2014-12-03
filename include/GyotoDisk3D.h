/**
 * \file GyotoDisk3D.h
 * \brief A geometrically thick, optically thin disk
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoDisk3D_H_ 
#define __GyotoDisk3D_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class Disk3D; }
}

/**
 * \class Gyoto::Astrobj::Disk3D
 * \brief Geometrically thick disk read from FITS file
 * 
 *   This class is the base class for thick disks.  The emitter's
 *   velocity is given in a FITS file, together with emission 
 *   related quantity (typically temperature).
 *   This class mainly implements the Impact() function.
 *   Emission() function is here left to its default, and should
 *   be implemented according to specific needs in heir classes.
 *   Here the disk is supposed not to evolve in time. The dynamical
 *   treatment is provided in heir classes.
 *
 *   The 3D disk is assumed to be described by a regular 
 *   (non adaptive) grid of cylindrical geometry. The disk
 *   is a slab from rin_ to rout_ and zmin_ (typically = -zmax_) to zmax_.
 *
 */
class Gyoto::Astrobj::Disk3D : public Gyoto::Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Disk3D>;
 private:
  std::string filename_; ///< Optional FITS file name containing the arrays
  /**
   * An array of dimensionality double[nr_][nz_][nphi_][nnu_]. In FITS
   * format, the first dimension is nu, the second phi, the third
   * z and last r. It typically contains temperature and is used only by
   * subclasses.
   */
  double * emissquant_; ///< Physical quantity yielding emission.

  double * opacity_; ///< Opacity, same dimensions as emissquant_

  /**
   * An array of dimensionality double[nr_][nz_][nphi_][3]. In FITS format,
   * the second dimension is phi, the third z and last r. The first plane in
   * the first FITS dimention is dphi/dt, the second dz/dt the last dr/dt.
   */
  double * velocity_; ///< Velocity(r, z, phi)

  double dnu_; ///< Frequency scale of PatternDisk::emission_ in Hz
  double nu0_; ///< Lowest frequency provided in PatternDisk::emission_ in Hz
  size_t nnu_; ///< Number of frequencies provided in PatternDisk::emission_

  double dphi_; ///< &delta;&phi; between two grid columns
  double phimin_;///< Minimum &phi; in grid
  size_t nphi_; ///< Grid size in the &phi; direction
  double phimax_; ///< Maximum &phi; in grid

  /**
   * XML elment: &lt;RepeatPhi&gt;.
   * FITS keyword: HIERARCH GYOTO Disk3D RepeatPhi
   */
  size_t repeat_phi_; ///< Number of times the pattern should be repeated to cover [0, 2&Pi;]
  //double phi0_==0, phi max is always 2*M_PI

  double dz_; ///< Altitude step
  double zmin_; ///< Minimum altitude
  size_t nz_; ///< Grid size in the altitude direction
  double zmax_; ///< Maximum altitude

  double dr_; ///< Radius step
  double rin_; ///< Inner radius of the grid
  size_t nr_; ///< Number of rows in the patternGrid size in the r direction
  double rout_; ///< Outer radius of the grid

  int zsym_; ///< 1 to symmetrize the grid z -> -z (default case)

  double tPattern_; ///< If the disk is being rotated (like a pattern disk) this is the origin of time for this rotation
  double omegaPattern_; ///< If the disk is being rotated (like a pattern disk) this is the rotation velocity dphi/dt

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  // fillProperty is overridden to remove leading "!" from FITS filename
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;

  Disk3D(); ///< Standard constructor.
  
  Disk3D(const Disk3D& ) ;///< Copy constructor.
  virtual Disk3D* clone () const; ///< Cloner.
  
  virtual ~Disk3D() ;                        ///< Destructor.
  
  // Accessors
  // ---------
 public:

#ifdef GYOTO_USE_CFITSIO
  /// Read parameters and arrays from FITS file.
  virtual void fitsRead(std::string filename_);

  /// Write parameters and arrays to FITS file.
  virtual void fitsWrite(std::string filename_);
#endif

  void file(std::string const &f);
  std::string file() const;
  void zsym(bool t);
  bool zsym() const;
  void tPattern(double t);
  double tPattern() const;
  void omegaPattern(double t);
  double omegaPattern() const;


  /// Set Disk3D::emissquant_.
  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setEmissquant(double * pattern);

  void opacity(double * pattern);

  /// Set Disk3D::velocity__.
  /**
   * The pointer is copied directly, not the array content.
   *
   * This is a low-level function. Beware that:
   *  - previously allocated array will not be freed automatically;
   *  - array attached when the destructor is called will be freed.
   */
  void setVelocity(double * pattern);

  /// Set Disk3D::emissquant_.
  /**
   * Disk3D::emissquant_ is freed if not NULL, reallocated, and
   * pattern is copied into emission_.
   *
   * If Disk3D::velocity_ or has been set previously with mismatching
   * sizes, it is deallocated too.
   *
   * Finally, Disk3D::nnu_, Disk3D::nphi_, Disk3D::nz_ and
   * Disk3D::nr_ are set according to naxes.
   *
   * \param pattern Array to copy as emission_. May be NULL in which
   * case emission_ is simply deallocated and set to NULL.
   *
   * \param naxes { nnu_, nphi_, nz_, nr_ }.
   */
  virtual void copyEmissquant(double const * const pattern = NULL,
			      size_t const naxes[4] = NULL);

  /// Get Disk3D::emissquant_.
  virtual double const * getEmissquant() const;

  /// Get { Disk3D::nnu_, Disk3D::nphi_, Disk3D::nz_, Disk3D::nr_ }.
  virtual void getEmissquantNaxes( size_t naxes[4] ) const ;

  virtual void copyOpacity(double const * const pattern = NULL,
			      size_t const naxes[4] = NULL);

  /// Get Disk3D::opacity_.
  virtual double const * opacity() const;

  /// Set Disk3D::velocity_.
  /**
   * Disk3D::velocity_ is first freed if not NULL and set to NULL.
   *
   * If pattern is not NULL, Disk3D::emissquant_ must have been set
   * previously with matching dimensions. Disk3D::velocity_ is
   * then reallocated, and pattern is copied into velocity_.
   *
   * \param pattern Array to copy as velocity_. May be NULL in which
   * case velocity_ is simply deallocated and set to NULL.
   *
   * \param naxes { nphi_, nz_, nr_ }.
   */
  virtual void copyVelocity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);
  /// Get Disk3D::velocity_.
  virtual double const * getVelocity() const;

  /// Set Disk3D::repeat_phi_.
  virtual void repeatPhi(size_t n);
  /// Get Disk3D::repeat_phi_.
  virtual size_t repeatPhi() const;

  /// Set Disk3D::nu0_.
  virtual void nu0(double freq);
  /// Get Disk3D::nu0_.
  virtual double nu0() const;

  /// Set Disk3D::dnu_.
  virtual void dnu(double dfreq);
  /// Get Disk3D::dnu_.
  virtual double dnu() const;

  /// Set Disk3D::rin_.
  void rin(double rrin);
  /// Get Disk3D::rin_.
  double rin() const;

  /// Set Disk3D::rout_.
  void rout(double rout);
  /// Get Disk3D::rout_.
  double rout() const;

  /// Set Disk3D::zmin_.
  void zmin(double zmin);
  /// Get Disk3D::zmin_.
  double zmin() const;

  /// Set Disk3D::zmax_.
  void zmax(double zmax);
  /// Get Disk3D::zmax_.
  double zmax() const;

  /// Set Disk3D::phimin_.
  void phimin(double phimin);
  /// Get Disk3D::phimin_.
  double phimin() const;

  /// Set Disk3D::phimax_.
  void phimax(double phimax);
  /// Get Disk3D::phimax_.
  double phimax() const;

 protected:
  void getIndices(size_t i[4], double const co[4], double nu=0.) const ;
  ///< Get emissquant_ cell corresponding to position co[4].

 public:
  int Impact(Photon *ph, size_t index, Astrobj::Properties *data);

  /// Get fluid 4-velocity at point.
  /**
   * Fill vel with the 4-vector velocity of the fluid at 4-position
   * pos.
   *
   * \param[in] pos 4-position at which to compute velocity;
   * \param[out] vel 4-velocity at pos.
   */
  virtual void getVelocity(double const pos[4], double vel[4])  ;

};

#endif
