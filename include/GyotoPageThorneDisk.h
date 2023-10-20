/**
 * \file GyotoPageThorneDisk.h
 * \brief A geometrically thin, optically thick disk
 *
 *  This class describes a disk contained in the z=0 (equatorial)
 *  plane, extending from r=r_ISCO to r=infinity.  The flux emitted
 *  at radius r is given by Page & Thorne (1974, ApJ 191:499,
 *  Eqs. 11b, 14, 15).
 *
 *  Only bolometric flux is implemented (as quantity User4), no
 *  spectral resolution.
 */

/*
  Copyright 2011-2016, 2018 Frederic Vincent, Thibaut Paumard
  
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

#ifndef __GyotoPageThorneDisk_H_ 
#define __GyotoPageThorneDisk_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class PageThorneDisk; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>
#include <GyotoBlackBodySpectrum.h>

/**
 * \class Gyoto::Astrobj::PageThorneDisk
 * \brief Geometrically thin disk in Kerr metric
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane, extending from r=r_ISCO to r=infinity.  The flux emitted
 *   at radius r is given by Page & Thorne (1974, ApJ 191:499,
 *   Eqs. 11b, 14, 15).
 *
 *   The metric must be either KerrBL or KerrKS. Emission, Spectrum
 *   and BinSpectrum are <STRONG>not</STRONG> provide, the only
 *   intensity provided is provided, as quantity User4 and it is the
 *   default quantity returned if nothing is requested. The other
 *   quantities implemented in ThinDisk are also provided.
 *
 */
class Gyoto::Astrobj::PageThorneDisk
: public Astrobj::ThinDisk,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PageThorneDisk>;
 private:
  double aa_; ///< Generic::gg_ spin parameter, monitored by tell()
  double aa2_; ///< aa_<SUP>2</SUP>
  double x0_; ///< Value cached for bolometricEmission()
  double x1_; ///< Value cached for bolometricEmission()
  double x2_; ///< Value cached for bolometricEmission()
  double x3_; ///< Value cached for bolometricEmission()
  double mdot_; ///< accretion rate (for BB spectrum computation)
  bool uniflux_; ///< Flag for uniform flux = 1
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< disk black body
  ///< emission law
  
  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

  PageThorneDisk(); ///< Standard constructor
  
  PageThorneDisk(const PageThorneDisk& ) ;///< Copy constructor
  virtual PageThorneDisk* clone () const; ///< Cloner
  
  virtual ~PageThorneDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  using ThinDisk::metric;
  virtual void metric(SmartPointer<Metric::Generic>);
  ///< Set metric, checking that it is either KerrBL or KerrKS

  /// Set #mdot_ to v
  void mdot(double v);
  double mdot() const;
  void uniFlux(bool t) ;
  bool uniFlux() const ;

 private:
  virtual void updateSpin() ;
  ///< Get spin from metric, which must be KerrBL or KerrKS

 public:
  using ThinDisk::emission;
  /**
   * \brief Not implemented
   * Throws a Gyoto::Error
   */
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;

  /**
   * \brief Bolometric emission
   *
   * Similar to Generic::emission(), but bolometric.
   */
  virtual double bolometricEmission(double nuem, double dsem,
				    double const c_obj[8]) const;

  /**
   * \brief 
   * processHitQuantities fills the requested data in Impact. For
   * PageThorneDisk, only fill User4, which corresponds to bolometric
   * intensity.
   */
  virtual void processHitQuantities(Photon* ph, state_t const &coord_ph_hit,
                                   double const *coord_obj_hit, double dt,
                                   Astrobj::Properties* data) const;

  Gyoto::Quantity_t getDefaultQuantities();

  // Hook::Listener API //
 public:
  /**
   * \brief Update PageThorneDisk::aa_
   *
   * Calls updateSpin().
   *
   * See Hook::Listener::tell()
   */
  virtual void tell(Gyoto::Hook::Teller *msg);

  virtual void radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
       Eigen::Matrix4d *Onu,
       double const *nuem , size_t nbnu, double dsem,
       state_t const &cph, double const *co) const;

};

#endif
