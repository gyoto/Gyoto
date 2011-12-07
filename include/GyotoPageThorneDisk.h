/**
 * \file GyotoThinInfiniteDiskBL.h
 * \brief A geometrically thin, optically thick disk
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

/**
 * \class Gyoto::Astrobj::PageThorneDisk
 * \brief Geometrically thin disk in Kerr metric
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane, extending from r=r_ISCO to r=infinity.  The flux emitted
 *   at radius r is given by Page & Thorne (1974, ApJ 191:499,
 *   Eqs. 11b, 14, 15).
 *
 *   The metric must be either KerrBL or KerrKS.
 *
 *   Warning: The spin value and the inner radius are cached when the
 *   metric is assigned using setMetric(). If the spin is later
 *   changed in the metric, updateSpin() must be called.
 *
 */
class Gyoto::Astrobj::PageThorneDisk : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PageThorneDisk>;
 private:
  double aa_; ///< Spin
  double aa2_;
  double x0_;
  double x1_;
  double x2_;
  double x3_;

  // Constructors - Destructor
  // -------------------------
 public:
  
  PageThorneDisk(); ///< Standard constructor
  
  PageThorneDisk(const PageThorneDisk& ) ;///< Copy constructor
  virtual PageThorneDisk* clone () const; ///< Cloner
  
  virtual ~PageThorneDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual void setMetric(SmartPointer<Metric::Generic>);
  ///< Set metric, checking that it is either KerrBL or KerrKS

  virtual void updateSpin() ;
  ///< Get spin from metric, which must be KerrBL or KerrKS

 public:
  virtual double emission(double nu_em, double dsem,
			  double c_ph[8], double c_obj[8]) const;
 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ; ///< called from Factory
  static Astrobj::Subcontractor_t Subcontractor;
  static void Init();
#endif

};

#endif
