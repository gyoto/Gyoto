/**
 * \file GyotoThinDiskPL.h
 * \brief A subclass of ThinDisk emitting according to a powerlaw
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

#ifndef __GyotoThinDiskPL_H_ 
#define __GyotoThinDiskPL_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class ThinDiskPL; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>
#include <GyotoBlackBodySpectrum.h>

/**
 * \class Gyoto::Astrobj::ThinDiskPL
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
class Gyoto::Astrobj::ThinDiskPL : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDiskPL>;
 private:
  double PLSlope_;
  double PLRho_;
  double PLRadRef_;
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< disk black body

  // Constructors - Destructor
  // -------------------------
 public:
  
  ThinDiskPL(); ///< Standard constructor
  
  ThinDiskPL(const ThinDiskPL& ) ;///< Copy constructor
  virtual ThinDiskPL* clone () const; ///< Cloner
  
  virtual ~ThinDiskPL() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual double emission(double nu_em, double dsem,
			  double c_ph[8],double c_obj[8]) const;

  double emissionBB(double nu, 
		    double co[8]) const;

  int setParameter(std::string name, std::string content);
 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ; ///< called from Factory
#endif

};

#endif
