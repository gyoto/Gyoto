/**
 * \file GyotoPatternDisk.h
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
  std::string filename_;
  /**
   * An array of dimensionality double[nr_][nphi_][nnu_]. In FITS
   * format, the first dimension is nu, the second phi, and the third
   * r.
   */
  double * emission_; ///< Inu(nu, r, phi)

  double * opacity_; ///< same dimenstions as emission, or NULL

  /**
   * An array of dimensionality double[nr_][nphi_][2]. In FITS format,
   * the second dimension is phi, and the third r. The first plane in
   * the first FITS dimention is dphi/dt, the second dr/dt.
   */
  double * velocity_; ///< velocity(, r, phi)

  /**
   * In case of adaptive grid.
   */
  double * radius_; ///< radius vector

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

  double dnu_;
  double nu0_;
  size_t nnu_;

  double dphi_;
  size_t nphi_;

  /**
   * XML elment: &lt;RepeatPhi&gt;.
   * FITS keyword: HIERARCH GYOTO PatternDisk RepeatPhi
   */
  size_t repeat_phi_;
  //double phi0_==0, phi max is always 2*M_PI

  double dr_;
  size_t nr_;
  //  double r0_; // this is rin_



  // Constructors - Destructor
  // -------------------------
 public:
  PatternDisk(); ///< Standard constructor
  
  PatternDisk(const PatternDisk& ) ;///< Copy constructor
  virtual PatternDisk* clone () const; ///< Cloner
  
  virtual ~PatternDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual void   setInnerRadius(double); ///< Set rin_
  virtual void   setOuterRadius(double); ///< Set rout_

  /**
   * Unit: radians per geometrical unit time.
   */
  virtual void   setPatternVelocity(double); ///< Set pattern angular velocity
  virtual double getPatternVelocity(); ///< Set pattern angular velocity
  virtual void fitsRead(std::string filename_);
  virtual void fitsWrite(std::string filename_);
  ///< Read data from file


  /**
   * \param pattern: new emission_ array.
   * \param dims[3] = { nnu_, nphi_, nr_ };
   */
  virtual void copyIntensity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);
  ///< attach emission_ array and set its size
  virtual double const * const getIntensity() const;
  virtual void getIntensityNaxes( size_t naxes[3] ) const ;

  virtual void copyOpacity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);
  virtual double const * const getOpacity() const;

  virtual void copyVelocity(double const * const pattern = NULL,
			      size_t const naxes[2] = NULL);
  virtual double const * const getVelocity() const;

  virtual void copyGridRadius(double const * const pattern = NULL,
			      size_t nr = 0 );
  virtual double const * const getGridRadius() const;

  virtual void repeatPhi(size_t n);
  virtual size_t repeatPhi() const;

  virtual void nu0(double freq);
  virtual double nu0() const;

  virtual void dnu(double dfreq);
  virtual double dnu() const;

  virtual int setParameter(std::string name, std::string content);

 protected:
  void getIndices(size_t i[3], double const co[4], double nu=0.) const ;
  ///< get emission_ cell corresponding to position co[4]

 public:
  virtual double emission(double nu_em, double dsem,
			  double c_ph[8], double c_obj[8]) const;
  virtual double transmission(double nu_em, double dsem, double coord[8]) const;

  virtual void getVelocity(double const pos[4], double vel[4])  ;

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
  virtual void setParameters(FactoryMessenger *fmp);
#endif

};

#endif
