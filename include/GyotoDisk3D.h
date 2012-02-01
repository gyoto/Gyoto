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
 *   This class describes a thick disk.  The flux emitted
 *   at radius r altitude z and longitude phi at frequency nu is given 
 *   in a FITS file.
 *
 */
class Gyoto::Astrobj::Disk3D : public Gyoto::Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Disk3D>;
 private:
  std::string filename_;
  /**
   * An array of dimensionality double[nr_][nz_][nphi_][nnu_]. In FITS
   * format, the first dimension is nu, the second phi, the third
   * z and last r.
   */
  double * emissquant_; ///< Inu(nu, r, z, phi)

  /**
   * An array of dimensionality double[nr_][nz_][nphi_][3]. In FITS format,
   * the second dimension is phi, the third z and last r. The first plane in
   * the first FITS dimention is dphi/dt, the second dz/dt the last dr/dt.
   */
  double * velocity_; ///< velocity(r, z, phi)

  double dnu_;
  double nu0_;
  size_t nnu_;

  double dphi_;
  size_t nphi_;

  /**
   * XML elment: &lt;RepeatPhi&gt;.
   * FITS keyword: HIERARCH GYOTO Disk3D RepeatPhi
   */
  size_t repeat_phi_;
  //double phi0_==0, phi max is always 2*M_PI

  double dz_;
  double zmin_;
  size_t nz_;
  double zmax_;

  double dr_;
  double rin_;
  size_t nr_;
  double rout_;



  // Constructors - Destructor
  // -------------------------
 public:
  Disk3D(); ///< Standard constructor
  
  Disk3D(const Disk3D& ) ;///< Copy constructor
  virtual Disk3D* clone () const; ///< Cloner
  
  virtual ~Disk3D() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  //  virtual void   setInnerRadius(double); ///< Set rin_
  //  virtual void   setOuterRadius(double); ///< Set rout_

  virtual void fitsRead(std::string filename_);
  virtual void fitsWrite(std::string filename_);
  ///< Read data from file

  void setEmissquant(double * pattern);
  void setVelocity(double * pattern);
  /**
   * \param pattern: new emission_ array.
   * \param dims[4] = { nnu_, nphi_, nz_, nr_ };
   */
  virtual void copyEmissquant(double const * const pattern = NULL,
			      size_t const naxes[4] = NULL);
  ///< attach emissquant_ array and set its size
  virtual double const * const getEmissquant() const;
  virtual void getEmissquantNaxes( size_t naxes[4] ) const ;

  virtual void copyVelocity(double const * const pattern = NULL,
			      size_t const naxes[3] = NULL);
  virtual double const * const getVelocity() const;

  virtual void repeatPhi(size_t n);
  virtual size_t repeatPhi() const;

  virtual void nu0(double freq);
  virtual double nu0() const;

  virtual void dnu(double dfreq);
  virtual double dnu() const;

  void rin(double rrin);
  double rin();

  void dr(double dr);
  double dr();

  void rout(double rout);
  double rout();

  void zmin(double zmin);
  double zmin();

  void dz(double dz);
  double dz();

  void zmax(double zmax);
  double zmax();

  virtual int setParameter(std::string name, std::string content);

 protected:
  void getIndices(size_t i[4], double const co[4], double nu=0.) const ;
  ///< get emission_ cell corresponding to position co[4]

 public:
  int Impact(Photon *ph, size_t index, Astrobj::Properties *data);
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
