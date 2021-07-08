/**
 * \file GyotoFlaredDiskSynchrotron.h
 * \brief A disk defined from a 2D grid in the equatorial plane
 * and extrapolated in the vertical direction with H/r<<1
 */

/*
    Copyright 2019-2021 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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

#ifndef __GyotoFlaredDiskSynchrotron_H_ 
#define __GyotoFlaredDiskSynchrotron_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class FlaredDiskSynchrotron; }
  class GridData2D;
}

#include <GyotoStandardAstrobj.h>
#include <GyotoGridData2D.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>

/**
 * \class Gyoto::Astrobj::FlaredDiskSynchrotron
 */

class Gyoto::Astrobj::FlaredDiskSynchrotron
: public Astrobj::Standard,
  public GridData2D,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::FlaredDiskSynchrotron>;
 private:
  SmartPointer<Spectrum::KappaDistributionSynchrotron> spectrumKappaSynch_;
  std::string filename_; ///< Optional FITS file name containing the arrays
  double hoverR_; ///< Value of aspect ratio H/R of flared disk, where R is the radius projected in the equatorial plane and H the altitude above the equatorial plane
  double numberDensityMax_cgs_; ///< Maximum cgs value of number density
  double temperatureMax_; ///< Maximum temperature in K
  double BMax_cgs_; ///< Maximun strenght of the 3 veceor magnetic field, defined by numberDensityMax_cgs_, temperatureMax_ and beta_
  double beta_;
  /**
   * An array of dimensionality double[nr_][nphi_][nt_]. In FITS
   * format, the first dimension is t, the second phi, and the third
   * r.
   */
  double * density_; ///< Surface density (&nu;, r, &phi;)
    /**
   * An array of dimensionality double[nr_][nphi_][nt_][2]. In FITS format,
   * the second dimension is phi, and the third r. The first plane in
   * the first FITS dimention is dr/dt, the second d&phi;/dt.
   */
  double * velocity_; ///< velocity(r, &phi;)
  double * Bvector_; ///<  4vector of the magnetic field
  double * time_array_; /// 1D Vector containing the times values of each time steps (dt not constant)
  double magnetizationParameter_; ///< (B<SUP>2</SUP>/(4 pi)) / (n<SUB>e</SUB> m<SUB>p</SUB> c<SUP>2</SUP>)
  double deltat_;///< time translation
  double gamm1_; /// polytropic index - 1
  bool flag_; /// flag for a fixed magnetic field or average

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  // Constructors - Destructor
  // -------------------------
  FlaredDiskSynchrotron(); ///< Standard constructor
  
  FlaredDiskSynchrotron(const FlaredDiskSynchrotron& ) ;///< Copy constructor
  virtual FlaredDiskSynchrotron* clone () const; ///< Cloner
  
  virtual ~FlaredDiskSynchrotron() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void file(std::string const &f) ;
  std::string file() const;
  void hoverR(double const hor) ;
  double hoverR() const;
  /*
    timeTranslation shifts the value of GridData2D::tmin_ and tmax_,
    allowing to scan the full simulation without having to change
    the value of the Screen observation time (which is typically
    not provided in M unit in the XML). 
    Choosing a negative timeTranslation, i.e. performing tmin_,tmax_-=dt, 
    amounts to increasing the Screen observation time by the same value, 
    tobs+=dt.

   */
  void timeTranslation_inMunit(double const dt) ;
  double timeTranslation_inMunit() const ;
  void magnetizationParameter(double rr);
  double magnetizationParameter() const;
  void kappaIndex(double index);
  double kappaIndex()const;
  double numberDensityMax() const;
  double numberDensityMax(std::string const &unit) const;
  void numberDensityMax(double dens) ;
  void numberDensityMax(double dens, std::string const &unit);
  void temperatureMax(double tt);
  double temperatureMax() const;
  void polytropicIndex(double gamma);
  double polytropicIndex() const;
  void betaAtMax(double beta);
  double betaAtMax() const;
  void copyDensity(double const *const density,
		   size_t const naxes[3]);
  double const * getDensity() const;
  void copyVelocity(double const *const velocity,
		    size_t const naxes[3]);
  double const * getVelocity() const;
  void copyBvector(double const *const Bvector,
        size_t const naxes[3]);
  double const * getBvector() const;
  void copyTimeArray(double const *const time_array, size_t const ntimes);
  double const * getTimeArray() const;
 public:
  using Generic::metric;
  std::vector<size_t> fitsRead(std::string filename) ;
  virtual double operator()(double const coord[4]) ;
  virtual void radiativeQ(double Inu[], 
			  double Taunu[],
			  double const nu_ems[], size_t nbnu, 
			  double dsem,
			  state_t const &coord_ph,
			  double const coord_obj[8]) const;
  virtual void getVelocity(double const pos[4], double vel[4]) ;



};

#endif
