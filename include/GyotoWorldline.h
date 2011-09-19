/**
 * \file GyotoWorldline.h 
 * \brief geodesics ?
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

#ifndef __GyotoWorldline_H_
#define __GyotoWorldline_H_ 

#include <iostream>
#include <fstream>
#include <string>
#include <GyotoDefs.h>

namespace Gyoto {
  class Worldline;
}

#include <GyotoSmartPointer.h>
#include <GyotoMetric.h>
#include <GyotoWorldlineIntegState.h>
#include <GyotoScreen.h>

/**
 * \class Gyoto::Worldline
 * \brief geodesic?
 * 
 */
class Gyoto::Worldline {

  // Data : 
  // -----
 protected:
  SmartPointer<Gyoto::Metric::Generic> metric_ ; ///< The Gyoto::Metric in this part of the universe
  double* x0_;///< t or T
  double* x1_;///< r or x
  double* x2_;///< theta or y
  double* x3_;///< phi or z
  double* x0dot_;///< tdot or Tdot
  double* x1dot_;///< rdot or xdot
  double* x2dot_;///< thetadot or ydot
  double* x3dot_;///< phidot or zdot
  size_t x_size_;///< Size of x0, x1... arrays
  size_t imin_;///< Minimum index for which x0, x1... have been computed
  size_t i0_;  ///< Index of initial condition in array
  size_t imax_;///< Maximum index for which x0, x1... have been computed
  double delta_;///< Initial integrating step ; defaults to 0.01
  double tlim_;///< Minimum time for integration, stop integration if t<tlim ; defaults to 0
  double * cst_; ///< Worldline's csts of motion (if any)
  size_t cst_n_; ///< Number of constants of motion

  // Constructors - Destructor
  // -------------------------
 public: 
  Worldline() ; ///< Default constructor
  Worldline(const size_t sz) ; ///< Default constructor
  
  Worldline(const Worldline& ) ;                ///< Copy constructor
  
  /// Constructor from a file (see \c sauve(FILE*) )
  //Worldline(FILE *) ;                    
  
  virtual ~Worldline() ;                        ///< Destructor

  int getImin() const;
  int getImax() const;
  int getI0() const;

  virtual double getMass() const = 0; ///< Get mass of particule.
  void   setMetric(SmartPointer<Metric::Generic>); ///< Set metric Smartpointer
  SmartPointer<Metric::Generic> getMetric() const; ///< Get metric
  void   setInitCoord(const double coord[8], const int dir = 0); ///< Set Initial coordinate
  void reset() ; ///< Forget integration, keeping initial contition

  virtual std::string className() const ; ///< "Worldline"
  virtual std::string className_l() const ; ///< "worldline"

  // Memory management
  // ----------------- 
 protected:
  /**
   * The default size is GYOTO_DEFAULT_X_SIZE
   */
  void xAllocate(); ///< Allocate x0, x1 etc. with default size

  /**
   * \param size : number of cells in each array x0, x1 etc.
   */
  void xAllocate(size_t size); ///< Allocate x0, x1 etc. with a specified size.

  /**
   * Double the size of arrays x0, x1 etc. and copy old version of the
   * array in the first half if dir =1 and in the second half if dir
   * =-1.
   *
   * \param dir : 1 to expand after last element, -1 to expand before
   * first element
   *
   * \return ind : if dir=1, new index of old last element, if dir=-1,
   * new index of old first element
   */
  size_t xExpand(int dir); ///< Expand x0, x1 etc... to hold more elements
 
  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another Worldline
  void operator=(const Worldline&) ;        
  void setDelta(const double delta); ///< Set delta
  double getTlim(); ///< Get tlim value
  void setTlim(double tlim); ///< Set tlim to a given value

  /**
   * Return pointer to array holding the previously set
   * Metric-specific constants of motion
   */
  double const * getCst() const ; ///< Returns the worldline's cst of motion (if any)
  /**
   * Set Metric-specific constants of motion
   */
  void setCst(double const * cst, size_t const ncsts) ;

  /**
   * Set initial condition for this Photon :
   *
   * \param gg : Gyoto::SmartPointer to the Gyoto::Metric in this universe;
   *
   * \param coord : 8 element array containing the initial condition,
   *        i.e. the 4-position and the 4-velocity of the Photon at
   *        the receiving end;
   *
   * \param sys : an integer stating in which coordinate system coord
   *        is given.
   */
  void setInitialCondition(SmartPointer<Metric::Generic> gg, 
			   const double coord[8],
			   const int dir) ;
  ///<Set or re-set the initial condition prior to integration.
 

  void getInitialCoord(double dest[8]) const; ///< get initial coordinate
  void getCoord(size_t index, double dest[8]) const; ///< get coordinates corresponding to index
  void getCartesianPos(size_t index, double dest[4]) const; ///< get Cartesian expression of 4-position at index.

  void xFill(double tlim) ; ///< Fill x0, x1... by integrating the Worldline from previously set inittial condition to time tlim

  // Accessors
  // ---------
 public:
  //virtual void position(double t, double* res) = 0 ;
  int get_nelements() const;
  void get_t(double *dest) const;

  /**
   * Get the 6 Cartesian coordinates (x, y, z, dx/dt, dy/dt, dz/dt)
   * for specific dates. The coordinates will be computed using the
   * integrator, so they will be as accurate as possible. Transforming
   * to Cartesian coordinates is not necessarily meaningful.
   *
   * \param dates: the list of dates for which the coordinates are to
   *                be computed;
   *
   * \param n_dates: the number of dates to compute ;
   *
   * \param x*: arrays in which to store the result. These pointer may
   *               be set to NULL to retrieve only part of the
   *               information. They must be pre-allocated.
   *
   */
  void getCartesian(double const * const dates, size_t const n_dates,
		double * const x, double * const y,
		double * const z, double * const xprime=NULL,
		double * const yprime=NULL,  double * const zprime=NULL) ;

  void get_xyz(double* x, double *y, double *z) const;

  /**
   * Get 8-coordinates for spcific dates. The coordinates will be
   * computed using the integrator, so they will be as accurate as
   * possible. Some heuristics are used to speed up the process and it
   * is presumably faster to call this routine with a sorted list of
   * dates. The line will be integrated further as required. An error
   * will be thrown if it is not possible to reach a certain date.
   *
   * \param dates: the list of dates for which the coordinates are to
   *                be computed;
   *
   * \param n_dates: the number of dates to compute ;
   *
   * \param x*: arrays in which to store the result. These pointer may
   *               be set to NULL to retrieve only part of the
   *               information. They must be pre-allocated.
   *
   */
  void getCoord(double const * const dates, size_t const n_dates,
		double * const x1dest,
		double * const x2dest, double * const x3dest,
		double * const x0dot=NULL,  double * const x1dot=NULL,
		double * const x2dot=NULL,  double * const x3dot=NULL) ;

  /**
   *  Get all the pre-computed 8 coordinates (e.g. thanks to a prior
   *  call to xFill()) of this worldline.
   */
  void getCoord(double *x0, double *x1, double *x2, double *x3) const ;
  void getSkyPos(SmartPointer<Screen> screen, double *dalpha, double *ddellta, double *dD) const;
  void get_dot(double *x0dot, double *x1dot, double *x2dot, double *x3dot) const ;
  void get_prime(double *x1prime, double *x2prime, double *x3prime) const ;
  
  // Outputs
  // -------
 public:
  //virtual void sauve(FILE *) const ;            ///< Save in a file
  void save_txyz(char * fichierxyz) const ;            ///< Save in a file
  void save_txyz(char* const filename, double const t1, double const  mass_sun,
		 double const distance_kpc, std::string const unit, SmartPointer<Screen> sc = NULL);///< Save, converted

  /// Display
  friend std::ostream& operator<<(std::ostream& , const Worldline& ) ;
  
  
};

#endif
