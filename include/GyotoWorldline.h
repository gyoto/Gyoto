/**
 * \file GyotoWorldline.h 
 * \brief Timelike or null geodesics
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
  class FactoryMessenger;
}

#include <GyotoSmartPointer.h>
#include <GyotoMetric.h>
#include <GyotoScreen.h>
#include <GyotoHooks.h>

/**
 * \class Gyoto::Worldline
 * \brief  Timelike or null geodesics
 *
 * Supported XML parameters:
 *  - InitialCoordinate or InitCoord: 8-element vector yielding the initial
 *    4-position and 4-velocity;
 *  - Only for massive particle (Gyoto::Astrobj::Star): Position
 *    (yielding initial 4-position) and Velocity (yielding initial
 *    3-velocity);
 *  - Delta: integration step, initial in case or adaptive step;
 *  - Adaptive or NonAdaptive: sets whether integration step should be
 *    adaptive; default: Adaptive.;
 *  - MaxIter: maximum number of iterations for the integration;
 *    default: 100000.
 * 
 */
class Gyoto::Worldline
:  protected Gyoto::Hook::Listener
{

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
  bool   adaptive_; ///< Whether integration should use adaptive delta
  bool secondary_; ///< choose 0 to compute only primary image
  double delta_;///< Initial integrating step ; defaults to 0.01
  double tmin_;///< Minimum time for integration, stop integration if t<tmin ; defaults to -DBL_MAX
  double * cst_; ///< Worldline's csts of motion (if any)
  size_t cst_n_; ///< Number of constants of motion
  int wait_pos_; ///< Hack in setParameters()
  double * init_vel_; ///< Hack in setParameters()
  size_t maxiter_ ; ///< Maximum number of iterations when integrating

  // Constructors - Destructor
  // -------------------------
 public: 
  Worldline() ; ///< Default constructor
  Worldline(const size_t sz) ; ///< Default constructor
  
  Worldline(const Worldline& ) ;                ///< Copy constructor
  
  /// Refine constructor
  /**
   * Meant to instanciate a copy of orig with a smaller step to refine
   * integration, for instance for more accurate radiative transfer
   * integration.
   *
   * See Photon::Photon(Photon* orig, size_t i0, int dir, double
   * step_max) and Photon::Refined.
   *
   * \param orig Worldline to refine
   * \param i0 Index of coordinate in orig to take as initial condition
   * \param dir Direction of integration
   * \param step_max Maximum integration step
   */
  Worldline(Worldline* orig, size_t i0, int dir, double step_max) ;

  virtual ~Worldline() ;                        ///< Destructor

  size_t getImin() const; ///< Get index of computed date furthest in the past
  size_t getImax() const; ///< Get index of computed date furthest in the future
  size_t getI0() const; ///< Get index of initial condition

  virtual double getMass() const = 0; ///< Get mass of particule.
  void   metric(SmartPointer<Metric::Generic>); ///< Set metric Smartpointer
  SmartPointer<Metric::Generic> metric() const; ///< Get metric
  virtual void   setInitCoord(const double coord[8], int dir = 0); ///< Set Initial coordinate

  /**
   * \brief Set initial coordinate
   *
   * \param pos initial 4-position
   * \param vel initial 3-velocity
   * \param dir direction of integration
   */
  virtual void setInitCoord(double pos[4], double vel[3], int dir=1);

  virtual void setPosition(double pos[4]); ///< Set initial 4-position
  virtual void setVelocity(double vel[3]); ///< Set initial 3-velocity

  void reset() ; ///< Forget integration, keeping initial contition
  void reInit() ; ///< Reset and recompute particle properties

  virtual std::string className() const ; ///< "Worldline"
  virtual std::string className_l() const ; ///< "worldline"

  // Memory management
  // ----------------- 
 protected:
  /**
   * The default size is GYOTO_DEFAULT_X_SIZE
   */
  virtual void xAllocate(); ///< Allocate x0, x1 etc. with default size

  /**
   * \param size : number of cells in each array x0, x1 etc.
   */
  virtual void xAllocate(size_t size); ///< Allocate x0, x1 etc. with a specified size.

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
  virtual size_t xExpand(int dir); ///< Expand x0, x1 etc... to hold more elements
 

  /**
   * If you need to expand more arrays than x0_ ... x3_ and the dots,
   * call this on your array before calling xExpand(int dir).
   *
   * \param[inout] x array to expand
   * \param[in] dir
   */
  virtual void xExpand(double * &x, int dir); ///< Expand one array to hold more elements

  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another Worldline
  void operator=(const Worldline&) ;        
  void delta(const double delta); ///< Set delta
  void delta(double, const std::string &unit);   ///< Set default step in specified units
  double delta() const ; ///< Get delta
  double delta(const std::string &unit) const ;  ///< Get default step in specified units
  double tMin() const ; ///< Get tmin value
  void tMin(double tlim); ///< Set tmin to a given value
  void adaptive (bool mode) ; ///< Set adaptive_
  bool adaptive () const ; ///< Get adaptive_
  void secondary (bool sec) ; ///< Set secondary_
  bool secondary () const ; ///< Get secondary_
  void maxiter (size_t miter) ; ///< Set maxiter_
  size_t maxiter () const ; ///< Get maxiter_

  /**
   * Return pointer to array holding the previously set
   * Metric-specific constants of motion
   */
  double const * getCst() const ; ///< Returns the worldline's cst of motion (if any)

  /// Set Metric-specific constants of motion
  /**
   * The will (re)allocate Worldline::cst_, copy cst into it, and set
   * Worldline::cst_n_.
   */
  void setCst(double const * cst, size_t const ncsts) ;

  /// Set or re-set the initial condition prior to integration.
  /**
   * \param gg    Gyoto::SmartPointer to the Gyoto::Metric in this universe;
   * \param coord 8 element array containing the initial condition,
   *        i.e. the 4-position and the 4-velocity of the Photon at
   *        the receiving end;
   * \param dir direction: 1 for future, -1 for past.
   */
  void setInitialCondition(SmartPointer<Metric::Generic> gg, 
			   const double coord[8],
			   const int dir) ;

  void getInitialCoord(double dest[8]) const; ///< Get initial coordinate
  void getCoord(size_t index, double dest[8]) const; ///< Get coordinates corresponding to index
  void getCartesianPos(size_t index, double dest[4]) const; ///< Get Cartesian expression of 4-position at index.


  virtual void xStore(size_t ind, double coord[8]) ; ///< Store coord at index ind
  virtual void xFill(double tlim) ; ///< Fill x0, x1... by integrating the Worldline from previously set inittial condition to time tlim

  /**
   * \brief Set parameter by name
   *
   * Assume MyKind is a subclass of Worldline which has two
   * members (a string StringMember and a double DoubleMember):
   * \code
   * int MyKind::setParameter(std::string name,
   *                          std::string content,
   *                          std::string unit) {
   *   if      (name=="StringMember") setStringMember(content);
   *   else if (name=="DoubleMember") setDoubleMember(atof(content.c_str()),
   *                                                  unit);
   *   else return Worldline::setParameter(name, content, unit);
   *   return 0;
   * }
   * \endcode
   *
   * \param name XML name of the parameter
   * \param content string representation of the value
   * \param unit string representation of the unit
   * \return 0 if this parameter is known, 1 if it is not.
   */
  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit) ;

#ifdef GYOTO_USE_XERCES
  /**
   * \brief Process XML entity
   * Uses wait_pos_ and init_vel_ to make sure setVelocity() is called
   * after setPosition().
   */
  virtual void setParameters(FactoryMessenger *fmp) ;
  /**
   * Derived classes implementations should implement fillElement to save their
   * parameters to XML and call the generic implementation to save
   * generic parts such as adaptive_: Worldline::fillElement(fmp).
   */
  virtual void fillElement(FactoryMessenger *fmp) const ;
                                             ///< XML output
#endif

  // Accessors
  // ---------
 public:
  /**
   * \brief Get number of computed dates
   */
  size_t get_nelements() const;

  /**
   * \brief Get computed dates
   */
  void get_t(double *dest) const;

  
  /// Get the 6 Cartesian coordinates for specific dates.
  /**
   * The 6 coordinates (x, y, z, dx/dt, dy/dt, dz/dt) will be computed
   * using the integrator and interpolated if necessary, so they will
   * be as accurate as possible. Transforming to Cartesian coordinates
   * is not necessarily meaningful.
   *
   * \param[in] dates List of dates for which the coordinates are to
   *                be computed;
   *
   * \param[in] n_dates Number of dates to compute ;
   *
   * \param[out] x, y, z, xprime, yprime, zprime Arrays in which to
   * store the result. These pointer may be set to NULL to retrieve
   * only part of the information. Else, they must be pre-allocated.
   *
   */
  void getCartesian(double const * const dates, size_t const n_dates,
		double * const x, double * const y,
		double * const z, double * const xprime=NULL,
		double * const yprime=NULL,  double * const zprime=NULL) ;

  /**
   * \brief Get 3-position in cartesian coordinates for computed dates
   */
  void get_xyz(double* x, double *y, double *z) const;

  /**
   * \brief Get 8-coordinates for specific dates.
   *
   * The coordinates will be
   * computed using the integrator, so they will be as accurate as
   * possible. Some heuristics are used to speed up the process and it
   * is presumably faster to call this routine with a sorted list of
   * dates. The line will be integrated further as required. An error
   * will be thrown if it is not possible to reach a certain date.
   *
   * \param dates the list of dates for which the coordinates are to
   *                be computed;
   * \param n_dates the number of dates to compute ;
   * \param x1dest, x2dest, x3dest, x0dot, x1dot, x2dot, x3dot arrays
   *               in which to store the result. These pointer may be
   *               set to NULL to retrieve only part of the
   *               information. They must be pre-allocated.
   *
   */
  void getCoord(double const * const dates, size_t const n_dates,
		double * const x1dest,
		double * const x2dest, double * const x3dest,
		double * const x0dot=NULL,  double * const x1dot=NULL,
		double * const x2dot=NULL,  double * const x3dot=NULL) ;

  /**
   * \brief Get all computed positions
   *
   *  Get all the pre-computed 8 coordinates (e.g. thanks to a prior
   *  call to xFill()) of this worldline.
   */
  void getCoord(double *x0, double *x1, double *x2, double *x3) const ;

  /**
   * \brief Bring &theta; in [0,&Pi;] and &phi; in [0,2&Pi;]
   *
   * checkPhiTheta() Modifies coord if the corrdinates are spherical-like
   * so that coord[2]=theta is in [0,pi] and coord[3]=phi is in [0,2pi].
   * Important to use in all astrobj in spherical coordinates
   * to prevent "z-axis problems".
   */
  void checkPhiTheta(double coord[8]) const;

  /**
   * \brief Get computed positions in sky coordinates
   */
  void getSkyPos(SmartPointer<Screen> screen, double *dalpha, double *ddellta, double *dD) const;

  /**
   * \brief Get computed 4-velocities
   */
  void get_dot(double *x0dot, double *x1dot, double *x2dot, double *x3dot) const ;

  /**
   * \brief Get computed 3-velocities
   */
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
  
 protected:
  virtual void tell(Gyoto::Hook::Teller*);
  class IntegState;
};


/**
 * \class Gyoto::Worldline::IntegState
 * \brief Current state of a geodesic integration
 */

class Gyoto::Worldline::IntegState : SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Worldline::IntegState>;

 private:
  /// Worldline that we are integrating.
  /**
   * Beware this is not a SmartPointer. Make sure line_ still exists
   * when calling nestStep().
   */
  Worldline * line_;

  /// The Metric in this end of the Universe.
  /**
   * Taken from Worldline::line_, never updated.
   */
  Gyoto::SmartPointer<Gyoto::Metric::Generic> gg_;
  
  double coord_[8]; ///< Previously determined coordinate.
  double norm_; ///< Current norm of the 4-velocity.
  double normref_; ///< Initial norm of the 4-velocity.
  double delta_; ///< Integration step (current in case of adaptive).

  /// Whether Worldline::delta_ is adaptive.
  /**
   * Taken from Worldline::line_, never updated.
   */
  bool adaptive_;

 public:
  /// Constructor
  /**
   * \param line The Worldline that we are integrating. Sets:
   * Worldline::line_, Worldline::gg_, Worldline::adaptive_.
   * \param coord Initial coordinate.
   * \param delta Integration step. Sign determines direction.
   */
  IntegState(Worldline * line, const double *coord, const double delta);
  
  /// Make one step.
  /**
   * \param[out] coord Next position-velocity;
   * \param[in] h1max maximum step in case of adaptive integration
   */
  virtual int nextStep(double *coord, double h1max=1e6);

  virtual ~IntegState();
};

#endif
