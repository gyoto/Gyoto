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

#ifdef HAVE_BOOST
# include <functional>
# include <array>
# include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#endif

namespace Gyoto {
  class Worldline;
  class FactoryMessenger;
}

#include <GyotoSmartPointer.h>
#include <GyotoMetric.h>
#include <GyotoScreen.h>
#include <GyotoHooks.h>

/// Define the bunch of Properties that make up the Worldline interface
/**
 * This macro, which is called automatically by
 * GYOTO_WORLDLINE_PROPERTY_END(c, a), must be inserted in the
 * definition of the Property list for any class derived from
 * Worldline.
 */
#define GYOTO_WORLDLINE_PROPERTIES(c)					\
  GYOTO_PROPERTY_BOOL(c, HighOrderImages, PrimaryOnly, _secondary)	\
  GYOTO_PROPERTY_DOUBLE(c, RelTol, _relTol)				\
  GYOTO_PROPERTY_DOUBLE(c, AbsTol, _absTol)				\
  GYOTO_PROPERTY_DOUBLE(c, DeltaMaxOverR, _deltaMaxOverR)		\
  GYOTO_PROPERTY_DOUBLE(c, DeltaMax, _deltaMax)				\
  GYOTO_PROPERTY_DOUBLE(c, DeltaMin, _deltaMin)				\
  GYOTO_PROPERTY_STRING(c, Integrator, _integrator)			\
  GYOTO_PROPERTY_SIZE_T(c, MaxIter, _maxiter)				\
  GYOTO_PROPERTY_BOOL(c, Adaptive, NonAdaptive, _adaptive)		\
  GYOTO_PROPERTY_DOUBLE_UNIT(c, Delta, _delta)				\
  GYOTO_PROPERTY_VECTOR_DOUBLE(c, InitCoord, _initCoord)		\
  GYOTO_PROPERTY_METRIC(c, Metric, _metric)

/// Define the wrapper accessors used in GYOTO_WORLDLINE_PROPERTIES(class)
/**
   This macro, which is called automatically by
   GYOTO_WORLDLINE_PROPERTY_END(c, a), must be called once with the
   definition of the methods (.C file) of any class that derives from
   Worldline. The corresponding macro GYOTO_WORLDLINE must be called
   in the corresponding class declaration (.h file).

   This is made necessary by how multiple inheritence works: directly
   using the accessors in the Worldline API leads to segfault at
   runtime (unless too much extra care is taken) and may go unnoticed.

   These accessors must be declared in the class declaration using the
   GYOTO_WORLDLINE macro.
*/
#define GYOTO_WORLDLINE_ACCESSORS(c)					\
  void c::_secondary(bool s) {secondary(s);}				\
  bool c::_secondary() const {return secondary();}			\
  void c::_adaptive(bool s) {adaptive(s);}				\
  bool c::_adaptive() const {return adaptive();}			\
  void c::_relTol(double f){relTol(f);}					\
  double c::_relTol()const{return relTol();}				\
  void c::_absTol(double f){absTol(f);}					\
  double c::_absTol()const{return absTol();}				\
  void c::_deltaMin(double f){deltaMin(f);}				\
  double c::_deltaMin()const{return deltaMin();}			\
  void c::_deltaMax(double f){deltaMax(f);}				\
  double c::_deltaMax()const{return deltaMax();}			\
  void c::_deltaMaxOverR(double f){deltaMaxOverR(f);}			\
  double c::_deltaMaxOverR()const{return deltaMaxOverR();}		\
  void c::_delta(double f){delta(f);}					\
  double c::_delta()const{return delta();}				\
  void c::_delta(double f, std::string const &u){delta(f, u);}		\
  double c::_delta(std::string const &u)const{return delta(u);}		\
  void c::_maxiter(size_t f){maxiter(f);}				\
  size_t c::_maxiter()const{return maxiter();}				\
  void c::_integrator(std::string const &f){integrator(f);}		\
  std::string c::_integrator() const {return integrator();}		\
  std::vector<double> c::_initCoord()const{return initCoord();}		\
  void c::_initCoord(std::vector<double> const &f){initCoord(f);}	\
  void c::_metric(SmartPointer<Metric::Generic>f){metric(f);}		\
  SmartPointer<Metric::Generic> c::_metric() const{return metric();}

/// Drop-in replacement for GYOTO_PROPERTY_END(), which adds the Worldline interface
/**
 * This macro replaces GYOTO_PROPERTY_END(c, a) for classes that
 * derive from Worldline. It calls GYOTO_WORLDLINE_PROPERTIES(a) and
 * GYOTO_WORLDLINE_ACCESSORS(c). If this macro is used,
 * GYOTO_WORLDLINE must be called in the class declaration (.h file).
 */
#define GYOTO_WORLDLINE_PROPERTY_END(c, a) \
  GYOTO_WORLDLINE_PROPERTIES(c)		   \
  GYOTO_PROPERTY_END(c, a)		   \
  GYOTO_WORLDLINE_ACCESSORS(c)

/// Declare the Worldline interface wrappers
/**
   This macro must be called in the class declaration (.h file), in a
   public section. Its sibling GYOTO_WORLDLINE_ACCESSORS(c) must be
   called with the class method definition (.C file). Note that
   GYOTO_WORLDLINE_PROPERTY_END(c, a) calls
   GYOTO_WORLDLINE_ACCESSORS(c).
*/
#define GYOTO_WORLDLINE					\
  void _delta(const double delta);			\
  void _delta(double, const std::string &unit);		\
  double _delta() const ;				\
  double _delta(const std::string &unit) const ;	\
  void _adaptive (bool mode) ;				\
  bool _adaptive () const ;				\
  void _secondary (bool sec) ;				\
  bool _secondary () const ;				\
  void _maxiter (size_t miter) ;			\
  size_t _maxiter () const ;				\
  void _integrator(std::string const & type);		\
  std::string _integrator() const ;			\
  double _deltaMin() const;				\
  void _deltaMin(double h1);				\
  void _absTol(double);					\
  double _absTol()const;				\
  void _relTol(double);					\
  double _relTol()const;				\
  void _deltaMax(double h1);				\
  double _deltaMax()const;				\
  double _deltaMaxOverR() const;			\
  void _deltaMaxOverR(double t);			\
  std::vector<double> _initCoord()const;		\
  void _initCoord(std::vector<double> const&f);		\
  void _metric(SmartPointer<Metric::Generic>);		\
  SmartPointer<Metric::Generic> _metric() const;

/**
 * \class Gyoto::Worldline
 * \brief  Timelike or null geodesics
 *
 * Their are two derived classes: Photon and Star. A Worldline can be
 * integrated from an initial condition either backward or forward in
 * time using xFill() (Photon::hit() also integrates the
 * Worldline). Member #state_ holds the integration state as well as
 * an integrator. There are several kinds of integration states, that
 * derive from IntegState::Generic.
 *
 * The coordinates of the Worldline are stored in #x0_, #x1_, #x2_,
 * #x3_, #x0dot_, #x1dot_, #x2dot_ ans #x3dot_. Those arrays are
 * extended as needed using xExpand(). These coordinates can be
 * retrieved using get_t(), get_xyz(), getCartesian(), getCoord() etc.
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
 *  - DeltaMin, DeltaMax, DeltaMaxOverR, AbsTol, RelTol: tuning
 *    parameters used by some of the integrators.
 * 
 */
class Gyoto::Worldline
:  protected Gyoto::Hook::Listener
{

  // Data : 
  // -----
 public:
  int stopcond; ///< Whether and why integration is finished

 protected:
  SmartPointer<Gyoto::Metric::Generic> metric_ ; ///< The Gyoto::Metric in this part of the universe
  double* x0_;///< t or T
  double* x1_;///< r or x
  double* x2_;///< &theta; or y
  double* x3_;///< &phi; or z
  double* x0dot_;///< tdot or Tdot
  double* x1dot_;///< rdot or xdot
  double* x2dot_;///< &theta;dot or ydot
  double* x3dot_;///< &phi;dot or zdot
  size_t x_size_;///< Size of #x0_, #x1_... arrays
  size_t imin_;///< Minimum index for which #x0_, #x1_... have been computed
  size_t i0_;  ///< Index of initial condition in array
  size_t imax_;///< Maximum index for which #x0_, #x1_... have been computed
  bool   adaptive_; ///< Whether integration should use adaptive delta

  /**
   * \brief Experimental: choose 0 to compute only primary image
   *
   * This feature is in development.
   */
  bool secondary_;

  /**
   * \brief Initial integrating step
   *
   * Default: #GYOTO_DEFAULT_DELTA
   */
  double delta_;


  /**
   * \brief Time limit for the integration (geometrical units)
   *
   * Computation does not go back before #tmin_. Default is -DBL_MAX. #tmin_ is
   * always expressed in geometrical units, it is essentially a tuning
   * parameter for the ray-tracing process. #tmin_ should be chosen to
   * always be longer than the distance between the screen and the
   * object.
   */
  double tmin_;

  double * cst_; ///< Worldline's csts of motion (if any)
  size_t cst_n_; ///< Number of constants of motion
  int wait_pos_; ///< Hack in setParameters()
  double * init_vel_; ///< Hack in setParameters()
  size_t maxiter_ ; ///< Maximum number of iterations when integrating

  /**
   * \brief Minimum integration step for the adaptive integrator
   *
   * The default (#GYOTO_DEFAULT_DELTA_MIN) is usually fine.
   *
   * For IntegState::Legacy, set it in the Metric instead!
   */
  double delta_min_;

  /**
   * \brief Maximum integration step for the adaptive integrator
   *
   * The default (#GYOTO_DEFAULT_DELTA_MAX) is usually fine.
   *
   * For IntegState::Legacy, set it in the Metric instead!
   */
  double delta_max_;

  /**
   * \brief Numerical tuning parameter
   *
   * For IntegState::Legacy, set it in the Metric instead!
   *
   * Ensure that delta (the numerical integration step) is never
   * larger than a fraction of the distance between the current
   * location and the center of the coordinate system.
   *
   * The default (#GYOTO_DEFAULT_DELTA_MAX_OVER_R) is usually fine.
   */
  double delta_max_over_r_;

  /**
   * \brief Absolute tolerance of the integrator
   *
   * Used by the adaptive integrators implemented in
   * IntegState::Boost. Refer to the boost::numeric::odeint
   * documentation for more details.
   */
  double abstol_;

  /**
   * \brief Absolute tolerance of the integrator
   *
   * Used by the adaptive integrators implemented in
   * IntegState::Boost. Refer to the boost::numeric::odeint
   * documentation for more details.
   */
  double reltol_;

  // Constructors - Destructor
  // -------------------------
 public: 
  Worldline() ; ///< Default constructor
  
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

  size_t getImin() const; ///< Get #imin_
  size_t getImax() const; ///< Get #imax_
  size_t getI0() const; ///< Get #i0_

  virtual double getMass() const = 0; ///< Get mass of particule.
  void   metric(SmartPointer<Metric::Generic>); ///< Set metric Smartpointer
  SmartPointer<Metric::Generic> metric() const; ///< Get metric

  void initCoord(std::vector<double> const&); 
  std::vector<double> initCoord() const; 

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

  /**
   * \brief Set the integrator
   *
   * Initialize #state_ to use the required integrator.
   *
   * \param[in] type Either "Legacy" or (if HAVE_BOOST) one of
   *                 "runge_kutta_cash_karp54",
   *                 "runge_kutta_fehlberg78", "runge_kutta_dopri5",
   *                 "runge_kutta_cash_karp54_classic"
   */
  void integrator(std::string const & type);

  /**
   * \brief Describe the integrator used by #state_
   */
  std::string integrator() const ;

  /**
   * \brief Get #delta_min_
   */
  double deltaMin() const;

  /**
   * \brief Set #delta_min_
   */
  void deltaMin(double h1);

  /**
   * \brief Get #delta_max_
   */
  double deltaMax() const;


  void absTol(double); ///< Set #abstol_
  double absTol()const; ///< Get #abstol_
  void relTol(double); ///< Set #reltol_
  double relTol()const; ///< Get #reltol_

  /**
   * Get delta max at a given position
   *
   * \param[in] pos 4-position
   * \param[in] delta_max_external external constraint on delta_max
   * \return the smallest value between #delta_max_,
   * delta_max_external, and R*#delta_max_over_r_ where R is pos[1] in
   * spherical coordinates and max(x1, x2, x3) in Cartesian
   * coordinates.
   */
  virtual double deltaMax(double const pos[8], double delta_max_external) const;

  /**
   * Set delta_max_
   */
  void deltaMax(double h1);

  double deltaMaxOverR() const; ///< Get #delta_max_over_r_
  void deltaMaxOverR(double t); ///< Set #delta_max_over_r_

  // Memory management
  // ----------------- 
 protected:
  /**
   * The default size is #GYOTO_DEFAULT_X_SIZE
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
  void delta(const double delta); ///< Set #delta_
  void delta(double, const std::string &unit);   ///< Set #delta_ in specified units
  double delta() const ; ///< Get #delta_
  double delta(const std::string &unit) const ;  ///< Get #delta_ in specified units
  double tMin() const ; ///< Get #tmin_
  double tMin(const std::string &unit) const ;  ///< Get #tmin_ in specified unit
  void tMin(double tlim); ///< Set #tmin_
  void tMin(double, const std::string &unit);   ///< Set #tmin_ in specified unit
  void adaptive (bool mode) ; ///< Set #adaptive_
  bool adaptive () const ; ///< Get #adaptive_
  void secondary (bool sec) ; ///< Set #secondary_
  bool secondary () const ; ///< Get #secondary_
  void maxiter (size_t miter) ; ///< Set #maxiter_
  size_t maxiter () const ; ///< Get #maxiter_

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

  class IntegState {
  public:
    class Generic;
    class Legacy;
#ifdef HAVE_BOOST
    class Boost;
#endif
  };


  /**
   * \brief An object to hold the integration state
   */
  SmartPointer<Worldline::IntegState::Generic> state_;
};


/**
 * \class Gyoto::Worldline::IntegState::Generic
 * \brief Current state of a geodesic integration
 */
class Gyoto::Worldline::IntegState::Generic : SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Worldline::IntegState::Generic>;
 protected:
  /// Worldline that we are integrating.
  /**
   * Beware this is not a SmartPointer. Make sure line_ still exists
   * when calling nestStep().
   */
  Worldline * line_;
  double delta_; ///< Integration step (current in case of #adaptive_).
  bool adaptive_; ///< Whether to use an adaptive step
  double norm_; ///< Current norm of the 4-velocity.
  double normref_; ///< Initial norm of the 4-velocity.
  /// The Metric in this end of the Universe.
  /**
   * Taken from Worldline::line_, never updated.
   */
  Gyoto::SmartPointer<Gyoto::Metric::Generic> gg_;

 public:
  /**
   * \brief Normal constructor
   *
   * Sets #line_
   */
  Generic(Worldline *parent);

  /**
   * \brief Virtual destructor
   */
  virtual ~Generic();

  /**
   * \brief Deep copy
   *
   * Derived classes must implement it
   */
  virtual Generic * clone(Worldline*newparent) const =0 ;

  /**
   * \param line The Worldline that we are integrating. Sets:
   * Worldline::line_, Worldline::gg_, Worldline::adaptive_.
   * \param coord Initial coordinate.
   * \param delta Integration step. Sign determines direction.
   */
  virtual void init(Worldline * line, const double *coord, const double delta);

  /**
   * \brief Cache whatever needs to be cached
   *
   * This is called by all the methods in Worldline each time an
   * member that could be cached in Worldline::state_
   * changes. Therefore, user code should normally not have to call
   * it.
   */
  virtual void init();

  /**
   * \brief Check norm
   *
   * Issue a warning using #GYOTO_SEVERE if norm is
   * drifting. nextStep() implementations should call it.
   */
  virtual void checkNorm(double coord[8]);

  /**
   * \brief Return the integrator kind
   */
  virtual std::string kind()=0;

  /// Make one step.
  /**
   * \param[out] coord Next position-velocity;
   * \param[in] h1max maximum step in case of adaptive integration
   */
  virtual int nextStep(double *coord, double h1max=GYOTO_DEFAULT_DELTA_MAX)=0;

  /// Make one step of exactly this size.
  /**
   * doStep() is meant to refine a computation made using
   * nextStep(). In particular, there is no checking for norm
   * conservation.
   *
   * \param[in] coordin current position-velocity;
   * \param[in] step  exact step to use.
   * \param[out] coordout next position-velocity;
   */
  virtual void doStep(double const coordin[8], 
		      double step,
		      double coordout[8])=0;
};

/**
 * \class Gyoto::Worldline::IntegState::Legacy
 * \brief Home-brewed integrator
 *
 * The integrator used by this IntegState::Generic implementation is
 * actually implemented in Metric::Generic::myrk4_adaptive(). It does
 * not use most of the tuning parameters Worldline, it uses the
 * homonym parameters in Metric::Generic instead. to use this
 * integrator, pass "Legacy" to Worldline::integrator(std::string
 * type).
 */
class Gyoto::Worldline::IntegState::Legacy : public Generic {
  friend class Gyoto::SmartPointer<Gyoto::Worldline::IntegState::Legacy>;

 private:
  double coord_[8]; ///< Previously determined coordinate.

 public:
  /// Constructor

  Legacy(Worldline *parent);
  Legacy * clone(Worldline*newparent) const ;
  using Generic::init;
  void init(Worldline * line, const double *coord, const double delta);
  virtual std::string kind();

  virtual int nextStep(double *coord, double h1max=1e6);

  virtual void doStep(double const coordin[8], 
		      double step,
		      double coordout[8]);

  virtual ~Legacy();
};

#ifdef HAVE_BOOST
/**
 * \class Gyoto::Worldline::IntegState::Boost
 * \brief Boost integrator
 *
 * This Worldline::IntegState::Generic implementation provides several
 * integrators from the boost::numeric::odeint library. To select it,
 * pass one of "runge_kutta_cash_karp54", "runge_kutta_fehlberg78",
 * "runge_kutta_dopri5", or "runge_kutta_cash_karp54_classic" to
 * Worldline::integrator(std::string type).
 */
class Gyoto::Worldline::IntegState::Boost : public Generic {
  friend class Gyoto::SmartPointer<Gyoto::Worldline::IntegState::Boost>;
 public:
  /**
   * \brief Enum to represent the integrator flavour
   */
  enum Kind {runge_kutta_cash_karp54,
	     runge_kutta_fehlberg78,
	     runge_kutta_dopri5,
	     runge_kutta_cash_karp54_classic };
 private:
  /// Integrator flavour
  Kind kind_;

  /// Stepper used by the adaptive-step integrator
  std::function<boost::numeric::odeint::controlled_step_result
    (std::array<double,8>&, double&, double&)> try_step_;

  /// Stepper used by the non-adaptive-step integrator
  std::function<void(std::array<double,8>&, double)> do_step_;
 public:
  /// Constructor
  /**
   * Since this IntegState::Generic implementation can actually be
   * used to implement several distinct integrators, it is necessary
   * to specify which one is meant.
   */
  Boost(Worldline* parent, std::string type);

  /// Constructor
  /**
   * Since this IntegState::Generic implementation can actually be
   * used to implement several distinct integrators, it is necessary
   * to specify which one is meant.
   */
  Boost(Worldline* parent, Kind type);
  Boost * clone(Worldline* newparent) const ;
  virtual ~Boost();
  virtual void init();
  virtual void init(Worldline * line, const double *coord, const double delta);
  virtual int nextStep(double *coord, double h1max=1e6);
  virtual void doStep(double const coordin[8], 
		      double step,
		      double coordout[8]);
  virtual std::string kind();
  
};
#endif /// HAVE_BOOST

#endif
