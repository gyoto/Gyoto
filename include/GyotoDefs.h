/**
 * \file GyotoDefs.h
 * \brief Gyoto ubiquitous macros and typedefs
 */
/*
    Copyright 2011 Thibaut Paumard

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

#ifndef __GyotoDefs_H_ 
#define __GyotoDefs_H_ 

#include "GyotoConfig.h"
#include <float.h>

/**
 * \brief Replacement for GNU extension sincos
 *
 * If #HAVE_SINCOS is undefined, Gyoto provides a trivial
 * implementation.
 * \param[in] t Angle in radian;
 * \param[out] s Adress where sin(t) should be stored; 
 * \param[out] c Adress where cos(t) should be stored.
 */
#if !HAVE_SINCOS
#define sincos(t, s, c) *s=sin(t); *c=cos(t)
#else
# ifdef DOXYGEN_RUN
#  define sincos(t, s, c) (undefined)
# endif
#endif

/* Typedef for various Gyoto data types */
namespace Gyoto {
  //\{
  /**
   * \name Observable quantities
   *
   * Individual quantities are represented as a variable of this
   * type. A combination of quantities is an ored list of Quantity_t,
   * e.g.
   *
   * \code
   * GYOTO_QUANTITY_INTENSITY | GYOTO_QUANTITY_EMISSIONTIME
   * \endcode
   * 
   * To check wheter a given quantity is listed in a Quantity_t
   * variable quant:
   * \code
   * if (quant & GYOTOQUANTITY_EMISSION) ...
   * \endcode
   *
   * List of all possible Quantity_t individual values and the
   * corresponding string (see Gyoto::Scenery) with corresponding XML
   * name:
   */

  /// Type for observabke quantities
  typedef unsigned int Quantity_t;

  /* Generic */

  /// Intensity: I<SUB>&nu;</SUB> at Scenery::freq_obs_.
#define GYOTO_QUANTITY_INTENSITY      1
  /// EmissionTime: Emission date.
#define GYOTO_QUANTITY_EMISSIONTIME   2
  /// MinDistance: Behaves like minimal distance between Photon and Astrobj.
  /**
   * Not always exactly a distance, though. 
   */
#define GYOTO_QUANTITY_MIN_DISTANCE   4
  /// FirstDmin: First Photon-Astrobj distance local minimum while integrating back in time.
#define GYOTO_QUANTITY_FIRST_DMIN     8
  /// Redshift: &nu;<SUB>obs</SUB>/&nu;<SUB>em</SUB>.
#define GYOTO_QUANTITY_REDSHIFT      16
  /// ImpactCoords: Astrobj and Photon 8-coordinates at emission.
  /**
   * A 16-element vector. See Gyoto::Quantity_t.
   */ 
#define GYOTO_QUANTITY_IMPACTCOORDS  32
  /// Spectrum: I<SUB>&nu;</SUB> at each frequency in Scenery::screen_->getMidpoints().
#define GYOTO_QUANTITY_SPECTRUM     512
  /// Spectrum: &int;<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP>I<SUB>&nu;</SUB> d&nu; in each frequency channel in Scenery::screen_.
#define GYOTO_QUANTITY_BINSPECTRUM 1024
  /* Astrobj-specific */
  /// User1: Gyoto::Astrobj specific Gyoto::Quantity_t
#define GYOTO_QUANTITY_USER1        32768
  /// User2: Gyoto::Astrobj specific Gyoto::Quantity_t
#define GYOTO_QUANTITY_USER2        16384
  /// User3: Gyoto::Astrobj specific Gyoto::Quantity_t
#define GYOTO_QUANTITY_USER3         8192
  /// User4: Gyoto::Astrobj specific Gyoto::Quantity_t
#define GYOTO_QUANTITY_USER4         4096
  /// User5: Gyoto::Astrobj specific Gyoto::Quantity_t
#define GYOTO_QUANTITY_USER5         2048
  //\}

  /**
   * \name Gyoto messages
   * \brief Controling which messages are shown to the user
   *
   * The user should be able to choose which messages are shown to
   * her. In Gyoto, this is determined by a user-settable verbosity
   * level (see Gyoto::verbose()) and a user-settable debug mode (see
   * Goyto::debug()).
   *
   * The following macros define various debug and verbosity level and
   * provide short-cuts to display formatted messages only at a given
   * verbosity level or in debug mode.
   */
  //\{
  /// Type for verbosity levels
  typedef unsigned int Verbosity_t;

  /// Default debug mode
#define GYOTO_DEFAULT_DEBUG_MODE 0

  /// Quiet Gyoto::Verbosity_t
  /**
   * Only very few messages may be output to stdout at this level.
   *
   * Use #GYOTO_QUIET to display messages at this level.
   */
#define GYOTO_QUIET_VERBOSITY   1

  /// Severe warnings
  /**
   * Severe warning messages are output to stderr if
   * Gyoto::verbose()&ge;GYOTO_SEVERE_VERBOSITY.
   *
   * Use #GYOTO_SEVERE to display messages at this level.
   */
#define GYOTO_SEVERE_VERBOSITY  3

  /// Warnings
  /**
   * Warning messages are output to stderr if
   * Gyoto::verbose()&ge;GYOTO_WARNING_VERBOSITY.
   *
   * Use #GYOTO_WARNING to display messages at this level.
   */
#define GYOTO_WARNING_VERBOSITY  GYOTO_SEVERE_VERBOSITY

  /// Default verbosity level
  /**
   * Normal messages are output to stdout if
   * Gyoto::verbose()&ge;GYOTO_DEFAULT_VERBOSITY.
   *
   * Use #GYOTO_MSG to display messages at this level.
   */
#define GYOTO_DEFAULT_VERBOSITY 5

  /// Informative messages
  /**
   * Informative messages are output to stderr if
   * Gyoto::verbose()&ge;GYOTO_INFO_VERBOSITY.
   *
   * Use #GYOTO_INFO to display messages at this level.
   */
#define GYOTO_INFO_VERBOSITY   10

  /// Maximum verbosity level
  /**
   * In debug mode, all messages are displayed in addition to specific
   * debug information.
   *
   * To display debug messages, check #GYOTO_DEBUG_MODE or use
   * #GYOTO_DEBUG for instance.
   */
#define GYOTO_DEBUG_VERBOSITY 3000

  /// Display a message to stdout even in quiet mode.
  /**
   * Only very few messages may be output to stdout at this
   * level. This should be reserved to messages shown at most once in
   * a run.
   *
   * \code
   * GYOTO_QUIET << "Important message displayed once" << std::endl;
   * \endcode
   */
#define GYOTO_QUIET  if (Gyoto::verbose() >= GYOTO_QUIET_VERBOSITY) std::cout 

  /// Display a severe level message to stderr.
  /**
   * \code
   * GYOTO_SEVERE << "Important warning" << std::endl;
   * \endcode
   */
#define GYOTO_SEVERE if(Gyoto::verbose()>=GYOTO_SEVERE_VERBOSITY) std::cerr<<"SEVERE: "

  /// Display a warning level message to stderr.
  /**
   * \code
   * GYOTO_WARNING << "Warning" << std::endl;
   * \endcode
   */
#define GYOTO_WARNING if(Gyoto::verbose()>=GYOTO_SEVERE_VERBOSITY) std::cerr<<"WARNING: "

  /// Display normal message to stdout.
  /**
   * The message will by default be shown to the user. To be reserved
   * to messages shown at most once per Gyoto::Screen row.
   *
   * \code
   * GYOTO_MSG << "Message" << std::endl;
   * \endcode
   */
#define GYOTO_MSG    if (Gyoto::verbose() >= GYOTO_DEFAULT_VERBOSITY) std::cout 

  /// Display informative message to stderr.
  /**
   * Message will be shown to the user only if Gyoto::verbose() has
   * been explicitely raised.
   *
   * \code
   * GYOTO_MSG << "Message" << std::endl;
   * \endcode
   */
#define GYOTO_INFO   if (Gyoto::verbose() >= GYOTO_INFO_VERBOSITY) std::cerr<<"INFO: "

  /// Unit ignored because libudunits2 was disabled.
  /**
   * Use this warning when a conversion has been dropped due to
   * #HAVE_UDUNITS being undefined.
   *
   * \param from From unit.
   * \param to To unit.
   */
#define GYOTO_WARNING_UDUNITS(from, to) \
  GYOTO_WARNING << "unit ignored (trying to convert from \"" << from \
		<< "\" to "					     \
		<< to \
		<< "\"), you may have more chance recompiling Gyoto with --with-udunits\n"

  /// Output expression value in debug mode
  /**
   * Output, only in debug mode both code and value that code
   * yield. For instance:
   * \code
   * int a=1, b=2;
   * GYOTO_DEBUG_EXPR(a+b);
   * \endcode
   * will essentially output, only in debug mode:
   * \code
   * DEBUG: <function signature>: a+b=3
   * \endcode
   */
#define GYOTO_DEBUG_EXPR(a) GYOTO_DEBUG << #a << "=" << a << std::endl

  /// Output array content in debug mode
  /**
   * Output, only in debug, name and content of array.
   * For instance:
   * \code
   * int a[]= {1, 2, 3};
   * GYOTO_DEBUG_ARRAY(a,3);
   * \endcode
   * will essentially output, only in debug mode:
   * \code
   * DEBUG: <function signature>: a=[1,2,3]
   * \endcode
   *
   * \param a Array
   * \param n Number of elements to show (array must be at least this
   * size).
   */
#define GYOTO_DEBUG_ARRAY(a,n) if (GYOTO_DEBUG_MODE) {            \
    std::cerr << "DEBUG: " << __PRETTY_FUNCTION__ << ": "         \
	      << #a << "=[" << a[0] ;				  \
    for (size_t _gyoto_debug_array_i=1; _gyoto_debug_array_i < n; ++_gyoto_debug_array_i) \
      std::cerr << "," << a[_gyoto_debug_array_i] ;			\
    std::cerr << "]" << std::endl ;}

  /// Display debug message
  /**
   * Message is deiplayed only if #GYOTO_DEBUG_MODE is true and is
   * prepended with the word "DEBUG:" and the signature of the
   * function.
   *
   * \code
   * GYOTO_DEBUG << "message" << endl;
   * \endcode
   */
#define GYOTO_DEBUG  if (GYOTO_DEBUG_MODE) std::cerr << "DEBUG: " << __PRETTY_FUNCTION__ << ": "

  /// Start debug-only block.
  /**
   * Code between GYOTO_IF_DEBUG and #GYOTO_ENDIF_DEBUG is only
   * executed in debug mode.
   */
#define GYOTO_IF_DEBUG if (GYOTO_DEBUG_MODE) {

  /// End debug-only block.
  /**
   * Code between #GYOTO_IF_DEBUG and GYOTO_ENDIF_DEBUG is only
   * executed in debug mode.
   */
#define GYOTO_ENDIF_DEBUG }

  /// Whether debug mode is activated (run-time).
#define GYOTO_DEBUG_MODE Gyoto::debug()

  //\}
  //\{
  /**
   * \name Coordinate system kind
   * GYOTO_COORDKIND_CARTESIAN or GYOTO_COORDKIND_SPHERICAL
   *
   * Every Gyoto::Metric has a coordinate system kind. It can be used
   * by functions which need to express coordinates always in
   * spherical or always in Cartesian form, with trivial conversion
   * between the two.
   */
  /**
   * \brief Type for coordinate system kinds
   */ 
  typedef unsigned int CoordKind_t;
#define GYOTO_COORDKIND_UNSPECIFIED 0 ///< Unspecified coordinate kind
#define GYOTO_COORDKIND_CARTESIAN 1 ///< Cartesian-like coordinate system
#define GYOTO_COORDKIND_SPHERICAL 2 ///< Spherical-like coordinate system
  //\}
}

//{
/**
 * \name Default values for various things
 */

#define GYOTO_DEFAULT_X_SIZE 1024 ///< Default size for arrays in a Worldline

/**
 * \brief Default value for the initial step in the integration loop.
 *
 * Since the step is (most of the time) adaptive, this default only
 * has little influence, but sometimes, it matters. Also used in
 * Scenery.
 */
#define GYOTO_DEFAULT_DELTA 0.01

/**
 * \brief Default value for the maximum step in the integration loop.
 */
#define GYOTO_DEFAULT_DELTA_MAX DBL_MAX

/**
 * \brief Default value for the minimum step in the integration loop.
 */
#define GYOTO_DEFAULT_DELTA_MIN 1e-2


/**
 * \brief Default value for Gyoto::Worldline::maxiter_
 */
#define GYOTO_DEFAULT_MAXITER 100000

/**
 * \brief Precision on the determination of a date
 *
 * E.g. in Gyoto::Photon::findMin(), Gyoto::Photon::findValue().
 */
#define GYOTO_T_TOL 1e-4


/// \brief Default value for Screen::dmax_
#define GYOTO_SCREEN_DMAX DBL_MAX

//For displays with setw and setprecision
/// Precision when outputting double values
#define GYOTO_PREC  15
/// \brief Field width when outputting double values
#define GYOTO_WIDTH 25

/* Plugins Stuff */
/// \brief Default list of default plug-ins to load...
#ifndef GYOTO_DEFAULT_PLUGINS
#define GYOTO_DEFAULT_PLUGINS "stdplug,nofail:lorene"
#endif

#ifndef GYOTO_PLUGIN_SFX
#define GYOTO_PLUGIN_SFX "so"
#endif

//\}

//\{
/**
 *\name Physical constants
 */
/// \brief Celerity of light (m/s)
#define GYOTO_C          299792458.
/// \brief Celerity of light (cm/s)
#define GYOTO_C_CGS 2.99792458e10
/// \brief Gravitational constant (SI = m^3 * kg^-1 * s-2)
#define GYOTO_G 6.67428e-11
/// \brief Gravitational constant (cgs: cm^3 * g^-1 * s-2)
#define GYOTO_G_CGS 6.67428e-8
/// \brief G/c^2=6.67428e-11/299792458.^2
#define GYOTO_G_OVER_C_SQUARE 7.426138e-28
/// \brief Planck's constant (h) in SI (J.s=kg.m^2/s) 
#define GYOTO_PLANCK 6.62606896e-34
/// \brief Planck's constant (h) in c.g.s (g.cm^2/s) 
#define GYOTO_PLANCK_CGS 6.62606896e-27
/// \brief h/c^2 in SI (kg.s)
#define GYOTO_PLANCK_OVER_C_SQUARE 7.372496e-51
/// \brief Boltzmann's constant (k) in SI (J/K)
#define GYOTO_BOLTZMANN 1.3806504e-23
/// \brief Boltzmann's constant (k) in cgs (erg/K)
#define GYOTO_BOLTZMANN_CGS 1.3806504e-16
/// \brief Stefan-Boltzmann's constant (sigma) in cgs (erg/cm2/s/K4)
#define GYOTO_STEFANBOLTZMANN_CGS 5.670373e-5
/// \brief h/k (K.s = K/Hz)
#define GYOTO_PLANCK_OVER_BOLTZMANN 4.7992373e-11
/// \brief ideal gas constant R in SI
#define GYOTO_GAS_CST 8.3144621
/// \brief ideal gas constant R in erg/(K mol)
#define GYOTO_GAS_CST_CGS 8.3144621e7
/// \brief Avogadro constant
#define GYOTO_AVOGADRO 6.0221413e23
/// \brief Thomson cross-section in cgs
#define GYOTO_THOMSON_CGS 6.6524e-25
/// \brief Fine structure constant (=1/137)
#define GYOTO_ALPHA_F 0.00729927
/// \brief proton mass in cgs
#define GYOTO_PROTON_MASS_CGS 1.67262158e-24
/// \brief electron mass in cgs
#define GYOTO_ELECTRON_MASS_CGS 9.10938188e-28
/// \brief electron classical radius in cgs
#define GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS 2.8179e-13
/// \brief elementary charge in cgs (erg^{1/2} cm^{1/2})
#define GYOTO_ELEMENTARY_CHARGE_CGS 4.80320427e-10
/// \brief Euler-Mascheroni constant
#define GYOTO_EULER_MASCHERONI  0.577216
/// \brief atomic mass unit in cgs
#define GYOTO_ATOMIC_MASS_UNIT_CGS 1.660537781e-24
/// \brief Factor to convert I<SUB>&nu;</SUB> from c.g.s. to SI
#define GYOTO_INU_CGS_TO_SI 0.001

/// \brief Sun mass (kg)
#define GYOTO_SUN_MASS    1.98843e30
/// \brief Sun mass (g)
#define GYOTO_SUN_MASS_CGS    1.98843e33
/// \brief Sun radius (m)
#define GYOTO_SUN_RADIUS     6.955e8
/// \brief Kiloparsec (m)
#define GYOTO_KPC        3.08568025e19
/// \brief Astronomical Unit (m)
#define GYOTO_ASTRONOMICAL_UNIT 1.49597870700e11
/// \brief Light-year (m)
#define GYOTO_LIGHT_YEAR 9.4607304725808e15

/// \brief Convert from radians to degrees
#define GYOTO_RADEG 57.2957795130823
/// \brief Convert from degrees to radians
#define GYOTO_DEGRAD 0.0174532925199433
/// \brief Convert from arcminutes to radians
#define GYOTO_MINRAD 2.908882086657216e-04
/// \brief Convert from arcseconds to radians
#define GYOTO_SECRAD 4.848136811095360e-06
/// \brief Convert from milliarcseconds to radians
#define GYOTO_MASRAD 4.848136811095360e-09
/// \brief Convert from microarcseconds to radians
#define GYOTO_MUASRAD 4.848136811095360e-12

/// \brief Convert from eV to Hz
#define GYOTO_eV2Hz 2.417989348e+14

//\}

#ifndef GYOTO_NO_DEPRECATED
#warning Using deprecated method names.\
  Define GYOTO_NO_DEPRECATED to disable.
//\{
/**
 *\name Renamed methods
 *
 * Define GYOTO_NO_DEPRECATED to disable these macros and the warning. For instance:
 * \code
 * make CPPFLAGS=-DGYOTO_NO_DEPRECATED
 * \endcode
 */
# define getMetric       metric
# define setMetric       metric 
# define setScreen       screen
# define getScreen       screen
# define getRmax         rMax
# define setRmax         rMax
# define getMass         mass
# define setMass         mass
# define getCoordKind    coordKind
# define setCoordKind    coordKind
# define getKind         kind
# define setKind         kind
# define getSpin         spin
# define setSpin         spin
# define getIntegKind    integKind
# define setIntegKind    integKind
# define getFileName     fileName
# define setFileName     fileName
# define getDistance     distance
# define setDistance     distance
# define getPALN         PALN
# define setPALN         PALN
# define getArgument     argument
# define setArgument     argument
# define getInclination  inclination
# define setInclination  inclination
# define getAstrobj      astrobj
# define setAstrobj      astrobj
# define getSpectrometer spectrometer
# define setSpectrometer spectrometer
# define getSpectrum     spectrum
# define setSpectrum     spectrum
# define getOpacity      opacity
# define setOpacity      opacity
# define setDelta        delta
# define getDelta        delta
# define setDelta        delta
# define getDelta        delta
# define setDelta0       delta0
# define getDelta0       delta0
# define setAlpha0       alpha0
# define getAlpha0       alpha0
//\}
#endif

#endif
