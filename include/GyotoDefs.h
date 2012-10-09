/**
 * \file GyotoDefs.h
 * \brief #defines and typedefs for Gyoto
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


/* Typedef for various Gyoto data types */
namespace Gyoto {
  /**
   * \brief Spectrometer kind
   * 
   * One of: GYOTO_SPECTRO_KIND_NONE, GYOTO_SPECTRO_KIND_FREQ,
   * GYOTO_SPECTRO_KIND_FREQLOG, GYOTO_SPECTRO_KIND_WAVE,
   * GYOTO_SPECTRO_KIND_WAVELOG
   *  
   */
  typedef unsigned int SpectroKind_t;
#define GYOTO_SPECTRO_KIND_NONE    0
#define GYOTO_SPECTRO_KIND_FREQ    1 ///< Spectrometer kind="freq"
#define GYOTO_SPECTRO_KIND_FREQLOG 2 ///< Spectrometer kind="freqlog"
#define GYOTO_SPECTRO_KIND_WAVE    3 ///< Spectrometer kind="wave"
#define GYOTO_SPECTRO_KIND_WAVELOG 4 ///< Spectrometer kind="wavelog"

  /**
   * \brief Observable quantities
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
   * corresponding string (see Scenery):
   *
   * - GYOTO_QUANTITY_INTENSITY   : Intensity
   * - GYOTO_QUANTITY_EMISSIONTIME: EmissionTime
   * - GYOTO_QUANTITY_MIN_DISTANCE: MinDistance
   * - GYOTO_QUANTITY_FIRST_DMIN  : FirstDmin
   * - GYOTO_QUANTITY_REDSHIFT    : RedShift
   * - GYOTO_QUANTITY_SPECTRUM    : Spectrum
   * - GYOTO_QUANTITY_BINSPECTRUM : BinSpectrum
   * - GYOTO_QUANTITY_IMPACTCOORDS: ImpactCoords
   * - GYOTO_QUANTITY_USER1       : User1
   * - GYOTO_QUANTITY_USER2       : User2
   * - GYOTO_QUANTITY_USER3       : User3
   * - GYOTO_QUANTITY_USER4       : User4
   * - GYOTO_QUANTITY_USER5       : User5
   */
  typedef unsigned int Quantity_t;
/* Generic */
#define GYOTO_QUANTITY_INTENSITY      1
#define GYOTO_QUANTITY_EMISSIONTIME   2
#define GYOTO_QUANTITY_MIN_DISTANCE   4
#define GYOTO_QUANTITY_FIRST_DMIN     8
#define GYOTO_QUANTITY_REDSHIFT      16
#define GYOTO_QUANTITY_IMPACTCOORDS  32
#define GYOTO_QUANTITY_SPECTRUM     512
#define GYOTO_QUANTITY_BINSPECTRUM 1024
/* Astrobj-specific */
#define GYOTO_QUANTITY_USER1        32768
#define GYOTO_QUANTITY_USER2        16384
#define GYOTO_QUANTITY_USER3         8192
#define GYOTO_QUANTITY_USER4         4096
#define GYOTO_QUANTITY_USER5         2048

  /**
   * \brief Verbosity levels
   */
  typedef unsigned int Verbosity_t;
#define GYOTO_DEFAULT_DEBUG_MODE 0
#define GYOTO_QUIET_VERBOSITY   1
#define GYOTO_SEVERE_VERBOSITY  3
#define GYOTO_WARNING_VERBOSITY  GYOTO_SEVERE_VERBOSITY
#define GYOTO_DEFAULT_VERBOSITY 5
#define GYOTO_INFO_VERBOSITY   10
#define GYOTO_DEBUG_VERBOSITY 3000

#define GYOTO_QUIET  if (Gyoto::verbose() >= GYOTO_QUIET_VERBOSITY) \
    std::cout 
#define GYOTO_SEVERE if(Gyoto::verbose()>=GYOTO_SEVERE_VERBOSITY) \
    std::cerr<<"SEVERE: "
#define GYOTO_WARNING if(Gyoto::verbose()>=GYOTO_SEVERE_VERBOSITY) \
    std::cerr<<"WARNING: "
#define GYOTO_MSG    if (Gyoto::verbose() >= GYOTO_DEFAULT_VERBOSITY) \
    std::cout 
#define GYOTO_INFO   if (Gyoto::verbose() >= GYOTO_INFO_VERBOSITY) \
    std::cerr<<"INFO: "

#define GYOTO_DEBUG_EXPR(a) GYOTO_DEBUG << #a << "=" << a << std::endl
#define GYOTO_DEBUG_ARRAY(a,n) if (GYOTO_DEBUG_MODE) {            \
    std::cerr << "DEBUG: " << __PRETTY_FUNCTION__ << ": "         \
	      << #a << "=[" << a[0] ;				  \
    for (size_t i=1; i < n; ++i) std::cerr << "," << a[i] ;	  \
    std::cerr << "]" << std::endl ;}
#define GYOTO_DEBUG  if (GYOTO_DEBUG_MODE) std::cerr << "DEBUG: "	\
					      << __PRETTY_FUNCTION__ << ": "
#define GYOTO_IF_DEBUG if (GYOTO_DEBUG_MODE) {
#define GYOTO_ENDIF_DEBUG }

# define GYOTO_DEBUG_MODE Gyoto::debug()

  /**
   * \brief Coordinate system kinds
   * GYOTO_COORDKIND_CARTESIAN or GYOTO_COORDKIND_SPHERICAL
   */ 
  typedef unsigned int CoordKind_t;
#define GYOTO_COORDKIND_UNSPECIFIED 0 ///< Unspecified coordinate kind
#define GYOTO_COORDKIND_CARTESIAN 1 ///< Cartesian-like coordinate system
#define GYOTO_COORDKIND_SPHERICAL 2 ///< Spherical-like coordinate system

}

/* Default values for various things */

#define GYOTO_DEFAULT_X_SIZE 1024 ///< Default size for arrays in a Worldline


/**
 * Default value for the initial step in the integration loop. Since
 * the step is (most of the time) adaptive, this default only has
 * little influence, but sometimes, it matters. Also used in Scenery.
 */
#define GYOTO_DEFAULT_DELTA 0.01

/**
 * Precision on the determination of a date (e.g. in
 * Photon::findMin(), Photon::findValue()).
 */
#define GYOTO_T_TOL 1e-4

/* Plugins Stuff */
#ifndef GYOTO_DEFAULT_PLUGINS
#define GYOTO_DEFAULT_PLUGINS "stdplug,nofail:lorene"
#endif

#ifndef GYOTO_PLUGIN_SFX
#ifdef __APPLE__
#define GYOTO_PLUGIN_SFX "dylib"
#else
#define GYOTO_PLUGIN_SFX "so"
#endif
#endif

/* Physical constants */

/// \brief Celerity of light (m/s)
#define GYOTO_C          299792458.
/// \brief Celerity of light (cm/s)
#define GYOTO_C_CGS 2.99792458e10
/// \Brief Gravitational constant (SI = m^3 * kg^-1 * s-2)
#define GYOTO_G 6.67428e-11
/// \Brief Gravitational constant (cgs: cm^3 * g^-1 * s-2)
#define GYOTO_G_CGS 6.67428e-8
/// \brief G/c^2=6.67428e-11/299792458.^2
#define GYOTO_G_OVER_C_SQUARE 7.426138e-28
/// \brief Planck's constant (h) in SI (J.s=kg.m^2/s) 
#define GYOTO_PLANCK 6.62606896e-34
#define GYOTO_PLANCK_CGS 6.62606896e-27
/// \brief h/c^2 in SI (kg.s)
#define GYOTO_PLANCK_OVER_C_SQUARE 7.372496e-51
/// \brief Boltzmann's constant (k) in SI (J/K)
#define GYOTO_BOLTZMANN 1.3806504e-23
/// \brief Boltzmann's constant (k) in cgs (erg/K)
#define GYOTO_BOLTZMANN_CGS 1.3806504e-16
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
#define GYOTO_DEGRAD 0.0174532925199433
#define GYOTO_MINRAD 2.908882086657216e-04
#define GYOTO_SECRAD 4.848136811095360e-06
#define GYOTO_MASRAD 4.848136811095360e-09
#define GYOTO_MUASRAD 4.848136811095360e-12

/// \brief Default value for Screen::dmax_
#define GYOTO_SCREEN_DMAX 1e7

//For displays with setw and setprecision
#define GYOTO_PREC  15
#define GYOTO_WIDTH 25

#endif
