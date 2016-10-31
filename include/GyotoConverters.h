/**
 * \file GyotoConverters.h
 * \brief GYOTO converters
 *
 *  As a (recommended) option, Gyoto can use the libudunits2 library
 *  by UCAR/Unidata
 *  (http://www.unidata.ucar.edu/software/udunits/udunits-2/udunits2.html)
 *  to perform conversions to and from various units. Since udunits
 *  parses units, the following are equivalent for instance:
 *  "mJy.sr-1", "mJy/sr", "1e-3Jy/sr". Gyoto considers every string as
 *  UTF-8: either use UTF-8 or stay with ASCII. This means "µ" can be
 *  used for "micro" and powers can be noted as superscripts: m³,
 *  pix².
 *
 *  In addition to the stock units known by udunits2, Gyoto registers
 *  the following (which may be used in any context): Jansky (symbol
 *  Jy), symbol "pc" for parsec, "sunradius", "sunmass", symbol "as"
 *  for "arcsec".
 *  
 *  Other units are context-sensitive: "geometrical" allows converting
 *  between geometrical units and other legnth units, but only when a
 *  Metric is defined, and may not be used (yet) in compound
 *  units. Likewise, "geometrical_time" can be used as duration unit
 *  whenever a Metric is defined, but not in a compound unit.
 *
 *  When a Screen is defined, "pix" can be used as an angle unit (you
 *  need to call Screen::mapPixUnit() and Screen::unmapScreenUnit(),
 *  which is done automatically in certain contexts).
 *
 *  Units can often be specified using the "unit" XML attribute:
 * \code
 * <Meric kind="KerrBL">
 *    <Mass unit="sunmass">
 *       4e6
 *    </Mass>
 * </Metric>
 * <Screen>
 *    <Distance unit="kpc">
 *       8
 *    </Distance>
 *    <FieldOfView unit="µas">
 *       150
 *    </FieldOfView>
 * </Metric>
 * \endcode
 *
 * Units for output quantities are specified after the name of the
 * quantity, in squared brackets, e.g.:
 * \code
 *   <Quantities>
 *     Spectrum[erg.s-1.cm-2.sr-1.Hz-1]
 *   </Quantities>
 * \endcode
 * or
 * \code
 *   <Quantities>
 *     Spectrum[mJy/pix²]
 *   </Quantities>
 * \endcode
 */

/*
    Copyright 2011-2016 Thibaut Paumard

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

#ifndef __GyotoConverters_H_ 
#define __GyotoConverters_H_ 

#include <GyotoUtils.h>
#include <GyotoSmartPointer.h>

#ifdef HAVE_UDUNITS
#include <udunits2.h>
#endif

#include <string>
#include <sstream>

namespace Gyoto {
  namespace Metric {
    class Generic;
  }
  namespace Units {
    /**
     * \namespace Gyoto::Units
     * \brief Units-related classes and functions
     */ 
#ifdef HAVE_UDUNITS
    /**
     * \class Gyoto::Units::Unit
     * \brief Wrapper around ut_unit from udunits
     *
     * Gyoto::Units::Unit objects usually cast seamlessly to and from
     * udunits2 ut_unit* and std::string.
     */
    class Unit;

    /**
     * \class Gyoto::Units::Converter
     * \brief Wrapper around ut_converter from udunits
     *
     * A Gyoto::Units::Converter object is a functor and can be used
     * to convert efficiently between the two units specified at
     * instantiation time:
     *
     * \code
     * Gyoto::Units::Unit from_unit ("erg.s-1.cm-2.sr-1.Hz-1"),
     *                    to_unit   ("mJy/µas²");
     * double data_in[1000], data_out[1000];
     * Converter conv(from_unit, to_unit);
     * for (size_t i=0; i<1000; ++i) data_out[i] = conv(data_in[i]);
     * \endcode
     *
     * Since std::string cast automatically to Gyoto::Units::Unit
     * object, this is equivalent:
     *
     * \code
     * Converter conv("erg.s-1.cm-2.sr-1.Hz-1", "mJy/µas²");
     * for (size_t i=0; i<1000; ++i) data_out[i] = conv(data_in[i]);
     * \endcode
     */
    class Converter;

    /**
     * \brief Retrieve the unit system used in all of Gyoto
     */
    ut_system * getSystem();
#endif

    /**
     * \brief Load and initialize all (non-context-sensitive) units
     *
     * If udunits is used (preprocessor macro HAVE_UDUNITS), Init()
     * initializes the ut_system used throughout Gyoto and maps a few
     * additional units to the unit system.
     */
    void Init();

    /**
     * \brief Convert from arbitrary length unit to meters
     *
     * Convert value from unit represented by "unit" to meters.
     *
     * If gg is provided (and not NULL), use it to interpret the
     * string "geometrical" as representing
     * gg->Gyoto::Metric::Generic::unitLength().
     *
     * ToMeters() will also convert time, frequency and energy units
     * to meters (as in frequency -> wavelength).
     *
     * \param value
     *   (double) the value to convert, expressed according to "unit"
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "km", "sunradius" or "geometrical". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &,
     *   NULL if not specified) optional metric to interpret
     *   "geometrical".
     *
     * \return value, expressed in meters.
     */
    double ToMeters(double value, const std::string &unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);

    /**
     * \brief Convert to arbitrary length unit from meters
     *
     * Convert value to unit represented by "unit" from meters.
     *
     * If gg is provided (and not NULL), use it to interpret the
     * string "geometrical" as representing
     * gg->Gyoto::Metric::Generic::unitLength().
     *
     * ToMeters() will also convert to time, frequency and energy
     * units (as in wavelength -> frequency).
     *
     * \param value
     *   (double) the value to convert, expressed in meters.
     * \param unit (std::string) the "unit" to which to convert,
     *   e.g. "km", "sunradius" or "geometrical". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &,
     *   NULL if not specified) optional metric to interpret
     *   "geometrical".
     *
     * \return value, expressed in "unit".
     */
    double FromMeters(double value, const std::string &unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);

    /**
     * \brief Convert from arbitrary time unit to seconds
     *
     * Convert value from unit represented by "unit" to seconds.
     *
     * If gg is provided (and not NULL), use it to interpret the
     * string "geometrical_time" as representing
     * gg->Gyoto::Metric::Generic::unitLength()/GYOTO_C.
     *
     * \param value
     *   (double) the value to convert, expressed according to "unit"
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "s", "yr" or "geometrical_time". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &,
     *   NULL if not specified) optional metric to interpret
     *   "geometrical".
     *
     * \return value, expressed in seconds.
     */
    double ToSeconds(double value, const std::string &unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);

    /**
     * \brief Convert to arbitrary time unit from seconds
     *
     * Convert value to unit represented by "unit" from seconds.
     *
     * If gg is provided (and not NULL), use it to interpret the
     * string "geometrical_time" as representing
     * gg->Gyoto::Metric::Generic::unitLength()/GYOTO_C.
     *
     * \param value
     *   (double) the value to convert, expressed in seconds.
     * \param unit (std::string) the "unit" to which to convert,
     *   e.g. "s", "yr" or "geometrical_time". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &,
     *   NULL if not specified) optional metric to interpret
     *   "geometrical".
     *
     * \return value, expressed in "unit".
     */
    double FromSeconds(double value, const std::string &unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);

    /**
     * \brief Convert from arbitrary mass unit to kilograms
     *
     * Convert value from unit represented by "unit" to kilograms.
     *
     * \param value
     *   (double) the value to convert, expressed according to "unit"
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "g", "kg" or "sunmass". More units are
     *   supported if Gyoto was compiled with udunits support.
     *
     * \return value, expressed in kilograms.
     */
    double ToKilograms(double value, const std::string & unit);

    /**
     * \brief Convert to arbitrary mass unit from kilograms
     *
     * Convert value from unit represented by "unit" from kilograms.
     *
     * \param value
     *   (double) the value to convert, expressed inkilograms.
     * \param unit (std::string) the "unit" to which to convert,
     *   e.g. "g", "kg" or "sunmass". More units are
     *   supported if Gyoto was compiled with udunits support.
     *
     * \return value, expressed in "unit".
     */
    double FromKilograms(double value, const std::string & unit);

    /**
     * \brief Convert from arbitrary length unit to geometrical units
     *
     * Convert value from unit represented by "unit" to geometrical units.
     *
     * \param value
     *   (double) the value to convert, expressed according to "unit".
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "km", "sunradius" or "geometrical". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &)
     *   metric to interpret "geometrical".
     *
     * \return value, expressed in geometrical units.
     */
    double ToGeometrical(double value, const std::string & unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> & gg);

    /**
     * \brief Convert to arbitrary length unit from geometrical units
     *
     * Convert value to unit represented by "unit" from geometrical units.
     *
     * \param value
     *   (double) the value to convert, expressed in geometrical units.
     * \param unit (std::string) the "unit" to which to convert,
     *   e.g. "km", "sunradius" or "geometrical". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &)
     *   metric to interpret "geometrical".
     *
     * \return value, expressed in "unit".
     */
    double FromGeometrical(double value, const std::string & unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> & gg);

    /**
     * \brief Convert from arbitrary time unit to geometrical units
     *
     * \param value
     *   (double) the value to convert, expressed according to "unit".
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "s", "kyr" or "geometrical_time". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &)
     *   metric to interpret "geometrical_time".
     *
     * \return value, expressed in geometrical (time) units.
     */
    double ToGeometricalTime(double value, const std::string & unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> & gg);

    /**
     * \brief Convert to arbitrary time unit from geometrical units
     *
     * \param value
     *   (double) the value to convert, expressed in geometrical units.
     * \param unit (std::string) the "unit" to which to convert,
     *   e.g. "yr", "s" or "geometrical_time". More units are
     *   supported if Gyoto was compiled with udunits support.
     * \param gg (const Gyoto::SmartPointer<Gyoto::Metric::Generic> &)
     *   metric to interpret "geometrical_time".
     *
     * \return value, expressed in "unit".
     */
    double FromGeometricalTime(double value, const std::string & unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg);


    /**
     * \brief Convert from arbitrary frequency unit to Herz
     *
     * ToHerz will also convert from length and energy units (such as
     * "eV").
     *
     * \param value
     *   (double) the value to convert, expressed according to "unit".
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "MHz", "keV"
     *
     * \return value, expressed in "Hz".
     */
    double ToHerz(double value, const std::string &unit);

    /**
     * \brief Convert to arbitrary frequency unit from Herz
     *
     * FromHerz will also convert to length and energy units (such as
     * "eV").
     *
     * \param value
     *   (double) the value to convert, expressed according in "Hz".
     * \param unit (std::string) the "unit" from which to convert,
     *   e.g. "MHz", "keV"
     *
     * \return value, expressed in "units".
     */
    double FromHerz(double value, const std::string &unit);

#   ifdef HAVE_UDUNITS
    /**
     * \brief Is it possible to convert between unit1 and unit2?
     *
     * e.g. areConvertible("m", "kg") == 0; areConvertible("m", "km")==1.
     *
     * Warning: angle units are dimensionless, therefore e.g.
     * areConvertible("Jy", "Jy/microacsec2")==1. Numerically, "Jy" is
     * the same as "Jy/sr2".
     *
     * \param unit1 (Gyoto::Units::Unit) first unit
     * \param unit2 (Gyoto::Units::Unit) second unit
     *
     * \return bool, True if it is possible to convert between the two
     * units, 0 otherwise.
     */
    bool areConvertible(const Unit &unit1, const Unit &unit2);
#   endif
  }
}

#ifdef HAVE_UDUNITS
class Gyoto::Units::Unit : public Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Units::Unit>;
  friend class Gyoto::Units::Converter;
 private:
  ut_unit * unit_; ///< the underlying ut_unit (from udunits)
  std::string kind_; ///< the string used to instantiate this unit
 public:

  /**
   * \brief Build Unit described by string
   *
   * Throws a Gyoto::Error if anything goes wrong.
   *
   * \param unit string description of the unit, e.g. "mJy/sr2" or "sunmass".
   */
  Unit(const std::string &unit);

  /**
   * \brief Build Unit described by C string
   *
   * Throws a Gyoto::Error if anything goes wrong.
   *
   * \param unit char const * const description of the unit,
   *        e.g. "mJy/sr2" or "sunmass".
   */
  Unit(char const * const unit);

  /**
   * \brief Destructor
   *
   * Frees unit_.
   */
  ~Unit();

  /**
   * \brief Convert to Unit
   *
   * \param val double to convert
   * \param from_unit Unit from which to convert
   * 
   * \return value converted to unit_.
   */
  double To (double val, const Unit &from_unit);

  /**
   * \brief Convert from Unit
   *
   * \param val double to convert
   * \param to_unit Unit to which to convert
   * 
   * \return value converted to "to_unit".
   */
  double From (double val, const Unit &to_unit);

  /**
   * \brief Cast to string
   *
   * \return kind_
   */
  operator std::string() const ;

  /**
   * \brief Cast to ut_unit*
   *
   * \return unit_
   */
  operator ut_unit*() const ;
};

class Gyoto::Units::Converter : public Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Units::Converter>;
 private:
  cv_converter * converter_; ///< Underlying ut_converter object from udunits

 public:
  Converter(); ///< Construct trivial Converter (Converter()(x)==x)
  Converter(const Gyoto::Units::Unit& from,
	    const Gyoto::Units::Unit& to);
  ///< Construct Converter from two Unit
  ~Converter();
  ///< Destruct converter, freeing converter_

  void reset(); ///< Reset to trivial Converter (Converter()(x)==x)
  void reset(const Gyoto::Units::Unit& from,
	const Gyoto::Units::Unit& to);
  ///< Reset to converter from "from" to "to"

  /**
   * \brief Actually convert data
   *
   * The entire Gyoto::Units::Converter class is there just for this
   * operator, which converts value from unit "from" to unit "to"
   * where "from" and "to" are the two Units passed to the constructor
   * Gyoto::Units::Converter::Converter(const Gyoto::Units::Unit&
   * from, const Gyoto::Units::Unit& to).
   *
   * \param value double expressed in Unit from
   *
   * \return converted value expressed in Unit to
   */
  double operator()(double value) const ;
};

#endif

#endif
