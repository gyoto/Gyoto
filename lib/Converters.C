/*
    Copyright 2012-2015 Thibaut Paumard

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

#include "GyotoConverters.h"
#include "GyotoError.h"
#include "GyotoMetric.h"

using namespace Gyoto::Units ;
using namespace std ;

#ifdef HAVE_UDUNITS
static ut_system * SI (NULL);
static Gyoto::SmartPointer<Gyoto::Units::Unit> Meter (NULL);
static Gyoto::SmartPointer<Gyoto::Units::Unit> Second (NULL);
static Gyoto::SmartPointer<Gyoto::Units::Unit> KiloGram (NULL);
#endif

void Gyoto::Units::Init() {
#ifdef HAVE_UDUNITS
  ut_set_error_message_handler(ut_ignore);
  SI =  ut_read_xml(NULL);

  ut_unit * tmpu = NULL, * tmpu2 = NULL;

  /* Map pc for parsec */
  tmpu = ut_get_unit_by_name(SI, "parsec");
  ut_map_symbol_to_unit("pc", UT_UTF8, tmpu);
  ut_free(tmpu); tmpu = NULL;

  /* astronomical_unit aliases */
  tmpu = ut_get_unit_by_name(SI, "astronomical_unit");
  ut_map_symbol_to_unit("au", UT_UTF8, tmpu);
  ut_map_symbol_to_unit("AU", UT_UTF8, tmpu);
  ut_map_symbol_to_unit("ua", UT_UTF8, tmpu);
  ut_map_symbol_to_unit("astronomicalunit", UT_UTF8, tmpu);
  ut_free(tmpu); tmpu = NULL;

  /* as symbol for arcsec */
  tmpu = ut_get_unit_by_name(SI, "arcsec");
  if (!tmpu) GYOTO_ERROR("error initializing arcsec unit");
  ut_map_symbol_to_unit("as", UT_UTF8, tmpu);
  ut_free(tmpu); tmpu = NULL;

  /* sunradius */
  tmpu = ut_get_unit_by_name(SI, "meter");
  tmpu2 = ut_scale(GYOTO_SUN_RADIUS, tmpu);
  ut_map_symbol_to_unit("sunradius", UT_UTF8, tmpu2);
  ut_free(tmpu2);
  ut_free(tmpu);

  /* ly */
  tmpu = ut_get_unit_by_name(SI, "light_year");
  ut_map_symbol_to_unit("ly", UT_UTF8, tmpu);
  ut_free(tmpu);

  /* sunmass */
  tmpu = ut_get_unit_by_symbol(SI, "kg");
  tmpu2 = ut_scale(GYOTO_SUN_MASS, tmpu);
  ut_map_symbol_to_unit("sunmass", UT_UTF8, tmpu2);
  ut_free(tmpu2);
  ut_free(tmpu);

  /* Jansky */
  tmpu = ut_parse(SI, "W.m-2.Hz-1", UT_UTF8);
  if (!tmpu) GYOTO_ERROR("Cannot initialize Jansky");
  tmpu2=ut_scale(1e-26, tmpu);
  ut_map_name_to_unit("Jansky", UT_UTF8, tmpu2);
  ut_map_symbol_to_unit("Jy", UT_UTF8, tmpu2);
  ut_free(tmpu2);
  ut_free(tmpu);

  Meter = new Unit("meter");
  Second = new Unit("second");
  KiloGram = new Unit("kilogram");
#endif  
}

#ifdef HAVE_UDUNITS
/* Unit class */
Unit::Unit(const string &unit) : unit_(NULL), kind_(unit) {
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(kind_);
# endif
  if (kind_!="") {
    unit_ = ut_parse(SI, unit.c_str(), UT_UTF8);
#   ifdef GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_EXPR(unit_);
#   endif
    if (!unit_) GYOTO_ERROR(string("Error initializing Unit: ")+unit);
  }
}

Unit::Unit(char const * const unit) : unit_(NULL), kind_(unit) {
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(kind_);
# endif
  if (kind_!="") {
    unit_ = ut_parse(SI, unit, UT_UTF8);
#   ifdef GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_EXPR(unit_);
#   endif
    if (!unit_) GYOTO_ERROR(string("Error initializing Unit: ")+unit);
  }
}

Unit::~Unit() {
  ut_free(unit_);
}

double Unit::To(double val, const Unit &from_unit) {
  return Converter(from_unit, *this)(val);
}

double Unit::From(double val, const Unit &to_unit) {
  return Converter(*this, to_unit)(val);
}

Unit::operator std::string() const { return kind_; }
Unit::operator ut_unit*() const { return unit_; }

/* Converter */
Converter::Converter() : converter_(cv_get_trivial()) {}

Converter::Converter(const Gyoto::Units::Unit &from,
		     const Gyoto::Units::Unit &to) :
  converter_(NULL)
{
  reset(from, to);
}

void Converter::reset()
{
  if (converter_) cv_free(converter_);
  converter_ = cv_get_trivial();
}

void Converter::reset(const Gyoto::Units::Unit &from,
		      const Gyoto::Units::Unit &to)
{
  if (converter_) { cv_free(converter_); converter_=NULL; }
  if ((ut_unit*)(from) && (ut_unit*)(to)) {
    if (areConvertible(from, to)) converter_ = ut_get_converter(from, to);
    else {
      stringstream ss;
      ss << "Unsupported conversion: from \"" << string(from)
	 << "\" to " << string(to) << "\"";
      GYOTO_ERROR(ss.str());
    }
  } else reset();
}

Converter::~Converter() {
  if (converter_) { cv_free(converter_); converter_ = NULL; }
}

double Converter::operator()(double val) const {
  return cv_convert_double(converter_, val);
}

ut_system * Gyoto::Units::getSystem() { return SI; }

#endif

double Gyoto::Units::ToMeters(double val, const string &unit,
			      const SmartPointer<Metric::Generic> &gg) {
  if (unit=="" || unit=="m") return val ;
  if (unit=="geometrical") {
    if (gg) return val * gg -> unitLength();
    else GYOTO_ERROR("Metric required for geometrical -> meter conversion");
  }
# ifdef HAVE_UDUNITS
  Unit from(unit);
  if (areConvertible(from, *Meter)) return Meter->To(val, from);
  if (areConvertible(from, *Second)) return Second->To(val, from)*GYOTO_C;
  return GYOTO_C/ToHerz(val, unit);
# else
  if (unit=="cm")          return val * 1e-2;
  if (unit=="km")          return val * 1e3;
  if (unit=="sunradius")   return val * GYOTO_SUN_RADIUS;
  if ((unit=="astronomicalunit") ||
  	   (unit=="astronomical_unit") ||
  	   (unit=="AU") ||
  	   (unit=="au") ||
  	   (unit=="ua"))
                           return val * GYOTO_ASTRONOMICAL_UNIT;
  if (unit=="ly")          return val * GYOTO_LIGHT_YEAR;
  if (unit=="pc")          return val * GYOTO_KPC * 1e-3;
  if (unit=="kpc")         return val * GYOTO_KPC;
  if (unit=="Mpc")         return val * GYOTO_KPC * 1e3;
  stringstream ss;
  ss << "Unsupported conversion: \"" << unit;
  GYOTO_ERROR(ss.str());
  return 0;
# endif
}

double Gyoto::Units::FromMeters(double val, const string &unit,
			       const SmartPointer<Metric::Generic> &gg) {
  if ((unit=="") || (unit=="m")) return val ;
  if (unit=="geometrical") {
    if (gg) return val / gg -> unitLength();
    else GYOTO_ERROR("Metric required for meter -> geometrical conversion");
  }
# ifdef HAVE_UDUNITS
  Unit to (unit);
  if (areConvertible(to, *Meter)) return Meter->From(val, to);
  if (areConvertible(to, *Second)) return Second->From(val*GYOTO_C, to);
  return FromHerz(GYOTO_C/val, unit);
# else
  if (unit=="cm")          return val * 1e2;
  if (unit=="km")          return val * 1e-3;
  if (unit=="sunradius")   return val / GYOTO_SUN_RADIUS;
  if ((unit=="astronomicalunit") ||
  	   (unit=="astronomical_unit") ||
  	   (unit=="AU") ||
  	   (unit=="au") ||
  	   (unit=="ua"))
                           return val / GYOTO_ASTRONOMICAL_UNIT;
  if (unit=="ly")          return val / GYOTO_LIGHT_YEAR;
  if (unit=="pc")          return val / GYOTO_KPC * 1e3;
  if (unit=="kpc")         return val / GYOTO_KPC;
  if (unit=="Mpc")         return val / GYOTO_KPC * 1e-3;
  stringstream ss;
  ss << "Unsupported conversion: \"" << unit;
  GYOTO_ERROR(ss.str());
  return 0;
# endif
}

double Gyoto::Units::ToSeconds(double val, const string &unit,
			       const SmartPointer<Metric::Generic> &gg) {
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_EXPR(val);
  GYOTO_DEBUG_EXPR(unit);
  GYOTO_DEBUG_EXPR(gg);
  GYOTO_ENDIF_DEBUG
# endif
  if (unit=="" || unit=="s") ;
  else if (unit=="geometrical_time" || unit=="geometrical") {
    if (unit=="geometrical")
      GYOTO_WARNING << "Please use \"geometrical_time\" instead of "
	"\"geometrical\" for time unit, \"geometrical\" in this context "
	"is deprecated and will be removed soon";
    if (gg) val *= gg -> unitLength() / GYOTO_C ;
    else
      GYOTO_ERROR("Metric required for geometrical_time -> second conversion");
  }
# ifdef HAVE_UDUNITS
  else {
    Unit from (unit);
    if (areConvertible(from, *Meter)) val=Meter->To(val, from)/GYOTO_C;
    else val = Second->To(val, from);
  }
# else
  else if (unit=="min") val *= 60. ;
  else if (unit=="h") val *= 3600. ;
  else if (unit=="d") val *= 86400. ;
  else if (unit=="yr") val *= 3.15576e+07;
  else if (unit=="kyr") val *= 3.15576e+10;
  else {
    stringstream ss;
    ss << "Screen::time(): unknown unit \"" << unit << "\". Accepted units: "
       << "[s] geometrical min h d y";
    GYOTO_ERROR (ss.str());
  }
# endif
  return val;
}

double Gyoto::Units::FromSeconds(double val, const string &unit,
				 const SmartPointer<Metric::Generic> &gg) {
  if (unit=="" || unit=="s") ;
  else if (unit=="geometrical_time" || unit=="geometrical") {
    if (unit=="geometrical")
      GYOTO_WARNING << "Please use \"geometrical_time\" instead of "
	"\"geometrical\" for time unit, \"geometrical\" in this context "
	"is deprecated and will be removed soon";
    if (gg) val *= GYOTO_C / gg -> unitLength() ;
    else
      GYOTO_ERROR("Metric required for second -> geometrical_time conversion");
  }
# ifdef HAVE_UDUNITS
  else {
    Unit to (unit);
    if (areConvertible(to, *Meter)) val=Meter->From(val*GYOTO_C, to);
    else val = Second->From(val, to);
  }
# else
  else if (unit=="min") val /= 60. ;
  else if (unit=="h") val /= 3600. ;
  else if (unit=="d") val /= 86400. ;
  else if (unit=="yr") val /= 3.15576e+07;
  else if (unit=="kyr") val /= 3.15576e+10;
  else {
    stringstream ss;
    ss << "Screen::time(): unknown unit \"" << unit << "\". Accepted units: "
       << "[s] geometrical min h d y";
    GYOTO_ERROR (ss.str());
  }
# endif
  return val;
}

double Gyoto::Units::ToKilograms(double val, const string &unit)
{
# ifdef HAVE_UDUNITS
  if (unit=="" || unit=="kg") return val;
  return KiloGram->To(val, unit);
# else
  if (unit=="" || unit=="kg") ; // do nothing !
  else if (unit=="sunmass") val *= GYOTO_SUN_MASS;
  else if (unit=="g") val*= 1e-3;
  else {
    stringstream ss;
    ss << "Unsupported mass unit: \"" << unit
       << "\". Supported units: [kg] g sunmass";
    GYOTO_ERROR(ss.str());
  }
  return val;
# endif
}

double Gyoto::Units::FromKilograms(double val, const string &unit)
{
  if (unit=="" || unit=="kg") return val;
# ifdef HAVE_UDUNITS
  return KiloGram->From(val, unit);
# else
  else if (unit=="sunmass") val /= GYOTO_SUN_MASS;
  else if (unit=="g") val*= 1e3;
  else {
    stringstream ss;
    ss << "Unsupported mass unit: \"" << unit
       << "\". Supported units: [kg] g sunmass";
    GYOTO_ERROR(ss.str());
  }
  return val;
# endif
}

double Gyoto::Units::ToGeometrical(double val, const string &unit,
				   const SmartPointer<Metric::Generic> &gg) {
  if (unit == "" || unit == "geometrical") return val;
  if (!gg) GYOTO_ERROR("Need Metric to convert to geometrical units");
  return ToMeters(val, unit) / gg->unitLength();
}

double Gyoto::Units::FromGeometrical(double val, const string &unit,
				     const SmartPointer<Metric::Generic> &gg) {
  if (unit == "" || unit == "geometrical") return val;
  if (!gg) GYOTO_ERROR("Need Metric to convert from geometrical units");
  return FromMeters(val * gg->unitLength(), unit);
}

double Gyoto::Units::ToGeometricalTime(double val, const string &unit,
				   const SmartPointer<Metric::Generic> &gg) {
  if (unit == "" || unit == "geometrical_time") return val;
  if (!gg) GYOTO_ERROR("Need Metric to convert to geometrical units");
  return ToSeconds(val, unit) / gg->unitLength() * GYOTO_C;
}

double Gyoto::Units::FromGeometricalTime(double val, const string &unit,
				     const SmartPointer<Metric::Generic> &gg) {
  if (unit == "" || unit == "geometrical_time") return val;
  if (!gg) GYOTO_ERROR("Need Metric to convert from geometrical units");
  return FromSeconds(val * gg->unitLength() / GYOTO_C, unit);
}

double Gyoto::Units::ToHerz(double value, const string &unit) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "converting " << value << string(unit) << "to Herz" << endl;
# endif
  if (unit == "" || unit == "Hz") return value;
# ifdef HAVE_UDUNITS
  Unit from (unit), Hz("Hz");
  if (areConvertible(from, Hz))
    return Hz.To(value, from);
  if (areConvertible(from, *Meter))
    return GYOTO_C/Meter->To(value, from);
  Unit eV("eV");
  if (areConvertible(from, eV))
    return eV.To(value, from) * GYOTO_eV2Hz;
  stringstream ss;
  ss << "Units::ToHerz(): unknown unit \""
     << unit << "\".";
  GYOTO_ERROR (ss.str());
# else
  GYOTO_WARNING_UDUNITS(unit, "Hz");
  return value;
# endif
  GYOTO_ERROR("unit conversion failed");
  return 0.;
}

double Gyoto::Units::FromHerz(double value, const string &unit) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "converting " << value << "Hz to " << string(unit) << endl;
# endif
# ifdef HAVE_UDUNITS
  if (unit == "" || unit == "Hz") return value;
  Unit to(unit), Hz("Hz");
  if (areConvertible(to, Hz))
    return Hz.From(value, to);
  if (areConvertible(to, *Meter))
    return Meter->From(GYOTO_C/value, to);
  Unit eV("eV");
  if (areConvertible(to, eV))
    return eV.From(value / GYOTO_eV2Hz, to);
  stringstream ss;
  ss << "Units::ToHerz(): unknown unit \""
     << unit << "\".";
  GYOTO_ERROR (ss.str());
# else
  GYOTO_WARNING_UDUNITS(unit, "Hz");
# endif
  return 0.;
}

# ifdef HAVE_UDUNITS
bool Gyoto::Units::areConvertible(const Unit &unit1,
				  const Unit &unit2) {
  return ut_are_convertible(unit1, unit2);
}
# endif
