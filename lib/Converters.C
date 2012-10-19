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
  ut_map_symbol_to_unit("pc", UT_ASCII, tmpu);
  ut_free(tmpu); tmpu = NULL;

  /* astronomical_unit aliases */
  tmpu = ut_get_unit_by_name(SI, "astronomical_unit");
  ut_map_symbol_to_unit("au", UT_ASCII, tmpu);
  ut_map_symbol_to_unit("AU", UT_ASCII, tmpu);
  ut_map_symbol_to_unit("ua", UT_ASCII, tmpu);
  ut_map_symbol_to_unit("astronomicalunit", UT_ASCII, tmpu);
  ut_free(tmpu); tmpu = NULL;

  /* sunradius */
  tmpu = ut_get_unit_by_name(SI, "meter");
  tmpu2 = ut_scale(GYOTO_SUN_RADIUS, tmpu);
  ut_map_symbol_to_unit("sunradius", UT_ASCII, tmpu);
  ut_free(tmpu2);
  ut_free(tmpu);

  /* ly */
  tmpu = ut_get_unit_by_name(SI, "light_year");
  ut_map_symbol_to_unit("ly", UT_ASCII, tmpu);
  ut_free(tmpu);

  /* sunmass */
  tmpu = ut_get_unit_by_symbol(SI, "kg");
  tmpu2 = ut_scale(GYOTO_SUN_MASS, tmpu);
  ut_map_symbol_to_unit("sunmass", UT_ASCII, tmpu2);
  ut_free(tmpu2);
  ut_free(tmpu);

  /* Jansky */
  tmpu = ut_parse(SI, "W.m-2.Hz-1", UT_ASCII);
  if (!tmpu) throwError("Cannot initialize Jansky");
  tmpu2=ut_scale(1e-26, tmpu);
  ut_map_name_to_unit("Jansky", UT_ASCII, tmpu2);
  ut_map_symbol_to_unit("Jy", UT_ASCII, tmpu2);
  ut_free(tmpu2);
  ut_free(tmpu);

  Meter = new Unit("meter");
  Second = new Unit("second");
  KiloGram = new Unit("kilogram");
#endif  
}

#ifdef HAVE_UDUNITS
/* Unit class */
Unit::Unit(string unit) : unit_(NULL), kind_(unit) {
  unit_ = ut_parse(SI, unit.c_str(), UT_ASCII);
  if (!unit_) throwError("Error initializing Unit");
}

Unit::~Unit() {
  ut_free(unit_);
}

double Unit::To(double val, std::string from_unit) {
  return Converter(from_unit, this)(val);
}

double Unit::From(double val, std::string to_unit) {
  return Converter(this, to_unit)(val);
}

Unit::operator std::string() { return kind_; }
Unit::operator ut_unit*() { return unit_; }

/* Converter */
Converter::Converter(string from, string to) : from_(NULL), to_(NULL), converter_(NULL) {
  from_ = new Unit(from);
  to_ = new Unit(to);
  resetConverter_();
}

Converter::Converter(Gyoto::SmartPointer<Gyoto::Units::Unit> from, std::string to) :
  from_(from), to_(NULL), converter_(NULL)
{
  to_ = new Unit(to);
  resetConverter_();
}

Converter::Converter(std::string from, Gyoto::SmartPointer<Gyoto::Units::Unit> to) :
  from_(NULL), to_(to), converter_(NULL)
{
  from_ = new Unit(from);
  resetConverter_();
}

Converter::Converter(Gyoto::SmartPointer<Gyoto::Units::Unit> from,
		     Gyoto::SmartPointer<Gyoto::Units::Unit> to) :
  from_(from), to_(to), converter_(NULL)
{
  resetConverter_();
}

Converter::~Converter() {
  if (converter_) { cv_free(converter_); converter_ = NULL; }
  from_ = NULL;
  to_ = NULL;
}

void Converter::resetConverter_() {
  if (!ut_are_convertible(*from_, *to_)) {
    ut_free(*from_);
    ut_free(*to_);
    stringstream ss;
    ss << "Unsupported conversion: from \"" << string(*from_) << "\" to " << string(*to_) << "\"";
    throwError(ss.str());
  }
  converter_ = ut_get_converter(*from_, *to_);
}
double Converter::operator()(double val) {
  return cv_convert_double(converter_, val);
}

#endif

double Gyoto::Units::ToMeters(double val, string unit) {
# ifdef HAVE_UDUNITS
  if ((unit=="") || (unit=="m")) return val ;
  return Meter->To(val, unit); 
# else
  if ((unit=="") || (unit=="m")) return val ;
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
  throwError(ss.str());
  return 0;
# endif
}

double Gyoto::Units::ToKilograms(double val, string unit)
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
    throwError(ss.str());
  }
  return val;
# endif
}

double Gyoto::Units::ToGeometrical(double val, string unit,
				   SmartPointer<Metric::Generic> gg) {
  if (unit == "" || unit == "geometrical") return val;
  if (!gg) throwError("Need Metric to convert to geometrical units");
  return ToMeters(val, unit) / gg->unitLength();
}

