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

double Gyoto::Units::ToMeters(double val, string unit,
			       const SmartPointer<Metric::Generic> &gg) {
  if ((unit=="") || (unit=="m")) return val ;
  if (unit=="geometrical") {
    if (gg) return val * gg -> unitLength();
    else throwError("Metric required for geometrical -> meter conversion");
  }
# ifdef HAVE_UDUNITS
  return Meter->To(val, unit); 
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
  throwError(ss.str());
  return 0;
# endif
}

double Gyoto::Units::FromMeters(double val, string unit,
			       const SmartPointer<Metric::Generic> &gg) {
  if ((unit=="") || (unit=="m")) return val ;
  if (unit=="geometrical") {
    if (gg) return val / gg -> unitLength();
    else throwError("Metric required for meter -> geometrical conversion");
  }
# ifdef HAVE_UDUNITS
  return Meter->From(val, unit); 
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
  throwError(ss.str());
  return 0;
# endif
}

double Gyoto::Units::ToSeconds(double val, const string &unit,
			       const SmartPointer<Metric::Generic> &gg) {
  if (unit=="" || unit=="s") ;
  else if (unit=="geometrical_time" || unit=="geometrical") {
    if (unit=="geometrical")
      GYOTO_WARNING << "Please use \"geometrical_time\" instead of "
	"\"geometrical\" for time unit, \"geometrical\" in this context "
	"is deprecated and will be removed soon";
    if (gg) val *= gg -> unitLength() / GYOTO_C ;
    else
      throwError("Metric required for geometrical_time -> second conversion");
  }
# ifdef HAVE_UDUNITS
  else val = Converter(unit, "s")(val);
# else
  else if (unit=="min") val *= 60. ;
  else if (unit=="h") val *= 3600. ;
  else if (unit=="d") val *= 86400. ;
  else if (unit=="yr") val *= 3.15576e+07;
  else {
    stringstream ss;
    ss << "Screen::setTime(): unknown unit \"" << unit << "\". Accepted units: "
       << "[s] geometrical min h d y";
    throwError (ss.str());
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
      throwError("Metric required for second -> geometrical_time conversion");
  }
# ifdef HAVE_UDUNITS
  else val = Converter("s", unit)(val);
# else
  else if (unit=="min") val /= 60. ;
  else if (unit=="h") val /= 3600. ;
  else if (unit=="d") val /= 86400. ;
  else if (unit=="yr") val /= 3.15576e+07;
  else {
    stringstream ss;
    ss << "Screen::setTime(): unknown unit \"" << unit << "\". Accepted units: "
       << "[s] geometrical min h d y";
    throwError (ss.str());
  }
# endif
  return val;
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

double Gyoto::Units::FromGeometrical(double val, string unit,
				     SmartPointer<Metric::Generic> gg) {
  if (unit == "" || unit == "geometrical") return val;
  if (!gg) throwError("Need Metric to convert from geometrical units");
  return FromMeters(val * gg->unitLength(), unit);
}

bool Gyoto::Units::areConvertible(std::string unit1, std::string unit2) {
# ifdef HAVE_UDUNITS
  return ut_are_convertible(Unit(unit1), Unit(unit2));
# else
  throwError("Gyoto::Units::areConvertible unimplemented without udunits");
  return 0;
# endif
}
