#include "GyotoValue.h"
#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoSpectrum.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"
#include <iostream>
using namespace Gyoto ;
using namespace std ;

/// Value

Value::Value() {}
Value::~Value() {}

Value::Value(double val) : Double(val) {}
Value::operator double() const {return Double;}

Value::Value(bool val) : Bool(val) {}
Value::operator bool() const {return Bool;}

Value::Value(long val) : Long(val) {}
Value::operator long() const {return Long;}

Value::Value(unsigned long val) : ULong(val) {}
Value::operator unsigned long() const {return ULong;}

Value::Value(std::string val) : String(val) {}
Value::operator std::string() const {return String;}

Value::Value(std::vector<double> val) : VDouble(val) {}
Value::operator std::vector<double>() const {return VDouble;}

Value::Value(Gyoto::SmartPointer<Gyoto::Metric::Generic> p)
  : Metric(p) {}
Value::operator Gyoto::SmartPointer<Gyoto::Metric::Generic>()
{ return Metric; }

Value::Value(Gyoto::SmartPointer<Gyoto::Astrobj::Generic> p)
  : Astrobj(p) {}
Value::operator Gyoto::SmartPointer<Gyoto::Astrobj::Generic>()
{ return Astrobj; }

Value::Value(Gyoto::SmartPointer<Gyoto::Spectrum::Generic> p)
  : Spectrum(p) {}
Value::operator Gyoto::SmartPointer<Gyoto::Spectrum::Generic>()
{ return Spectrum; }

Value::Value(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> p)
  : Spectrometer(p) {}
Value::operator Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>()
{ return Spectrometer; }

Value::Value(Gyoto::SmartPointer<Gyoto::Screen> p)
  : Screen(p) {}
Value::operator Gyoto::SmartPointer<Gyoto::Screen>()
{ return Screen; }

Value& Value::operator=(Value const &right) {
# define ___local_case(member) member = right.member
  ___local_case(Double);
  ___local_case(Bool);
  ___local_case(Long);
  ___local_case(ULong);
  ___local_case(String);
  ___local_case(VDouble);
  ___local_case(Metric);
  ___local_case(Astrobj);
  ___local_case(Spectrum);
  ___local_case(Spectrometer);
  ___local_case(Screen);
  return *this;
# undef ___local_case
}
