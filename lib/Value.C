#include "GyotoValue.h"
#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoSpectrum.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"
#include "GyotoProperty.h"
#include <iostream>
using namespace Gyoto ;
using namespace std ;

/// Value

Value::Value() {}
Value::~Value() {}

#define ___local_stuff(t, t_t, T)			\
  Value::Value(t val) : type(Property::t_t), T(val){}	\
  Value::operator t() const {				\
  if (type!=Property::t_t)				\
    throwError("This Value does not hold a " #t);	\
  return T;						\
  }

Value::Value(long val) : type(Property::long_t), Long(val){}
Value::Value(unsigned long val) : type(Property::unsigned_long_t), ULong(val){}
Value::operator long() const {
  switch (type) {
  case Property::long_t:
    return Long;
  case Property::unsigned_long_t:
    return long(ULong);
  default:
    throwError("This Value does not hold a long (or unsigned long)");
  }
  return 0;
}

Value::operator unsigned long() const {
  switch (type) {
  case Property::long_t:
    return (unsigned long)(Long);
  case Property::unsigned_long_t:
    return ULong;
  default:
    throwError("This Value does not hold a long (or unsigned long)");
  }
  return 0;
}

___local_stuff(double, double_t, Double)
___local_stuff(bool, bool_t, Bool)
___local_stuff(std::string, string_t, String)
___local_stuff(std::vector<double>, vector_double_t, VDouble)
___local_stuff(std::vector<unsigned long>, vector_unsigned_long_t, VULong)
___local_stuff(Gyoto::SmartPointer<Gyoto::Metric::Generic>, metric_t, Metric)
___local_stuff(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>, astrobj_t, Astrobj)
___local_stuff(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>, spectrum_t, Spectrum)
___local_stuff(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>, spectrometer_t, Spectrometer)
___local_stuff(Gyoto::SmartPointer<Gyoto::Screen>, screen_t, Screen)

Value& Value::operator=(Value const &right) {
# define ___local_case(member) member = right.member
  ___local_case(type);
  ___local_case(Double);
  ___local_case(Bool);
  ___local_case(Long);
  ___local_case(ULong);
  ___local_case(String);
  ___local_case(VDouble);
  ___local_case(VULong);
  ___local_case(Metric);
  ___local_case(Astrobj);
  ___local_case(Spectrum);
  ___local_case(Spectrometer);
  ___local_case(Screen);
  return *this;
# undef ___local_case
}
