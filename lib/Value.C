#include "GyotoValue.h"
#include "GyotoMetric.h"
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

Value::Value(std::string val) : String(val) {}
Value::operator std::string() const {return String;}

Value::Value(std::vector<double> val) : VDouble(val) {}
Value::operator std::vector<double>() const {return VDouble;}

Value::Value(Gyoto::SmartPointer<Gyoto::Metric::Generic> p)
  : Metric(p) {cerr << "In Value constructor: Metric==" << Metric() << endl;}
Value::operator Gyoto::SmartPointer<Gyoto::Metric::Generic>()
{ return Metric; }

Value& Value::operator=(Value const &right) {
# define ___local_case(member) member = right.member
  ___local_case(Double);
  ___local_case(Bool);
  ___local_case(Long);
  ___local_case(String);
  ___local_case(VDouble);
  ___local_case(Metric);
  return *this;
# undef ___local_case
}
