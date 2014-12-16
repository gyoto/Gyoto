#include <Python.h>
#include <iostream>

#define GYOTO_NO_DEPRECATED
#include "GyotoConfig.h"
#include "GyotoFactory.h"
#include "GyotoMetric.h"

Gyoto::SmartPointer<Gyoto::Metric::Generic> pyGyotoMetric(std::string const&s) {
  return Gyoto::Metric::getSubcontractor(s.c_str())(NULL);
}
Gyoto::SmartPointer<Gyoto::Astrobj::Generic> pyGyotoAstrobj(std::string const&s) {
  return Gyoto::Astrobj::getSubcontractor(s.c_str())(NULL);
}
Gyoto::SmartPointer<Gyoto::Spectrum::Generic> pyGyotoSpectrum(std::string const&s) {
  return Gyoto::Spectrum::getSubcontractor(s.c_str())(NULL);
}
Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> pyGyotoSpectrometer(std::string const&s) {
  return Gyoto::Spectrometer::getSubcontractor(s.c_str())(NULL);
}
