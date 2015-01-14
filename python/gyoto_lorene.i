%module(docstring="The Gyoto Lorene plug-in") gyoto_lorene
%import gyoto.i

%{

#define GYOTO_NO_DEPRECATED
#include "GyotoConfig.h"

#include "GyotoWorldline.h"
#include "GyotoPhoton.h"
#include "GyotoUniformSphere.h"
#include "GyotoPhoton.h"
#include "GyotoScenery.h"
#include "GyotoSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoUniformSpectrometer.h"
#include "GyotoValue.h"
#include "GyotoThinDisk.h"

#include "GyotoRotStar3_1.h"
#include "GyotoNumericalMetricLorene.h"

  using namespace Gyoto;

%}

%array_class(double, array_double)
%array_class(unsigned long, array_unsigned_long)
%array_class(size_t, array_size_t)

GyotoSmPtrClassDerived(Metric, RotStar3_1)
GyotoSmPtrClassDerived(Metric, NumericalMetricLorene)

// Workaround cvar bug in Swig which makes help(module) fail:
%inline {
  namespace GyotoLorene {
    extern int __class__=0;
  }
}
