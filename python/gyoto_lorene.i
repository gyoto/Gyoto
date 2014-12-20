%module(docstring="The Gyoto Lorene plug-in") gyoto_lorene
%import gyoto.i

%{

#define GYOTO_NO_DEPRECATED
#include "GyotoWorldline.h"
#include "GyotoPhoton.h"
#include "GyotoUniformSphere.h"
#include "GyotoPhoton.h"
#include "GyotoScenery.h"
#include "GyotoSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoValue.h"
#include "GyotoThinDisk.h"

#include "GyotoRotStar3_1.h"
#include "GyotoNumericalMetricLorene.h"

%}

%array_class(double, array_double)
%array_class(double, array_unsigned_long)

GyotoSmPtrClassDerived(Metric, RotStar3_1)
GyotoSmPtrClassDerived(Metric, NumericalMetricLorene)

