/*
    Copyright 2014-2015 Thibaut Paumard

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
