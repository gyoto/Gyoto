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
%module(docstring="The Gyoto Lorene plug-in", package="gyoto") lorene
%import gyoto.i

%{

#define GYOTO_NO_DEPRECATED
#include "gyoto_swig.h"

#include "GyotoRotStar3_1.h"
#include "GyotoNumericalMetricLorene.h"

#include "GyotoNeutronStar.h"
#include "GyotoNeutronStarAnalyticEmission.h"
#include "GyotoNeutronStarModelAtmosphere.h"

using namespace Gyoto;

%}

%array_class(double, array_double)
%array_class(unsigned long, array_unsigned_long)
%array_class(size_t, array_size_t)

// This will be called upon extension initialization
%init {
#ifdef SWIGPYTHON
  import_array();
#endif
 }

GyotoSmPtrClassDerivedMetric(RotStar3_1)
GyotoSmPtrClassDerivedMetric(NumericalMetricLorene)

GyotoSmPtrClassDerived(Astrobj, NeutronStar)
GyotoSmPtrClassDerived(Astrobj, NeutronStarAnalyticEmission)
GyotoSmPtrClassDerived(Astrobj, NeutronStarModelAtmosphere)

// Workaround cvar bug in Swig which makes help(module) fail:
%inline {
  namespace GyotoLorene {
    extern int __class__;
  }
  int GyotoLorene::__class__ = 0;
}
