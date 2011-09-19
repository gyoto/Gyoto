/*
    Copyright 2011 Thibaut Paumard

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

#include "ygyoto.h"
#include <GyotoAnalyticMetric.h>
#include "yapi.h"

using namespace Gyoto;

extern "C" {
  // ANALYTICMETRIC CLASS
  // Constructor
  void
  Y_gyoto_AnalyticMetric_new(int n)
  {
    double mm=1.0;
    if (n>1) y_error("gyoto_Kerr_new takes at most 1 arguments");
    if (n==1) mm=ygets_d(0);

    SmartPointer<Metric::Generic> *gg = ypush_Metric();
    try { *gg = new AnalyticMetric(mm); }
    YGYOTO_STD_CATCH ;

    //    strcpy(obj->type,"AnalyticMetric");
  }

}
