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

#include <GyotoThinInfiniteDiskKS.h>
#include "ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

extern "C" {
  // THININFINITEDISK CLASS
  // Constructor
  void
  Y_gyoto_ThinInfiniteDiskKS(int n)
  {
    if (n!=2) y_error("gyoto_ThinInfiniteDisk_new takes exactly 2 arguments");
    SmartPointer<Metric::Generic> *gg = yget_Metric(1);
    if ((*gg)->getKind() != "KerrKS") y_error("Metric must be KerrKS");

    SmartPointer<Astrobj::Generic> *astrobj=ypush_Astrobj();
    try {
      *astrobj=new ThinInfiniteDiskKS(*gg);
    } YGYOTO_STD_CATCH ;
    //strcpy(obj->type,"ThinInfiniteDisk");
  }
}
