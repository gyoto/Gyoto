/*
    Copyright 2014 Thibaut Paumard

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

#include "GyotoWIP.h"
#include "GyotoUtils.h"

#include <iostream>

using namespace std ; 
using namespace Gyoto ;

WIP::WIP()
{
  GYOTO_WARNING << "**************************************************" << endl;
  GYOTO_WARNING << "Initializing a class marked *work in progress*."    << endl;
  GYOTO_WARNING << "It may be completely buggy."                        << endl;
  GYOTO_WARNING << "Please test extensively and/or contact the author." << endl;
  GYOTO_WARNING << "**************************************************" << endl;
}

WIP::WIP(string classname)
{
  if (classname == "") return;
  GYOTO_WARNING << "**************************************************" << endl;
  GYOTO_WARNING
    << "Class \"" << classname << "\" is marked *work in progress*."    << endl;
  GYOTO_WARNING << "It may be completely buggy."                        << endl;
  GYOTO_WARNING << "Please test extensively and/or contact the author." << endl;
  GYOTO_WARNING << "**************************************************" << endl;
}
