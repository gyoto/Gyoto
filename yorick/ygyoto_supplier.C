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

#include <ygyoto.h>
#include <yapi.h>

#define JOIN(x, y) JOIN_AGAIN(x, y)
#define JOIN_AGAIN(x, y) x ## y
#define Y___SET_YGYOTO_LOCAL_SUPPLIER JOIN(Y___set_, YGYOTO_LOCAL_SUPPLIER)

YGyotoSupplier_t* YGYOTO_LOCAL_SUPPLIER = 0;

extern "C" {

  void
  Y___SET_YGYOTO_LOCAL_SUPPLIER
  (int argc)
  {
    YGYOTO_LOCAL_SUPPLIER = (YGyotoSupplier_t*) ygets_l(0);
    ypush_nil();
  }

}
