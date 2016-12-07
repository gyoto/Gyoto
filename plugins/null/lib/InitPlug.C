/*
    Copyright © 2016 Thibaut Paumard

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

// #include all headers necessary to register subcontractors below
#include "GyotoNull.h"

using namespace Gyoto;

// Rename replace "null" with actual name of plug-in in the init
// function name below:
extern "C" void __GyotonullInit() {
  // Register subcontractors for all new classes
  Astrobj::Register("Null", &Astrobj::Subcontractor<Astrobj::Null>);
}
