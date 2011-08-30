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

#include <cstddef>
#include <GyotoSmartPointer.h>
#include <GyotoError.h>
#include <string>
#include <cstring>

Gyoto::SmartPointee::SmartPointee() : refCount (0) { }
Gyoto::SmartPointee::SmartPointee(const SmartPointee&) : refCount (0) { }
void Gyoto::SmartPointee::incRefCount () { refCount++; }
int Gyoto::SmartPointee::decRefCount () { return --refCount; }
int Gyoto::SmartPointee::getRefCount() { return refCount; }
