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

#include <GyotoError.h>
#include <iostream>
using namespace Gyoto;
using namespace std;

//Error::Error( const char* m ) : message(m), errcode(0) { }
Error::Error( const std::string m ) : message(m), errcode(0) { }

void Error::Report() const { cerr << message << endl; }

int Error::getErrcode() const { return errcode ; }

//char const * const Error::get_message() const { return message; }
std::string Error::get_message() const { return message; }

static GyotoErrorHandler_t * GyotoErrorHandler = NULL;

void Gyoto::setErrorHandler( GyotoErrorHandler_t* handler )
{ GyotoErrorHandler = handler ; }

void Gyoto::throwError( const std::string m ) {
  if (GyotoErrorHandler) (*GyotoErrorHandler)(m.c_str());
  else throw Error(m);
}
