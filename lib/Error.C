/*
    Copyright 2011, 2013 Thibaut Paumard

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
#include <cstdlib>
using namespace Gyoto;
using namespace std;

Error::Error( const std::string m ) : message(m), errcode(EXIT_FAILURE) { }

Error::Error( const Gyoto::Error &o): message(o.message), errcode(o.errcode) {}

void Error::Report() const { cerr << message << endl; }

int Error::getErrcode() const { return errcode ; }

//char const * const Error::get_message() const { return message; }
std::string Error::get_message() const { return message; }

static Gyoto::Error::Handler_t * GyotoErrorHandler = NULL;

void Gyoto::Error::setHandler( Gyoto::Error::Handler_t* handler )
{ GyotoErrorHandler = handler ; }

void Gyoto::throwError( const std::string m ) {
  if (GyotoErrorHandler) (*GyotoErrorHandler)(Error(m));
  else throw Error(m);
}

Gyoto::Error::operator const char * () const { return message.c_str(); }
