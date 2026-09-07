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
#include <iomanip>
#include <sstream>
#include <cstdlib>
using namespace Gyoto;
using namespace std;

Error::Error( const std::string &m, boost::stacktrace::stacktrace const &st) :
  message(m),
  stacktrace(),
  errcode(EXIT_FAILURE) {
  std::ostringstream oss;
  for (size_t i = 0; i < st.size(); ++i) {
    auto frame = st[i];

    oss << "\n";
    oss << setw(std::to_string(st.size()-1).length()) << i << "# ";
    oss << (i == 0 ? GYOTO_ANSI_BOLD GYOTO_ANSI_FG_RED : "")
        << frame;
    oss << GYOTO_ANSI_RESET;

  }
  stacktrace += oss.str();
  fullmessage = stacktrace + "\n\n" GYOTO_ANSI_BOLD GYOTO_ANSI_FG_RED + message;
}

Error::Error( const Gyoto::Error &o):
  message(o.message),
  stacktrace(o.stacktrace),
  fullmessage(o.fullmessage),
  errcode(o.errcode) {}

void Error::Report() const {
  cerr << fullmessage << endl;
}

int Error::getErrcode() const { return errcode ; }

std::string Error::get_message() const { return message; }
std::string Error::get_stacktrace() const { return stacktrace; }

static Gyoto::Error::Handler_t * GyotoErrorHandler = NULL;

void Gyoto::Error::setHandler( Gyoto::Error::Handler_t* handler )
{ GyotoErrorHandler = handler ; }

void Gyoto::throwError(const std::string &m,
		       boost::stacktrace::stacktrace const &trace) {
  if (GyotoErrorHandler) (*GyotoErrorHandler)(Error(m, trace));
  else throw Error(m, trace);
}

Gyoto::Error::operator const char * () const {
  return fullmessage.c_str();
}
