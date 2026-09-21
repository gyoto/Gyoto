/*
    Copyright 2011, 2013, 2026 Thibaut Paumard

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

// GyotoConfig.h must be included before boost/stacktrace.hpp
// It is included indirectly from GyotoError.h which include GyotoDefs.h
#include <boost/stacktrace.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
using namespace Gyoto;
using namespace std;

Error::Error(const std::string &m, size_t skip) :
  message(m),
  stacktrace(),
  errcode(EXIT_FAILURE) {
  boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace();
  std::ostringstream oss;
  oss << GYOTO_ANSI_ERROR_TAG << "Backtrace from error frame:"
      << GYOTO_ANSI_RESET;
  for (size_t i = skip; i < st.size(); ++i) {
    auto frame = st[i];

    oss << "\n";
    oss << setw(std::to_string(st.size()-1-skip).length()) << i-skip << "# ";
    oss << (i == skip ? GYOTO_ANSI_BOLD : "")
        << frame;
    oss << GYOTO_ANSI_RESET;

  }
  #if BOOST_STACKTRACE_USE_BACKTRACE
  GYOTO_DEBUG_THIS_EXPR(BOOST_STACKTRACE_USE_BACKTRACE);
  #elif BOOST_STACKTRACE_USE_ADDR2LINE
  GYOTO_DEBUG_THIS_EXPR(BOOST_STACKTRACE_USE_ADDR2LINE);
  #else
  GYOTO_DEBUG_THIS << "Boost.stacktrace uses basic backend" << endl;
  #endif
  stacktrace += oss.str();
  fullmessage = stacktrace + "\n\n" GYOTO_ANSI_ERROR + message
    + "\n\n" + GYOTO_ANSI_RESET;
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

void Gyoto::throwError(const std::string &m) {
  Error e(m, 2); // skip 2 frames: throwError and Error::Error()
  if (GyotoErrorHandler) (*GyotoErrorHandler)(e);
  else throw e;
}

Gyoto::Error::operator const char * () const {
  return fullmessage.c_str();
}
