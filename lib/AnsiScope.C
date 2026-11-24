/*
    Copyright 2026 Thibaut Paumard & Frédéric Vincent

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

#include <iostream>

#include "GyotoAnsi.h"

Gyoto::AnsiScope::AnsiScope(std::ostream& s)
  : stream_(s), reset_(true)
{}

Gyoto::AnsiScope::~AnsiScope() {
  if (reset_) stream_ << "\033[0m"; // ANSI code for RESET
}

void Gyoto::AnsiScope::no_reset() { reset_ = false; }

Gyoto::AnsiScope Gyoto::AnsiScope::cout() {
  return Gyoto::AnsiScope(std::cout);
}

Gyoto::AnsiScope Gyoto::AnsiScope::cerr() {
  return Gyoto::AnsiScope(std::cerr);
}

void Gyoto::AnsiScope::append(std::string txt) {
  stream_ << txt;
}

Gyoto::AnsiScope& Gyoto::AnsiScope::operator<<
(std::ostream& (*manip)(std::ostream&)) {
  manip(stream_);
  return *this;
}

void Gyoto::AnsiScope::flush() {stream_.flush();}
