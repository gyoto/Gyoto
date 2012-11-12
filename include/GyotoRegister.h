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

#ifndef __GyotoRegister_H_
#define __GyotoRegister_H_

#include <string>
#include "GyotoSmartPointer.h"

namespace Gyoto {
#ifdef GYOTO_USE_XERCES
  namespace Register {
    class Entry;
    void init( char const * pluglist = NULL );
    void list();
  }
#endif
  void loadPlugin(   char const * const plugname, int nofail = 0);
}

#ifdef GYOTO_USE_XERCES
class Gyoto::Register::Entry {
  friend void Register::list ();
protected:
  std::string name_;
  Gyoto::SmartPointee::Subcontractor_t* subcontractor_;
  int type_;
  Register::Entry* next_;
public:
  Entry(std::string name,
		Gyoto::SmartPointee::Subcontractor_t* subcontractor,
		Entry* next);
  ~Entry();
  Gyoto::SmartPointee::Subcontractor_t*
    getSubcontractor(std::string, int errmode=0);
};

#endif
#endif
