/**
 * \file GyotoWIP.h
 * \brief Work in progress class
 * 
 * A lot of work in Gyoto is still in progress. When created, every
 * classe should be inherit from this class, untill considered stable
 * by the developpers.
 *
 */

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

#ifndef __GyotoWIP_H_
#define __GyotoWIP_H_ 

#include <string>

namespace Gyoto {
    class WIP;
}

/**
 * \class Gyoto::WIP
 * \brief Base class for work in progress
 *
 * The constructors of this class simply issue a warning that the
 * (derived) class is work in progress.
 */
class Gyoto::WIP
{
 public:
  /// Issue a warning
  WIP();
  /// Issue a warning specifying the name of the derived class
  /**
   * If classname is the empty string (""), the warning is not issued.
   * Use this to mark that a class is no more work in progress without
   * breaking the ABI (i.e. in the Gyoto stable branch).
   */
  WIP(std::string classname);
};

#endif
