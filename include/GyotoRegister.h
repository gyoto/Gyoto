/*
    Copyright 2011-2015 Thibaut Paumard

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

/**
 * \file GyotoRegister.h
 * \brief Gyoto registers
 *
 * Gyoto::Register::Entry instances are used to map kind names to
 * Gyoto::SmartPointee::Subcontractor_t functions used to instanciate
 * objects from XML files through the Gyoto::Factory.
 */

namespace Gyoto {
  /**
   * \namespace Gyoto::Register
   * \brief Gyoto registers
   * 
   * Gyoto::Register::Entry instances are used to map kind names to
   * Gyoto::SmartPointee::Subcontractor_t functions used to
   * instanciate objects from XML files through the Gyoto::Factory.
   */
  namespace Register {

    /* Documented below */
    class Entry;

    /**
     * \brief Initialise the various registers
     *
     * Normally called once at application start-up, Register::init()
     * initiaizes the registers, loads the plug-ins, and fills the
     * registers as appropriate.
     *
     * \param pluglist Coma-separated list of plug-ins to load. If
     * NULL, default to the environment variable GYOTO_PLUGINS, if it
     * exists. Else use GYOTO_DEFAULT_PLUGINS. Failing to load a
     * plug-in prepended with "nofail:" is not fatal.
     */
    void init( char const * pluglist = NULL );

    /**
     * \brief List the various registers
     */
    void list();
  }

  /**
   * \brief Load a plugin by name
   *
   * Uses dlopen to load the file libgyoto-&lt;plugname&gt;.so, looks for
   * the function __Gyoto&lt;plugname&gt;Init inside it and run it.
   * Plug-ins must be located in the runtime link search path, or in
   * GYOTO_PKGLIBDIR, or in
   * GYOTO_PKGLIBDIR/GYOTO_SOVERS/.
   *
   * \param plugname Plug-in name.
   * \param nofail Unless nofail evals to true, the inability to find
   * a plug-in or to run the initialization function inside it throws
   * an Gyoto::Error.
   */
  void loadPlugin(   char const * const plugname, int nofail = 0);
}

/**
 * \brief Entry in a register (or a full register)
 *
 * A register is actually a chained list of Register::Entry
 * instances.
 */
class Gyoto::Register::Entry {
  /**
   * \brief List the various registers
   */
  friend void Gyoto::Register::list();
protected:
  std::string name_;
    ///< Kind name for the entry, as found in the "kind" XML attribute
  Gyoto::SmartPointee::Subcontractor_t* subcontractor_;
    ///< Pointer to the Gyoto::SmartPointee::Subcontractor_t function that produces an object of this kind
  Register::Entry* next_;
    ///< Next entry in the register, or NULL
public:
  /**
   * \brief Constructor
   */
  Entry(std::string name,
		Gyoto::SmartPointee::Subcontractor_t* subcontractor,
		Entry* next);
  ~Entry(); ///< Destructor

  /**
   * \brief Get subcontractor for a given name
   * 
   * Search through the register for an Entry matching name and return
   * the corresponding subcontractor.
   *
   * \param name Name of the kind to look for.
   * \param errmode 1 if getSubContractor() should return NULL upon
   * failure. Else a Gyoto::Error is thrown.
   * \return Pointer to subcontractor function.
   */
  Gyoto::SmartPointee::Subcontractor_t*
    getSubcontractor(std::string name, int errmode=0);
};

#endif
