/*
    Copyright 2011-2025 Thibaut Paumard

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
#include <glob.h>
#include "GyotoSmartPointer.h"

/**
 * \file GyotoRegister.h
 * \brief Gyoto registers
 *
 * Gyoto::Register::Entry instances are used to map kind names to
 * Gyoto::SmartPointee::Subcontractor_t functions used to instantiate
 * objects from XML files through the Gyoto::Factory.
 */

namespace Gyoto {
  /**
   * \namespace Gyoto::Register
   * \brief Gyoto registers
   * 
   * Gyoto::Register::Entry instances are used to map kind names to
   * Gyoto::SmartPointee::Subcontractor_t functions used to
   * instantiate objects from XML files through the Gyoto::Factory.
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
     * plug-in prepended with "nofail:" is not fatal but issues a
     * warning. Failing to load a plug-in prepended with "nowarn:" is
     * silently ignored.
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
   * \param[in] plugname C string Plug-in name.
   *
   * \param[in] nofail int Unless nofail evals to true, the inability
   *            to find a plug-in or to run the initialization
   *            function inside it throws an Gyoto::Error. If nofail
   *            is 2 or more, such conditions are silently ignored. If
   *            nofail is one, those conditions trigger a warning.
   *
   * \return void* handle to the dlopen'ed plug-in.
   */
  void * loadPlugin(   char const * const plugname, int nofail = 0);

  /**
   * \brief Check whether a given plug-in has already been loaded
   *
   * \param[in] plugname std::string Plug-in name.
   */
  bool havePlugin(std::string plugname);

  /**
   * \brief Load a plugin by name, only if not loaded yet
   *
   * \param[in] plugname std::string Plug-in name.
   *
   * \param[in] nofail int Unless nofail evals to true, the inability
   *            to find a plug-in or to run the initialization
   *            function inside it throws an Gyoto::Error. If nofail
   *            is 2 or more, such conditions are silently ignored. If
   *            nofail is one, those conditions trigger a warning.
   */
  void requirePlugin(std::string plugname, int nofail = 0);

  /**
   * \brief Get a copy of the plug-in path
   */
  std::vector<std::string> pluginPath();
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
  const std::string plugin_;
    ///< Plug-in from which this Entry was loaded
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
   * Search through the register (i.e. this Entry and its descendants)
   * for an Entry matching \p name and return the corresponding
   * subcontractor. If \p plugin is specified, only a subcontractor
   * matching both \p name and \p plugin will be returned. If \p
   * plugin is the empty string, then the first subcontractor matching
   * \p name will be returned, and the name of the plug-in it belongs
   * to will be returned in \p plugin upon output.
   *
   * \note
   * Gyoto::Register::Entry::getSubcontractor() has two notable
   * differences compared to the getSubcontractor() methods defined
   * using #GYOTO_GETSUBCONTRACTOR
   * (e.g. Gyoto::Metric::getSubcontractor()):
   * - it will not load \p plugin;
   * - \p plugin is a scalar string and not a vector of strings.
   *
   * \param[in] name Name of the kind to look for.
   * \param[inout] plugin e.g. "stdplug".
   * \param[in] errmode 1 if getSubContractor() should return NULL upon
   * failure. Else a Gyoto::Error is thrown.
   * \return Pointer to subcontractor function.
   */
  Gyoto::SmartPointee::Subcontractor_t*
    getSubcontractor(std::string name, std::string &plugin, int errmode=0);

  /**
   * \brief Get name
   */
  std::string name();

  /**
   * \brief Get plugin
   */
  std::string plugin();

  /**
   * \brief Get next
   */
  Register::Entry* next();
};

/**
 * \def GYOTO_GETSUBCONTRACTOR(space)
 *
 * \brief Defines the getSubcontractor() function for namespace \p
 * space.
 *
 * This macro is called for instance to define
 * Gyoto::Metric::getSubcontractor(). A function defined this way will:
 * - load any plug-in listed in \p plugin using requirePlugin();
 * - if \p plugin is empty, look for a subcontractor matching \p name
 *   in the relevant register,
 *   - if a match is found, insert the name of the corresponding plug-in
 *     in \p plugin and return a pointer to the matching
 *     subcontractor;
 *   - else,
 *     - if \p errmode is 1, return NULL;
 *     - if \p errmode is 0, throw a Gyoto::Error;
 * - else, look for a subcontractor matching \p name and one element of
 *   \p plugin in the relevant register,
 *   - if a match is found, return the subcontractor;
 *   - else,
 *     - if \p errmode is 1, return NULL;
 *     - if \p errmode is 0, throw a Gyoto::Error.
 *
 * \param[in] space a Gyoto namespace such as Metric, Astrobj,
 * Spectrum or Spectrometer.
 */
#define GYOTO_GETSUBCONTRACTOR(space)					\
  Gyoto::space::Subcontractor_t*					\
  Gyoto::space::getSubcontractor(std::string name,			\
				 std::vector<std::string> &plugin,	\
				 int errmode) {				\
    std::vector<std::string> mandatory;					\
    std::vector<std::string> fallback;					\
    std::vector<std::string> plug_path(Gyoto::pluginPath());		\
    GYOTO_DEBUG << "loading non-fallback plug-ins..." << std::endl;	\
    for (std::vector<std::string>::iterator plg = plugin.begin();	\
	 plg != plugin.end();						\
	 ++plg) {							\
      GYOTO_DEBUG_EXPR(*plg);						\
      if (plg -> rfind("fallback:", 0) != 0) {				\
	Gyoto::requirePlugin(*plg);					\
	mandatory.emplace_back(*plg);					\
      }									\
      else fallback.emplace_back(plg->substr(9));			\
    }									\
    GYOTO_DEBUG <<							\
      "loading fallback plug-ins until the Register is not empty"	\
		<< std::endl;						\
    for (std::vector<std::string>::iterator plg = fallback.begin();	\
	 plg != fallback.end() && !Gyoto::space::Register_;		\
	 ++plg) {							\
      GYOTO_DEBUG_EXPR(*plg);						\
      Gyoto::requirePlugin(*plg, 1);					\
    }									\
    for (std::vector<std::string>::iterator plg = fallback.begin();	\
	 plg != fallback.end() && !Gyoto::space::Register_;		\
	 ++plg) {							\
      GYOTO_DEBUG_EXPR(*plg);						\
      for (std::vector<std::string>::iterator path = plug_path.begin();	\
	   path != plug_path.end();					\
	   ++path) {							\
	std::string pattern = (*path + "libgyoto-" + *plg)		\
		     + "." GYOTO_PLUGIN_SFX;				\
	std::vector<std::string> files =				\
	  Gyoto::glob(pattern, GLOB_NOCHECK | GLOB_NOSORT |		\
		      GLOB_TILDE | GLOB_BRACE);				\
	for (std::vector<std::string>::iterator file = files.begin();	\
	     file != files.end();					\
	     ++file) {							\
	  GYOTO_DEBUG << "Trying " << *file << std::endl;		\
	  Gyoto::requirePlugin(*file, 1);				\
	}								\
      }									\
    }									\
    if (!Gyoto::space::Register_)					\
      throwError("No " GYOTO_STRINGIFY(space) " kind registered!");	\
    Subcontractor_t* sctr= NULL;					\
    GYOTO_DEBUG << "looking for " << name				\
                << " in non-fallback plug-ins..."  << std::endl;	\
    for (std::vector<std::string>::iterator plg = mandatory.begin();	\
	 plg != mandatory.end();					\
	 ++plg) {							\
      GYOTO_DEBUG_EXPR(*plg);						\
      sctr=(Subcontractor_t*)Gyoto::space::Register_			\
	-> getSubcontractor(name, *plg, 1);				\
      if (sctr) {							\
	GYOTO_DEBUG << "found " << name << " in plug-in "		\
		    << *plg << std::endl;				\
	return sctr;							\
      }									\
    }									\
    if (!mandatory.size()) {						\
      GYOTO_DEBUG << "looking for " << name				\
		  << " in registered plug-ins..."  << std::endl;	\
      std::string plg("");						\
      sctr = (Subcontractor_t*)Gyoto::space::Register_			\
	-> getSubcontractor(name, plg, 1);				\
      if (sctr){							\
	GYOTO_DEBUG << "found " << name << " in plug-in "		\
		    << plg << std::endl;				\
	plugin.emplace_back(plg);					\
	return sctr;							\
      } else if (!fallback.size() && !errmode)				\
	throwError ("Kind not found in any plug-in: "+name);		\
    }									\
    GYOTO_DEBUG << "looking for " << name				\
		<< " in fallback plug-ins..."  << std::endl;		\
    for (std::vector<std::string>::iterator plg = fallback.begin();	\
	 plg != fallback.end();						\
	 ++plg) {							\
      GYOTO_DEBUG_EXPR(*plg);						\
      Gyoto::requirePlugin(*plg, 1);					\
      sctr =								\
	(Subcontractor_t*)Gyoto::space::Register_			\
	-> getSubcontractor(name, *plg, 1);				\
      if (sctr) {							\
	GYOTO_DEBUG << "found " << name << " in plug-in "		\
		    << *plg << std::endl;				\
	return sctr;							\
      }									\
    }									\
    GYOTO_DEBUG << "looking for " << name				\
		<< " in fallback plug-ins, "				\
		<< "allowing for glob expansion..."  << std::endl;	\
    for (std::vector<std::string>::iterator plg = fallback.begin();	\
	 plg != fallback.end();						\
	 ++plg) {							\
      GYOTO_DEBUG_EXPR(*plg);						\
      for (std::vector<std::string>::iterator path = plug_path.begin();	\
	   path != plug_path.end();					\
	   ++path) {							\
	std::string pattern = (*path + "libgyoto-" + *plg)		\
		     + "." GYOTO_PLUGIN_SFX;				\
	std::vector<std::string> files =				\
	  Gyoto::glob(pattern, GLOB_NOCHECK | GLOB_NOSORT |		\
		      GLOB_TILDE | GLOB_BRACE);				\
	for (std::vector<std::string>::iterator file = files.begin();	\
	     file != files.end();					\
	     ++file) {							\
	  GYOTO_DEBUG << "Trying " << *file << std::endl;		\
	  Gyoto::requirePlugin(*file, 1);				\
	  sctr =							\
	    (Subcontractor_t*)Gyoto::space::Register_			\
	    -> getSubcontractor(name, *file, 1);			\
	  if (sctr) {							\
	    GYOTO_DEBUG << "found " << name << " in plug-in "		\
			<< *file << std::endl;				\
	    return sctr;						\
	  }								\
	}								\
      }									\
    }									\
    GYOTO_DEBUG << name << " not found anywhere, error?"		\
                << std::endl;						\
    if (!errmode)							\
      throwError("Kind not found in the specified plug-ins: "+name);	\
    return sctr;							\
  }

#endif
