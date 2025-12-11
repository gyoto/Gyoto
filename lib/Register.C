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
#include "GyotoRegister.h"
#include "GyotoUtils.h"
#include "GyotoAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoSpectrum.h"
#include "GyotoSpectrometer.h"
#include "GyotoConverters.h"

#include <dlfcn.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <cstring>

using namespace Gyoto;
using namespace std;

static std::string GyotoRegisterCurrentPlugin ("built-in");

static std::vector<std::string> GyotoRegisteredPlugins;

static std::vector<std::string> GyotoPluginPath;

static const std::vector<std::string> GyotoDefaultPluginPath =
  {
#   if defined GYOTO_LOCALPKGLIBDIR
    GYOTO_LOCALPKGLIBDIR GYOTO_SOVERS "/",
    GYOTO_LOCALPKGLIBDIR,
#   endif
    GYOTO_PKGLIBDIR "/" GYOTO_SOVERS "/",
    GYOTO_PKGLIBDIR "/"
  };


typedef void GyotoInitFcn();

bool Gyoto::havePlugin(std::string name) {
  for (size_t i=0; i < GyotoRegisteredPlugins.size(); ++i)
    if (GyotoRegisteredPlugins[i]==name) return true;
  return false;
}

void Gyoto::requirePlugin(std::string name, int nofail) {
  if (!havePlugin(name)) loadPlugin(name.c_str(), nofail);
}

void * Gyoto::loadPlugin(char const*const nam, int nofail) {
  string name(nam);

  // Determine file name
  string dlfile = "libgyoto-" ;
  dlfile += name ;
  dlfile += "." ;
  dlfile += GYOTO_PLUGIN_SFX ;

  // If nam _is_ a file name, retrieve plug-in name
  string marker="/libgyoto-";
  if (strstr(nam, marker.c_str())) {
    dlfile = name;
    name=name.substr(0, name.rfind("." GYOTO_PLUGIN_SFX));
    name=name.substr(name.rfind(marker)+marker.size());
    GYOTO_DEBUG << name << endl;
  }

  // Prepare name of init function
  string dlfunc = "__Gyoto";
  dlfunc += name ;
  dlfunc += "Init";
  void* handle = NULL;
  GyotoInitFcn* initfcn = NULL;
  char * err = NULL;

  // Keep track of all dlopen() errors
  string errors("Failed loading plug-in '");
  errors += name;
  errors += "'.\n  The following attempts were made:";

  // Try first without path (i.e. in default linker locations)
  GYOTO_DEBUG << "Loading plug-in: " << name <<endl;
  GYOTO_DEBUG << "Trying to dlopen " << dlfile << "...\n";
  handle = dlopen(dlfile.c_str(), RTLD_LAZY | RTLD_GLOBAL);
  if (!handle && (err=dlerror())) {
    errors += "\n  * Error loading ";
    errors += dlfile + ":\n    ";
    errors += err;
  }

  // Then try the various hard-coded locations
  std::vector<std::string>::iterator cur = GyotoPluginPath.begin();
  std::vector<std::string>::iterator end = GyotoPluginPath.end();
  std::string dlfull= dlfile;
  while (!handle && cur != end) {
    dlfull = *cur + dlfile;
    GYOTO_DEBUG << "Trying to dlopen " << dlfull << "...\n";
    handle = dlopen(dlfull.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    ++cur;
    if (!handle && (err=dlerror())) {
      errors += "\n  * Error loading ";
      errors += dlfull + ":\n    ";
      errors += err;
    }
  }

  // Check for success
  if (handle) {
    GYOTO_DEBUG << "Successfully loaded " << dlfull << ".\n";
  } else {
    GYOTO_DEBUG << "Failed loading " << dlfull << ".\n";
    if (nofail) {
      if (nofail == 1 && verbose() >= GYOTO_DEFAULT_VERBOSITY)
	cerr << "WARNING: unable to load optional plug-in " << dlfile << endl;
      return NULL;
    }
    GYOTO_ERROR(errors.c_str());
  }

  // Find and execute init function
  GYOTO_DEBUG << "Searching plug-in init function " << dlfunc << endl;
  initfcn = (GyotoInitFcn*)dlsym(handle, dlfunc.c_str());
  if ( (err=dlerror()) || !initfcn) {
    dlfunc = "__GyotoPluginInit";
    initfcn = (GyotoInitFcn*)dlsym(handle, dlfunc.c_str());
  }
  if ( (err=dlerror()) ) GYOTO_ERROR(err);
  GYOTO_DEBUG << "Calling plug-in init function " << dlfunc << endl;
  std::string tmp_name(GyotoRegisterCurrentPlugin);
  // In case nam is a file name, that's what we wan't to store
  GyotoRegisterCurrentPlugin = nam;
  (*initfcn)();
  GyotoRegisterCurrentPlugin = tmp_name;
  // In case nam is a file name, that's what we wan't to store
  GyotoRegisteredPlugins.insert(GyotoRegisteredPlugins.begin(), nam);
  GYOTO_DEBUG << "Done." << endl;

  return handle;
}

std::vector<std::string> Gyoto::pluginPath() {
  return GyotoPluginPath;
}

void Gyoto::pluginPath(const std::vector<std::string> &v) {
  GyotoPluginPath = v;
}

void Gyoto::Register::init(char const *  cpluglist) {

  // Initialize plug-in path if not already set
  if (!GyotoPluginPath.size()) {
    const char* GYOTO_PLUGIN_PATH = std::getenv("GYOTO_PLUGIN_PATH");
    if (GYOTO_PLUGIN_PATH) {
      std::istringstream iss(GYOTO_PLUGIN_PATH);
      std::string path;
      while (std::getline(iss, path, ':')) {
        if (!path.empty()) {
	  if (path.back() != '/') path += '/';
	  GyotoPluginPath.push_back(path);
	}
      }
    } else GyotoPluginPath = GyotoDefaultPluginPath;
  }

  // Clean registers
  Metric::initRegister();
  Astrobj::initRegister();
  Spectrum::initRegister();
  // This cleans and fills Spectometer::Register_
  Spectrometer::initRegister();

  GyotoRegisteredPlugins.push_back(GyotoRegisterCurrentPlugin);

  // Init units system
  Units::Init();

  // Init built-in plug-ins
#ifdef GYOTO_BUILTIN_STDPLUG
  __GyotostdplugInit();
#endif
#ifdef GYOTO_BUILTIN_LORENEPLUG
  __GyotoloreneInit();
#endif

  // Load DL plug-ins

  if (!cpluglist) cpluglist = getenv("GYOTO_PLUGINS");
  if (!cpluglist) cpluglist = GYOTO_DEFAULT_PLUGINS;

  std::string pluglist = cpluglist;

  if (pluglist.length()) {
    size_t first=0, last=0;
    string curplug="";
    int nofail=0;
    while (pluglist.length()) {
      last=pluglist.find(",");
      nofail=0;
      curplug=pluglist.substr(0, last);
      if (debug())
	cerr << "DEBUG: first: " << first << ", last: " << last
	     << ", pluglist: |" << pluglist << "|"
	     << ", curplug: |" << curplug << "|" << endl;
      if (last <= pluglist.length()) pluglist=pluglist.substr(last+1);
      else pluglist="";
      if (!curplug.compare(0, 7, "nofail:")) {
	curplug = curplug.substr(7);
	nofail=1;
      } else if (!curplug.compare(0, 7, "nowarn:")) {
	curplug = curplug.substr(7);
	nofail=2;
      }

      Gyoto::loadPlugin(curplug.c_str(), nofail);
      nofail=0;
    }
  }

  if (debug()) Register::list();

}

Register::Entry::Entry(std::string name,
		       Gyoto::SmartPointee::Subcontractor_t* subcontractor,
		       Register::Entry* next)
  : name_(name), subcontractor_(subcontractor), next_(next), plugin_(GyotoRegisterCurrentPlugin)
{}

Register::Entry::~Entry() { if (next_) delete next_; }


#ifndef GYOTO_NO_DEPRECATED
#warning Embedding deprecated method.\
  Define GYOTO_NO_DEPRECATED to disable.
Gyoto::SmartPointee::Subcontractor_t*
Register::Entry::getSubcontractor(std::string name, int errmode) {
  std::string plugin("");
  return getSubcontractor(name, plugin, errmode);
}
#endif

Gyoto::SmartPointee::Subcontractor_t*
Register::Entry::getSubcontractor(std::string name, std::string &plugin, int errmode) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(name);
    GYOTO_DEBUG_EXPR(name_);
    GYOTO_DEBUG_EXPR(plugin);
    GYOTO_DEBUG_EXPR(plugin_);
    GYOTO_DEBUG_EXPR(errmode);
  GYOTO_ENDIF_DEBUG
# endif
  bool any_plugin = (plugin == "");
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(any_plugin);
  GYOTO_ENDIF_DEBUG
# endif
  if (name_==name && (any_plugin || (plugin_ == plugin))) {
    if (any_plugin) plugin=plugin_;
    return subcontractor_;
  }
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(next_);
  GYOTO_ENDIF_DEBUG
# endif
  if (next_) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
      GYOTO_DEBUG << "recursing\n" ;
    GYOTO_ENDIF_DEBUG
#   endif
    return next_ -> getSubcontractor(name, plugin, errmode);
  }
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG << name << " not found in plug-in \"" << plugin << "\"\n" ;
  GYOTO_ENDIF_DEBUG
# endif
  if (errmode) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
      GYOTO_DEBUG << "returning\n" ;
    GYOTO_ENDIF_DEBUG
#   endif
    return NULL;
  }
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG << "throwing error\n" ;
  GYOTO_ENDIF_DEBUG
# endif
  GYOTO_ERROR ("Unregistered kind: "+name);
  return NULL; // will never get there, avoid compilation warning
}

std::string Register::Entry::name() {return name_;}
std::string Register::Entry::plugin() {return plugin_;}
Register::Entry* Register::Entry::next() {return next_;}

void Gyoto::Register::list() {
  Register::Entry* entry = NULL;

  cout <<
"Gyoto will look for plug-ins first in the run-time linker default locations\n"
"(typically includes directories listed in e.g. $LD_LIBRARY_PATH), then in the\n"
"following locations:" << endl;

  for (const auto &path : GyotoPluginPath)
    cout << path << endl;
  cout << endl;

  cout << "List of loaded plug-ins:" << endl;
  for (size_t i=0; i < GyotoRegisteredPlugins.size(); ++i)
    cout << "    " << GyotoRegisteredPlugins[i] << endl;

  cout << "List of available Metrics:" << endl;
  for (entry = Metric::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << " (in plug-in: " << entry -> plugin_ << ")" << endl;
  
  cout << "List of available Astrobjs:" << endl;
  for (entry = Astrobj::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << " (in plug-in: " << entry -> plugin_ << ")" << endl;
  
  cout << "List of available Spectra:" << endl;
  for (entry = Spectrum::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << " (in plug-in: " << entry -> plugin_ << ")" << endl;
    
  
  cout << "List of available Spectrometers:" << endl;
  for (entry = Spectrometer::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << " (in plug-in: " << entry -> plugin_ << ")" << endl;
    
}
