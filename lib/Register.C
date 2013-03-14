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
#include "GyotoRegister.h"
#include "GyotoUtils.h"
#include "GyotoAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoSpectrum.h"
#include "GyotoConverters.h"

#include <dlfcn.h>
#include <cstdlib>
#include <string>
#include <iostream>

using namespace Gyoto;
using namespace std;

typedef void GyotoInitFcn();

void Gyoto::loadPlugin(char const*const name, int nofail) {
  string dlfile = "libgyoto-" ;
  dlfile += name ;
  dlfile += "." ;
  dlfile += GYOTO_PLUGIN_SFX ;
  string dlfunc = "__Gyoto";
  dlfunc += name ;
  dlfunc += "Init";
  void* handle = NULL;
  GyotoInitFcn* initfcn = NULL;
  char * err = NULL;

  if (debug()) cerr << "DEBUG: loading plug-in: " << name
		    << " from file: " << dlfile << endl;
  handle = dlopen(dlfile.c_str(), RTLD_LAZY | RTLD_GLOBAL);
  if (!handle) {
    string dlpath = GYOTO_PREFIX ;
    dlpath += "/lib/gyoto/" ;
    string dlfull = dlpath + dlfile;
    handle = dlopen(dlfull.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    if (!handle) {
      dlfull = dlpath ;
      dlfull += GYOTO_SOVERS ;
      dlfull += "/" ;
      dlfull += dlfile ;
      handle = dlopen(dlfull.c_str(), RTLD_LAZY | RTLD_GLOBAL);
      if (!handle && nofail) {
	if (verbose() >= GYOTO_DEFAULT_VERBOSITY)
	  cerr << "WARNING: unable to load optional plug-in " << dlfile << endl;
	return;
      }
    }
  }
  if ( (err=dlerror()) ) throwError(err);
  if (!handle) throwError((string("Failed to load plug-in ")+dlfile).c_str());
  if (debug()) cerr << "DEBUG: calling plug-in init function " << dlfunc << endl;
  initfcn = (GyotoInitFcn*)dlsym(handle, dlfunc.c_str());
  if ( (err=dlerror()) ) throwError(err);
  (*initfcn)();
}

#if defined GYOTO_USE_XERCES
void Gyoto::Register::init(char const *  cpluglist) {

  // Clean registers
  Metric::initRegister();
  Astrobj::initRegister();
  Spectrum::initRegister();
  // This cleans and fills Spectometer::Register_
  Spectrometer::initRegister();

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
  : name_(name), subcontractor_(subcontractor), next_(next)
{}

Register::Entry::~Entry() { if (next_) delete next_; }


Gyoto::SmartPointee::Subcontractor_t*
Register::Entry::getSubcontractor(std::string name, int errmode) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(name);
    GYOTO_DEBUG_EXPR(errmode);
  GYOTO_ENDIF_DEBUG
# endif
  if (name_==name) return subcontractor_;
  if (next_) return next_ -> getSubcontractor(name, errmode);
  if (errmode) return NULL;
  throwError ("Unregistered kind: "+name);
  return NULL; // will never get there, avoid compilation warning
}

void Gyoto::Register::list() {
  Register::Entry* entry = NULL;

  cout << "List of available Metrics:" << endl;
  for (entry = Metric::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << endl;
  
  cout << "List of available Astrobjs:" << endl;
  for (entry = Astrobj::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << endl;
  
  cout << "List of available Spectra:" << endl;
  for (entry = Spectrum::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << endl;
    
  
  cout << "List of available Spectrometers:" << endl;
  for (entry = Spectrometer::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << endl;
    
}
#endif
