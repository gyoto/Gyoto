#include "GyotoNull.h"

using namespace Gyoto;

extern "C" void __GyotoPluginInit() {
  Astrobj::Register("Null", &Astrobj::Subcontractor<Astrobj::Null>);
}
