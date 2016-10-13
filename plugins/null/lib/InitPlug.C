// #include all headers necessary to register subcontractors below
#include "GyotoNull.h"

using namespace Gyoto;

// Rename replace "null" with actual name of plug-in in the init
// function name below:
extern "C" void __GyotonullInit() {
  // Register subcontractors for all new classes
  Astrobj::Register("Null", &Astrobj::Subcontractor<Astrobj::Null>);
}
