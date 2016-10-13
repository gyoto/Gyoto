#include "GyotoNull.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

Null::Null() : Generic("Null") {}
Null::Null(const Null& o) : Generic(o){}
Null * Null::clone() const { return new Null(*this); }
Null::~Null() {}
int Null::Impact(Photon *, size_t , Astrobj::Properties *) {return 0;}
