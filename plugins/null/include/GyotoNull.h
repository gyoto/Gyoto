/**
 * \file GyotoNull.h
 * \brief Null Astrobj, just for investigating geodesics
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
 *   Copyright (c) year  your_name
 *
 *
 */


#ifndef __GyotoNullAstrobj_H_ 
#define __GyotoNullAstrobj_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class Null; }
  class FactoryMessenger;
  namespace Spectrum {
    class Generic;
  }
}

#include <GyotoAstrobj.h>
#include <GyotoMetric.h>

/**
 * \class Gyoto::Null. 
 * \brief Empty Astrobj
 *
 *  The target of ray-traced Gyoto::Photon
 */
class Gyoto::Astrobj::Null : public Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Null>;
 public:
  Null();
  Null(const Null& o);
  virtual Null * clone() const ;
  virtual ~Null() ;
  virtual int Impact(Photon *ph, size_t index, Astrobj::Properties *data=NULL);
};

#endif
