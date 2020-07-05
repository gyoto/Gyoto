/**
 * \file GyotoShift.h
 * \brief Shift a metric
 * 
 */

/*
    Copyright 2020 Thibaut Paumard & Frédéric Vincent

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


#ifndef __GyotoShift_H_
#define __GyotoShift_H_

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric { class Shift; }
}

/**
 * \class Gyoto::Metric::Shift
 * \brief The Shift flat-space metric
 * 
 * Use &lt;Cartesian&gt; or &lt;/Spherical&gt; to select the coordinate system
 * kind.
 */

class Gyoto::Metric::Shift
  : public Gyoto::Metric::Generic,
    public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::Shift>;

 protected:
  Gyoto::SmartPointer<Gyoto::Metric::Generic> submet_;
  double offset_[4];
  
 public:
  // This is the bare minimum of what a Metric class must implement:
  GYOTO_OBJECT;
  Shift();
  virtual ~Shift();
  virtual Shift* clone() const ;

  virtual SmartPointer<Metric::Generic> subMetric() const;
  virtual void subMetric(SmartPointer<Metric::Generic>) ;
  virtual std::vector<double> offset() const; ///< Get vector copy of #pos_
  virtual void offset(std::vector<double> const&); ///< Set #pos_ from vector

  virtual void mass(const double);        ///< Set mass used in unitLength()
  using Generic::gmunu;
  virtual void gmunu(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;
  using Generic::gmunu_up;
  virtual void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;
  virtual void jacobian(double ARGOUT_ARRAY3[4][4][4], const double IN_ARRAY1[4]) const;
  virtual int isStopCondition(double const coord[8]) const;

# ifdef GYOTO_USE_XERCES
  virtual void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const ;
  virtual void setParameters(FactoryMessenger *fmp) ;
# endif

  virtual void tell(Gyoto::Hook::Teller *msg);
  
};

#endif
