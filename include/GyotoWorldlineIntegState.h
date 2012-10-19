/**
 * \file GyotoWorldlineIntegState.h
 * \brief Geodesic integration state
 * 
 *
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoWorldlineIntegState_H_
#define __GyotoWorldlineIntegState_H_ 

#include <iostream>
#include <fstream>

namespace Gyoto {
  class WorldlineIntegState;
}

#include <GyotoUtils.h>
#include <GyotoSmartPointer.h>
#include <GyotoAstrobj.h>
#include <GyotoMetric.h>
#include <GyotoWorldline.h>

/**
 * \class Gyoto::WorldlineIntegState
 * \brief Current state of a geodesic integration
 */

class Gyoto::WorldlineIntegState : SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::WorldlineIntegState>;

 private:
  Gyoto::SmartPointer<Gyoto::Metric::Generic> gg_;
  
 protected:
  double coord_[8];
  double coordnew_[8];
  double norm_;
  double normref_;
  double delta_;
  double h1_;
  double deltainit_;

 public:
  WorldlineIntegState();
  WorldlineIntegState(Gyoto::SmartPointer<Gyoto::Metric::Generic> gg,
		      const double *coord, const double delta);
  
  /**
   * \param coord[8] on input: old position-velocity, on output: new
   *                 position-velocity;
   * \param delta    integrating step in proper time. Set to 0 to use
   *                 default step, presumably adaptive.
   */
  virtual int nextStep(Gyoto::Worldline* line, double *coord, double delta=0.);
  double get_delta();
  void set_delta(double delta);
  void setCoord(double coord[8]);
  virtual ~WorldlineIntegState();
};

#endif
