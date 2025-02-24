/*
    Copyright 2025 Filipe Costa, Frédéric Vincent, Thibaut Paumard

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

#ifndef __BalasinGrumillerheader_H_
#define __BalasinGrumillerheader_H_

#include <iostream> // for basic I/O
#include <cmath>    // for mathematical functions
#include <fstream>
#include <string>
#include <vector>

#include <GyotoSmartPointer.h>
#include <GyotoObject.h>
#include <GyotoAstrobj.h>
#include <GyotoHooks.h>
#include <GyotoDefs.h>
#include "GyotoWorldline.h"
#include <GyotoMetric.h> // Include Gyoto Metric base class

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif


namespace Gyoto {
  namespace Metric { class BalasinGrumiller; }
}

/**
 * \class Gyoto::Metric::BalasinGrumiller
 * \brief Balasin-Grumiller solution in spherical coordinates, under the suggested approximation
 * \nu=0.
 *
 * This metric was originally proposed by in Int. J. Mod. Phys. D 17 (2008) 475 [astro-ph/0602519]
 * as a galactic model attempting to replace dark matter with relativistic effects in 
 * explaining galactic rotation curves.
 *
 * It has subsequently been shown in Costa et al. Phys. Rev. D 108, 044056 (2023) [arXiv:2303.17516] 
 * (Sec. IV.E) to be unsuitable for describing galaxies, consisting of a combination of two pathologies:
 * - Unphysical singularities
 * - An unsuitable choice of reference observers
 *  (the model being actually static with respect to asymptotic inertial frames)
 *
 * Further investigation in Costa-Natário Phys. Rev. D 110, 064056 (2024) [arXiv:2312.12302] showed 
 * that the metric produces also lensing effects starkly at odds with observations.
 * This implementation allows users to verify this conclusion independently. The metric parameters 
 * V0, R, r0 can be specified through the associated XML file. 
 */

class Gyoto::Metric::BalasinGrumiller : public Gyoto::Metric::Generic 
{
friend class Gyoto::SmartPointer<Gyoto::Metric::BalasinGrumiller>;
protected:
 double V0value_ ; 
 double Rvalue_ ;
 double r0value_ ;
public:
GYOTO_OBJECT;
BalasinGrumiller(); 
virtual BalasinGrumiller* clone() const ;

 // Mutators / assignment
  // ---------------------
 public:
  // default operator= is fine
  void V0value(const double V0value); ///< Set V0
  void Rvalue(const double Rvalue); ///< Set R
  void r0value(const double r0value); ///< Set r0

  // Accessors
  // ---------
 public:
  double V0value() const ; ///< Returns V0
  double Rvalue() const ; ///< Returns R
  double r0value() const ; ///< Returns r0

  void gmunu(double g[4][4], const double x[4]) const ;
  int christoffel(double dst[4][4][4], const double x[4]) const ;

  // Those two are implemented as examples.
 // double gmunu(const double x[4], int mu, int nu) const ;
 // double christoffel(const double coord[4],
//  const int alpha, const int mu, const int nu) const ;
 };

#endif
