/*
    Copyright 2013-2015, 2018 Frederic Vincent, Thibaut Paumard

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

/** \file GyotoThinDiskIronLine.h */

#ifndef __GyotoThinDiskIronLine_h
#define __GyotoThinDiskIronLine_h

#include <GyotoThinDisk.h>

namespace Gyoto {
  namespace Astrobj {
    class ThinDiskIronLine;
  }
}

class Gyoto::Astrobj::ThinDiskIronLine
: public Gyoto::Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDiskIronLine>;
 private:
  double plindex_; ///< power law index for line emission
  double linefreq_; ///< intrinsic line frequency (Hz)
  double cutradius_; ///< r<cutradius_ -> emission = 0
 public:
  GYOTO_OBJECT;
  ThinDiskIronLine();
  ThinDiskIronLine(const ThinDiskIronLine &o);
  virtual ~ThinDiskIronLine();
  virtual ThinDiskIronLine * clone() const ;

  using ThinDisk::emission;
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;
  void getVelocity(double const pos[4], double vel[4]);


  // standard pairs of accessors
  GYOTO_OBJECT_ACCESSORS(double, PowerLawIndex);
  GYOTO_OBJECT_ACCESSORS(double, LineFreq);
  void LineFreq(double v, std::string const &u);
  double LineFreq(std::string const &u)const;
  GYOTO_OBJECT_ACCESSORS(double, CutRadius);
  void CutRadius(double v, std::string const &u);
  double CutRadius(std::string const &u)const;

};
#endif
