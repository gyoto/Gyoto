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
			  double c_ph[8], double c_obj[8]) const;
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
