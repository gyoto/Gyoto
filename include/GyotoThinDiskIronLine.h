#ifndef __GyotoThinDiskIronLine_h
#define __GyotoThinDiskIronLine_h

#include <GyotoThinDisk.h>

namespace Gyoto {
  namespace Astrobj {
    class ThinDiskIronLine;
  };
};

class Gyoto::Astrobj::ThinDiskIronLine
: public Gyoto::Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDiskIronLine>;
 private:
  double plindex_; ///< power law index for line emission
  double linefreq_; ///< intrinsic line frequency (Hz)
  double cutradius_; ///< r<cutradius_ -> emission = 0
 public:
  ThinDiskIronLine();
  ThinDiskIronLine(const ThinDiskIronLine &o);
  virtual ~ThinDiskIronLine();
  virtual ThinDiskIronLine * clone() const ;
  virtual double emission(double nu_em, double dsem,
			  double c_ph[8], double c_obj[8]) const;
  void getVelocity(double const pos[4], double vel[4]);

  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit);
#ifdef GYOTO_USE_XERCES
  //virtual void fillElement(FactoryMessenger *fmp) const ;
#endif
#endif
};
