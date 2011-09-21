/**
 * \file GyotoThinInfiniteDiskKS.h
 * \brief A geometrically thin, optically thick disk
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2011 Frederic Vincent

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

#ifndef __GyotoThinInfiniteDiskKS_H_ 
#define __GyotoThinInfiniteDiskKS_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class ThinInfiniteDiskKS; }
}

//#include <GyotoMetric.h>
#include <GyotoKerrKS.h>
#include <GyotoAstrobj.h>

/**
 * \class Gyoto::Astrobj::ThinInfiniteDiskKS
 * \brief Geometrically thin disk in Metric::KerrKS metric
 * 
 *   This class describes a disk contained in the z=0 (equatorial) plane, extending
 *   from r=r_ISCO to r=infinity.  The flux emitted at radius r is given. 
 *   The metric is supposed to be KerrKS.
 * 
 */
  class Gyoto::Astrobj::ThinInfiniteDiskKS : public Astrobj::Generic {
    friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinInfiniteDiskKS>;

  /*
   */

 // Data : 
  // -----
 protected:
  
  SmartPointer<Gyoto::Metric::KerrKS> gg_ ; 

  double Lr_;//luminosity at the position where the disk has been hit (beaming taken into account)
  double rmin_;//ISCO radius
    
  // Constructors - Destructor
  // -------------------------
 public:
  
  ThinInfiniteDiskKS(const SmartPointer<Metric::KerrKS>& metric);
                   ///< Standard constructor

  ThinInfiniteDiskKS(const ThinInfiniteDiskKS& ) ;  ///< Copy constructor
  virtual ThinInfiniteDiskKS* clone () const; ///< Cloner
 
  /// Constructor from a file (see \c sauve(FILE*) )
  //ThinInfiniteDiskKS(FILE *) ;                    
  
  virtual ~ThinInfiniteDiskKS() ;                        ///< Destructor
  
  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another ThinInfiniteDiskKS
  void operator=(const ThinInfiniteDiskKS&) ;        
  
  // Accessors
  // ---------
 public:
  virtual int Impact(Photon *ph, size_t index, Astrobj::Properties *data=NULL);
  virtual double emission(double nu_em, double dsem,
			  double c_ph[8], double c_obj[8]) const;
  // Outputs
  // -------
 public:
  //virtual void sauve(FILE *) const ;            ///< Save in a file
  
  /// Display
  friend std::ostream& operator<<(std::ostream& , const ThinInfiniteDiskKS& ) ;
  
 
 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
  ///< called from Factory
  static Astrobj::Subcontractor_t Subcontractor;
  static void Init();
#endif
 
};

#endif
