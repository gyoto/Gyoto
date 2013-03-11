/**
* \file GyotoPolishDoughnut.h
* \brief A toroïdal accretion structure
*
*  Beware: the PolishDoughnut now in SI
*
*/

/*
    Copyright (c) 2012 Frederic Vincent, Odele Straub, Thibaut Paumard

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

#ifndef __GyotoPolishDoughnut_H_ 
#define __GyotoPolishDoughnut_H_ 

namespace Gyoto{
  namespace Astrobj { class PolishDoughnut; }
 class FactoryMessenger;
}

//#include <GyotoMetric.h>
#include <GyotoKerrBL.h>
#include <GyotoStandardAstrobj.h>
#include <GyotoFunctors.h>
#include <GyotoHooks.h>
//#include <GyotoPolishDoughnutCst.h>

/**
* \class Gyoto::PolishDoughnut
* \brief 
* 
 */
class Gyoto::Astrobj::PolishDoughnut
: public Astrobj::Standard,
  protected Gyoto::Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PolishDoughnut>;

 // Data : 
 // -----
private:
 SmartPointer<Gyoto::Metric::KerrBL> gg_;
 double l0_; ///< torus angular momentum. Tied to lambda.
 double lambda_; ///< torus adimentionned angular momentum
 double W_surface_; ///< torus potential surface value. Tied to lambda.
 double W_centre_; ///< torus potential central value. Tied to lambda.
 double r_cusp_; ///< torus cusp radius in geometrical units. Tied to lambda.
 double r_centre_; ///< torus central radius in geometrical units. Tied to lambda.
 double DeltaWm1_; ///< 1./(W_centre_ - W_surface_);
 double temperature_ratio_; ///< central ion/elec temperature ratio
 double central_density_; ///< central density in kg/L (same as g cm^-3)
 double centraltemp_over_virial_; ///< T_center/Tvirial
 double beta_; ///< Pmagn/Pgas
 int use_specific_impact_ ;///< PolishDoughnut::Impact_() or Standard::Impact()
 double aa_; ///< Metric spin, cached when setting lambda_
 double aa2_; ///< aa_^2
 size_t spectral_oversampling_;///< oversampling used in integrateEmission()

 // Constructors - Destructor
 // -------------------------
public:
 PolishDoughnut() ; ///< Default constructor
 //PolishDoughnut(SmartPointer<Metric::KerrBL> gg, double lambda) ; ///< Standard constructor
 PolishDoughnut(const PolishDoughnut& ) ;                ///< Copy constructor
 virtual  PolishDoughnut* clone() const;
 virtual ~PolishDoughnut() ;                        ///< Destructor


 // Mutators / assignment
 // ---------------------
public:
 /// Assignment to another PolishDoughnut
 void operator=(const PolishDoughnut&) ;        

 // Accessors
 // ---------
public:
 double getL0() const;
 // void   setL0(double l0); set by lambda_

 double getLambda() const;
 void   setLambda(double lambda);

 double getTemperatureRatio() const;
 void   setTemperatureRatio(double temp);

 double getCentralDensity() const;
 double getCentralDensity(std::string unit) const;
 void   setCentralDensity(double density);
 void   setCentralDensity(double density, std::string unit);

 double getCentralTempOverVirial() const;
 void   setCentralTempOverVirial(double val);

 double getBeta() const;
 void   setBeta(double beta);

 void   setSpectralOversampling(size_t); ///< set spectral_oversampling_
 size_t getSpectralOversampling() const ; ///< get spectral_oversampling_

 // Read only members, depend on lambda
 double getWsurface() const;
 double getWcentre() const;
 double getRcusp() const;
 double getRcentre() const;


 virtual Gyoto::SmartPointer<Gyoto::Metric::Generic> getMetric() const;
 virtual void setMetric(Gyoto::SmartPointer<Gyoto::Metric::Generic>);
 void useSpecificImpact(int yes=1);

 // ASTROBJ API
 // -----------
 /**
  * Depending on use_specific_impact_. See useSpecificImpact().
  */
 int Impact(Photon *ph, size_t index,
			    Astrobj::Properties *data);
 ///< Call either PolishDoughnut::Impact() or Standard::Impact()

 /**
  * This should not be better than Standard::Impact(). Call
  * useSpecificImpact() or Set &lt;UseSpecifictImpact/&gt; in the XML
  * file to use it.
  */
 int Impact_(Photon *ph, size_t index,
			    Astrobj::Properties *data);
 ///< A specific implementation of Generic::Impact()
 virtual double operator()(double const coord[4]) ;
 ///< Called by Astrobj::Generic::Impact()

 virtual int setParameter(std::string name,
			  std::string content,
			  std::string unit) ;
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
#endif

  // ASTROBJ processHitQuantities API
  // --------------------------------
protected:
  virtual void tell(Gyoto::Hook::Teller * msg);
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  using Standard::integrateEmission;
  virtual void integrateEmission(double * I, double * boundaries,
				 size_t* chaninds, size_t nbnu,
				 double dsem, double *cph, double *co) const;
  virtual void emission(double Inu[], double nu_em[], size_t nbnu,
			  double dsem, double coord_ph[8],
			  double coord_obj[8]) const;
  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]) const;
  double emissionBrems(double nu_em, double nu_crit, 
		       double numax, double T_electron,
		       double n_e, double n_j,
		       double amplification,
		       double Cbrems,
		       int comptonorder) const;
  double emissionSynch(double nu_em, double nu_crit,
		       double numax, double nu_0,
		       double T_electron,
		       double amplification,
		       double Csynch, 
		       double alpha1, double alpha2,
		       double alpha3, double preff,
		       int comptonorder) const;
  double transmission(double nuem, double dsem, double coord_ph[8]) const ;
  double BBapprox(double nuem, double Te) const;
  static double funcxM(double alpha1, double alpha2, double alpha3, double xM);

 // PURELY INTERNAL TO ASTROBJ
 // --------------------------
private:
 double potential(double r, double theta) const;

 /**
  * \class Gyoto::Astrobj::PolishDoughnut::intersection_t
  * \brief double intersection(double) Functor class
  *
  * Implement former double intersection(double) function as a
  * Gyoto::Functor::Double_Double_const subclass to access generic
  * root-finding methods.
  *
  * This class is instanciated in a single
  * PolishDoughnut::intersection member.
  */
  class intersection_t : public Gyoto::Functor::Double_Double_const {
  public:
    double aa_;
    double aa2_;
    double l0_;
    virtual double operator() (double) const;
  };
 /**
  * \class Gyoto::Astrobj::PolishDoughnut::transcendental_t
  * \brief double transcendental(double) Functor class
  *
  * Implement former double transcendental(double, double*) function
  * as a Gyoto::Functor::Double_Double_const subclass to access
  * generic root-finding methods.
  *
  * This class is as a local variable in PolishDoughnut::emission()
  */
  class transcendental_t : public Gyoto::Functor::Double_Double_const {
  public:
    double const * par;
    virtual double operator() (double) const;
  };
 intersection_t intersection; ///< double intersection(double) Functor

 public:
 static double bessi0(double xx);
 static double bessi1(double xx);
 static double bessk0(double xx);
 static double bessk1(double xx);
 static double bessk(int nn, double xx);

 // Outputs
 // -------
 public:

 /// Display
 friend std::ostream& operator<<(std::ostream& , const PolishDoughnut& ) ;

 public:
 

};

#endif
