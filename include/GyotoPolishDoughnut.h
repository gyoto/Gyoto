/**
 * \file GyotoPolishDoughnut.h
 * \brief A toroïdal accretion structure
 *
 *  Reference: Straub, O.; Vincent, F. H.; Abramowicz, M. A.;
 *  Gourgoulhon, E.; &amp; Paumard, T. 2012, <STRONG>Modelling the
 *  black hole silhouette in Sagittarius A* with ion tori</STRONG>,
 *  A&amp;A 543:83.
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
 * \brief A toroïdal accretion structure
 *
 *  Reference: Straub, O.; Vincent, F. H.; Abramowicz, M. A.;
 *  Gourgoulhon, E.; &amp; Paumard, T. 2012, <STRONG>Modelling the
 *  black hole silhouette in Sagittarius A* with ion tori</STRONG>,
 *  A&amp;A 543:83.
 * 
 */
class Gyoto::Astrobj::PolishDoughnut
: public Astrobj::Standard,
  protected Gyoto::Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PolishDoughnut>;

 // Data : 
 // -----
protected:
  /**
   * \brief KerrBL metric
   *
   * This Astrobj::Generic subclass only works (so far?) in KerrBL
   * metric.  We store a copy of the Astrobj::Generic::gg_
   * SmartPointer appropriately cast. The two pointers point to the
   * same instance.
   */
 SmartPointer<Gyoto::Metric::KerrBL> gg_; 
 double l0_; ///< Angular momentum. Tied to PolishDoughnut::lambda_.
 double lambda_; ///< Adimentionned angular momentum
 double W_surface_; ///< Potential surface value. Tied to PolishDoughnut::lambda_.
 double W_centre_; ///< Potential central value. Tied to PolishDoughnut::lambda_.
 double r_cusp_; ///< Cusp radius in geometrical units. Tied to PolishDoughnut::lambda_.
 double r_centre_; ///< Central radius in geometrical units. Tied to PolishDoughnut::lambda_.
 double DeltaWm1_; ///< 1./(W_centre_ - W_surface_);
 double central_density_; ///< Central density in kg/L (same as g cm^-3)
 /*
   WARNING! so far (jan. 2014) central_density_ is density_central
   in standard Polish doughnut model, but it is
   density_central*c2+pressure_central in Komissarov model
  */
 double centraltemp_over_virial_; ///< T<SUB>center</SUB>/T<SUB>virial</SUB>
 double beta_; ///< P<SUB>magn</SUB>/P<SUB>gas</SUB>
 double aa_; ///< PolishDoughnut::gg_ spin, cached when setting PolishDoughnut::lambda_
 double aa2_; ///< aa_<SUP>2</SUP>
 size_t spectral_oversampling_;///< Oversampling used in integrateEmission()
 bool komissarov_; ///< 1 if Komissarov model is integrated
 bool angle_averaged_; ///< 1 if Komissarov model should be angle averaged

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
 double getL0() const; ///< Get PolishDoughnut::l0_
 // void   setL0(double l0); set by lambda_

 double getLambda() const; ///< Get PolishDoughnut::lambda_
 void   setLambda(double lambda); ///< Set PolishDoughnut::lambda_

 double getCentralDensity() const; ///< Get PolishDoughnut::central_density_
 double getCentralDensity(std::string unit) const; ///< Get PolishDoughnut::central_density_ in specified unit
 void   setCentralDensity(double density); ///< Set PolishDoughnut::central_density_
 void   setCentralDensity(double density, std::string unit); ///< Set PolishDoughnut::central_density_ in specified unit

 double getCentralTempOverVirial() const; ///< Get PolishDoughnut::centraltemp_over_virial_
 void   setCentralTempOverVirial(double val); ///< Set PolishDoughnut::centraltemp_over_virial_

 double getBeta() const; ///< Get PolishDoughnut::beta_
 void   setBeta(double beta);///< Set PolishDoughnut::beta_

 void   setSpectralOversampling(size_t); ///< Set PolishDoughnut::spectral_oversampling_
 size_t getSpectralOversampling() const ; ///< Get PolishDoughnut::spectral_oversampling_

 bool komissarov() const; ///< Get PolishDoughnut::komissarov_
 void komissarov(bool komis); ///< Set PolishDoughnut::komissarov_

 // Read only members, depend on lambda
 double getWsurface() const; ///< Get PolishDoughnut::W_surface_
 double getWcentre() const; ///< Get PolishDoughnut::W_centre_
 double getRcusp() const; ///< Get PolishDoughnut::r_cusp_
 double getRcentre() const; ///< Get PolishDoughnut::r_centre_

 virtual Gyoto::SmartPointer<Gyoto::Metric::Generic> metric() const;
 virtual void metric(Gyoto::SmartPointer<Gyoto::Metric::Generic>);

 // ASTROBJ API
 // -----------
 int Impact(Photon *ph, size_t index,
			    Astrobj::Properties *data);

 virtual double operator()(double const coord[4]) ;

 virtual int setParameter(std::string name,
			  std::string content,
			  std::string unit) ;
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
#endif

  // ASTROBJ processHitQuantities API
  // --------------------------------
protected:
  /**
   * \brief Update PolishDoughnut::aa_
   *
   * See Hook::Listener::tell().
   */
  virtual void tell(Gyoto::Hook::Teller * msg);
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  using Standard::integrateEmission;

  /**
   * \brief &int;<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP> I<SUB>&nu;</SUB> d&nu; (or j<SUB>&nu;</SUB>)
   *
   * PolishDought::emission() is slow. Integrating it is very
   * slow. This implementation simply oversamples the spectrum by
   * PolishDoughnut::spectral_oversampling_ and sums the trapezoids.
   *
   * For general documentation, see Astrobj::Generic::integrateEmission(double * I, double const * boundaries, size_t const * chaninds, size_t nbnu, double dsem, double *cph, double *co) const .
   */
  virtual void integrateEmission(double * I, double * boundaries,
				 size_t* chaninds, size_t nbnu,
				 double dsem, double *cph, double *co) const;

  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]) const;
  virtual void emission(double Inu[], double nu_em[], size_t nbnu,
			double dsem, double coord_ph[8],
			double coord_obj[8]=NULL) const ;

  void emission_komissarov(double Inu[], double nu_em[], size_t nbnu,
			double dsem, double coord_ph[8],
			double coord_obj[8]=NULL) const ;

  double emissionBrems(double nu_em, double nu_crit, 
		       double numax, double T_electron,
		       double n_e, double n_j,
		       double amplification,
		       double Cbrems,
		       int comptonorder) const;
  ///< Bremsstrahlung proxy for emission()
  double emissionSynch(double nu_em, double nu_crit,
		       double numax, double nu_0,
		       double T_electron,
		       double amplification,
		       double Csynch, 
		       double alpha1, double alpha2,
		       double alpha3, double preff,
		       int comptonorder) const;

  double emissionSynchro_komissarov_direction(double Theta_elec, 
					      double number_density,
					      double nuem,
					      double nuc, 
					      double theta
					      ) const;
  double emissionSynchro_komissarov_averaged(double Theta_elec, 
					     double number_density,
					     double nuem,
					     double nuc 
					     ) const;
  ///< Synchrotron proxy for emission()
  double transmission(double nuem, double dsem, double coord_ph[8]) const ;
  double BBapprox(double nuem, double Te) const; ///< Approximated Black-Body function
  static double funcxM(double alpha1, double alpha2, double alpha3, double xM);
  ///< Mahadevan 96 fit function
 // PURELY INTERNAL TO ASTROBJ
 // --------------------------
 private:
  double potential(double r, double theta) const;
  ///< Potential defining shape, used by operator()()

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
    intersection_t(PolishDoughnut* parent);
    PolishDoughnut * papa;
    virtual double operator() (double) const;
  };
  friend class intersection_t;

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
   /**
    * \brief Parameter array
    *
    * \code
    *   double       rr = par[0] ;
    *   double      n_e = par[1] ;
    *   double       BB = par[2] ;
    *   double       Te = par[3] ;
    *   double   alpha1 = par[4] ;
    *   double   alpha2 = par[5] ;
    *   double   alpha3 = par[6] ;
    * \endcode
    */
    double const * par;
    virtual double operator() (double) const;
  };
 intersection_t intersection; ///< double intersection(double) Functor

 public:
 static double bessi0(double xx);///< Modified Bessel function I<SUB>0</SUB>
 static double bessi1(double xx);///< Modified Bessel function I<SUB>1</SUB>
 static double bessk0(double xx);///< Modified Bessel function K<SUB>0</SUB>
 static double bessk1(double xx);///< Modified Bessel function K<SUB>1</SUB>
 static double bessk(int nn, double xx);///< Modified Bessel function

 // Outputs
 // -------
 public:

 /// Display
 friend std::ostream& operator<<(std::ostream& , const PolishDoughnut& ) ;

 public:
 

};

#endif
