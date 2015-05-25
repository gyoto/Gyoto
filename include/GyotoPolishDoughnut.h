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
    Copyright (c) 2012-2015 Frederic Vincent, Odele Straub, Thibaut Paumard

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

#include <GyotoStandardAstrobj.h>
#include <GyotoFunctors.h>
#include <GyotoHooks.h>
#include <GyotoBlackBodySpectrum.h>

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
 SmartPointer<Spectrum::BlackBody> spectrumBB_;
 double l0_; ///< Angular momentum. Tied to PolishDoughnut::lambda_.
 double lambda_; ///< Adimentionned angular momentum
 double W_surface_; ///< Potential surface value. Tied to PolishDoughnut::lambda_.
 double W_centre_; ///< Potential central value. Tied to PolishDoughnut::lambda_.
 double r_cusp_; ///< Cusp radius in geometrical units. Tied to PolishDoughnut::lambda_.
 double r_centre_; ///< Central radius in geometrical units. Tied to PolishDoughnut::lambda_.
 double r_torusouter_ ; ///< Torus outer coordinate radius. Tied to PolishDoughnut::lambda_.
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

 /// fraction of thermal energy in non-thermal electrons
 /**
  * Obsiously, 0 means no non-thermal electrons at all. In this case,
  * the (trivial) non-thermal computations are skipped. Ther is thus
  * non need for an additional "nonthermal_" flag.
  */
 double deltaPL_;
 double expoPL_; ///< exponent of the non-thermal powerlaw = -expoPL_

 bool adaf_; ///< true to switch to an ADAF model rather tha Polish doughnut
 double ADAFtemperature_; ///< ADAF central temperature
 double ADAFdensity_; ///< ADAF central density

 bool changecusp_; ///< true to apply the fishy rcusp_ change (to be changed)
 bool rochelobefilling_; ///< true if torus filling its Roche lobe
 bool defangmomrinner_; ///< true if torus defined from l0 and rin
 double rintorus_; ///< Inner radius of the doughnut

 // Constructors - Destructor
 // -------------------------
public:
 GYOTO_OBJECT; // This object has Properties
#ifdef GYOTO_USE_XERCES
 // We need to filter some properties when writing XML
 void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;
#endif
 PolishDoughnut() ; ///< Default constructor
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

 double lambda() const; ///< Get PolishDoughnut::lambda_
 void   lambda(double lambda); ///< Set PolishDoughnut::lambda_

 double centralDensity() const; ///< Get PolishDoughnut::central_density_
 double centralDensity(std::string const &unit) const; ///< Get PolishDoughnut::central_density_ in specified unit
 void   centralDensity(double density); ///< Set PolishDoughnut::central_density_
 void   centralDensity(double density, std::string const &unit); ///< Set PolishDoughnut::central_density_ in specified unit

 double centralTempOverVirial() const; ///< Get PolishDoughnut::centraltemp_over_virial_
 void   centralTempOverVirial(double val); ///< Set PolishDoughnut::centraltemp_over_virial_

 double beta() const; ///< Get PolishDoughnut::beta_
 void   beta(double beta);///< Set PolishDoughnut::beta_

 void   spectralOversampling(size_t); ///< Set PolishDoughnut::spectral_oversampling_
 size_t spectralOversampling() const ; ///< Get PolishDoughnut::spectral_oversampling_

 bool changeCusp() const; ///< Get PolishDoughnut::komissarov_
 void changeCusp(bool t); ///< Set PolishDoughnut::komissarov_
 bool komissarov() const; ///< Get PolishDoughnut::komissarov_
 void komissarov(bool komis); ///< Set PolishDoughnut::komissarov_
 bool angleAveraged() const; ///< Get PolishDoughnut::angle_averaged_

 /**
  * if komis, also sets komissarov_ to true
  */
 void angleAveraged(bool komis); ///< Set PolishDoughnut::angle_averaged_

 void nonThermalDeltaExpo(std::vector<double> const &v);
 std::vector<double> nonThermalDeltaExpo() const;
 void angmomrinner(std::vector<double> const &v);
 std::vector<double> angmomrinner() const;
 void adafparams(std::vector<double> const &v);
 std::vector<double> adafparams() const;
 void adaf(bool t);
 bool adaf() const;
 void setParameter(Gyoto::Property const &p,
		   std::string const & name,
		   std::string const & content,
		   std::string const & unit);


 // Read only members, depend on lambda
 double getWsurface() const; ///< Get PolishDoughnut::W_surface_
 double getWcentre() const; ///< Get PolishDoughnut::W_centre_
 double getRcusp() const; ///< Get PolishDoughnut::r_cusp_
 double getRcentre() const; ///< Get PolishDoughnut::r_centre_

 using Standard::metric;
 virtual void metric(Gyoto::SmartPointer<Gyoto::Metric::Generic>);

 // ASTROBJ API
 // -----------
 int Impact(Photon *ph, size_t index,
			    Astrobj::Properties *data);

 virtual double operator()(double const coord[4]) ;

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

  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double nu_em[], size_t nbnu,
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
					      )  const;
  double emissionSynchro_komissarov_averaged(double Theta_elec, 
					     double number_density,
					     double nuem,
					     double nuc 
					     ) const;
  double emissionSynchro_komissarov_averaged_integ(double Theta_elec, 
						   double number_density,
						   double nuem,
						   double nuc 
						   ) const;
  double emissionSynchro_komissarov_PL_direction(
						 double number_density_PL,
						 double nuem, double nuc,
						 double theta_mag) const;
  double emissionSynchro_komissarov_PL_averaged(
						 double number_density_PL,
						 double nuem, double nuc
						) const;
  double absorptionSynchro_komissarov_PL_direction(
						   double number_density_PL,
						   double nuem, double nuc,
						   double theta_mag) const ;
  double absorptionSynchro_komissarov_PL_averaged(
						  double number_density_PL,
						  double nuem, double nuc
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
    const PolishDoughnut * papa;
    virtual double operator() (double) const;
  };
 intersection_t intersection; ///< double intersection(double) Functor

/**
  * \class Gyoto::Astrobj::PolishDoughnut::outerradius_t
  * \brief double outerradius(double) Functor class
  *
  * To find the outer doughnut radius.
  * This class is as a local variable in PolishDoughnut::emission()
  */
  class outerradius_t : public Gyoto::Functor::Double_Double_const {
  public:
    const PolishDoughnut * papa;
    virtual double operator() (double) const;
  };

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
