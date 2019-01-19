/**
 * \file GyotoPolishDoughnut.h
 * \brief A magnetized toroidal accretion structure
 *
 *  Latest reference: Vincent, F. H.; Yan, W.; Straub, O.;
 *  Zdziarski, A. A.; Abramowicz, M. A. 2015, <STRONG>A magnetized torus
 *  for modeling Sagittarius A* millimeter images and spectra</STRONG>,
 *  A&amp;A 574:A48.
 *
 *  First reference: Straub, O.; Vincent, F. H.; Abramowicz, M. A.;
 *  Gourgoulhon, E.; &amp; Paumard, T. 2012, <STRONG>Modelling the
 *  black hole silhouette in Sagittarius A* with ion tori</STRONG>,
 *  A&amp;A 543:83.
 *
 */

/*
    Copyright (c) 2012-2016, 2018-2019 Frederic Vincent, Odele Straub,
                                       Thibaut Paumard

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
#include <GyotoThermalBremsstrahlungSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>
#include <GyotoPowerLawSynchrotronSpectrum.h>

/**
 * \brief A toroidal accretion structure
 *
 *  Latest reference: Vincent, F. H.; Yan, W.; Straub, O.;
 *  Zdziarski, A. A.; Abramowicz, M. A. 2015, <STRONG>A magnetized torus
 *  for modeling Sagittarius A* millimeter images and spectra</STRONG>,
 *  A&amp;A 574:A48.
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
 SmartPointer<Spectrum::ThermalBremsstrahlung> spectrumBrems_;
 SmartPointer<Spectrum::ThermalSynchrotron> spectrumSynch_;
 SmartPointer<Spectrum::PowerLawSynchrotron> spectrumPLSynch_;
 double l0_; ///< Angular momentum. Tied to PolishDoughnut::lambda_.
 double lambda_; ///< Adimentionned angular momentum
 double W_surface_; ///< Potential surface value. Tied to PolishDoughnut::lambda_.
 double W_centre_; ///< Potential central value. Tied to PolishDoughnut::lambda_.
 double r_cusp_; ///< Cusp radius in geometrical units. Tied to PolishDoughnut::lambda_.
 double r_centre_; ///< Central radius in geometrical units. Tied to PolishDoughnut::lambda_.
 double r_torusouter_ ; ///< Torus outer coordinate radius. Tied to PolishDoughnut::lambda_.
 double DeltaWm1_; ///< 1./(W_centre_ - W_surface_);
 double central_enthalpy_cgs_; ///< Central enthalpy per unit volume in erg/cm3
 double central_temperature_; ///< T<SUB>center</SUB> in K
 double beta_; ///< P<SUB>gas</SUB>/P<SUB>magn</SUB> (careful not standard)
 double magnetizationParameter_; ///< P<SUB>magn</SUB>/(n<SUB>e</SUB> m<SUB>p</SUB> c<SUP>2</SUP>) (careful, very different from above)
 double aa_; ///< PolishDoughnut::gg_ spin, cached when setting PolishDoughnut::lambda_
 double aa2_; ///< aa_<SUP>2</SUP>
 size_t spectral_oversampling_;///< Oversampling used in integrateEmission()
 bool angle_averaged_; ///< 1 if Komissarov model should be angle averaged
 bool bremsstrahlung_; ///< 1 if Komissarov model should compute Brems radiation

 /// fraction of thermal energy in non-thermal electrons
 /**
  * Obsiously, 0 means no non-thermal electrons at all. In this case,
  * the (trivial) non-thermal computations are skipped. Ther is thus
  * non need for an additional "nonthermal_" flag.
  */
 double deltaPL_; ///< Fraction of electrons emitting powerlaw synchro

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
 GYOTO_OBJECT_THREAD_SAFETY;
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
 // Accessors
 // ---------
public:
 double getL0() const; ///< Get PolishDoughnut::l0_
 // void   setL0(double l0); set by lambda_

 double lambda() const; ///< Get PolishDoughnut::lambda_
 void   lambda(double lambda); ///< Set PolishDoughnut::lambda_

 double centralEnthalpyPerUnitVolume() const; ///< Get PolishDoughnut::central_enthalpy_cgs_
 double centralEnthalpyPerUnitVolume(std::string const &unit) const; ///< Get PolishDoughnut::central_enthalpy_cgs_ in specified unit
 void   centralEnthalpyPerUnitVolume(double density); ///< Set PolishDoughnut::central_enthalpy_cgs_
 void   centralEnthalpyPerUnitVolume(double density, std::string const &unit); ///< Set PolishDoughnut::central_enthalpy_cgs_ in specified unit

 double centralTemp() const; ///< Get PolishDoughnut::central_temperature_
 void   centralTemp(double val); ///< Set PolishDoughnut::central_temperature_

 double beta() const; ///< Get PolishDoughnut::beta_
 void   beta(double beta);///< Set PolishDoughnut::beta_

 void magnetizationParameter(double rr);
 double magnetizationParameter()const;

 void   spectralOversampling(size_t); ///< Set PolishDoughnut::spectral_oversampling_
 size_t spectralOversampling() const ; ///< Get PolishDoughnut::spectral_oversampling_

 bool changeCusp() const; ///< Get PolishDoughnut::komissarov_
 void changeCusp(bool t); ///< Set PolishDoughnut::komissarov_
 bool bremsstrahlung() const; ///< Get PolishDoughnut::bremsstrahlung_
 void bremsstrahlung(bool brems); ///< Set PolishDoughnut::bremsstrahlung_
 bool angleAveraged() const; ///< Get PolishDoughnut::angle_averaged_
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
  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, state_t const &cph, double const *co) const;

  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;

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
  * This class is instantiated in a single
  * PolishDoughnut::intersection member.
  */
  class intersection_t : public Gyoto::Functor::Double_Double_const {
  public:
    intersection_t(PolishDoughnut* parent);
    PolishDoughnut * papa;
    virtual double operator() (double) const;
  };
  friend class intersection_t;

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

 // Outputs
 // -------
 public:

 /// Display
 friend std::ostream& operator<<(std::ostream& , const PolishDoughnut& ) ;

 public:


};

#endif
