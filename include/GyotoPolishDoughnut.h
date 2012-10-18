/**
* \file GyotoPolishDoughnut.h
* \brief A toroïdal accretion structure
*
*  Beware: the PolishDoughnut currently uses c.g.s. units, it _will_
*  be converted to use S.I. soon.
*
*/

/*
*   Copyright (c) year  your_name
*
 *
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
//#include <GyotoPolishDoughnutCst.h>

/**
* \class Gyoto::PolishDoughnut
* \brief 
* 
 */
class Gyoto::Astrobj::PolishDoughnut : public Astrobj::Standard {
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
 double central_density_; ///< central density in g cm^-3
 double centraltemp_over_virial_; ///< T_center/Tvirial
 double beta_; ///< Pmagn/Pgas
 int use_specific_impact_ ;///< PolishDoughnut::Impact_() or Standard::Impact()
 double aa_; ///< Metric spin, cached when setting lambda_
 double aa2_; ///< aa_^2

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
 void   setCentralDensity(double density);

 double getCentralTempOverVirial() const;
 void   setCentralTempOverVirial(double val);

 double getBeta() const;
 void   setBeta(double beta);

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

 virtual int setParameter(std::string name, std::string content) ;
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
#endif

  // ASTROBJ processHitQuantities API
  // --------------------------------
protected:
  virtual void getVelocity(double const pos[4], double vel[4]) ;

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
  double funcxM(double alpha1, double alpha2, double alpha3, double xM) const;

 // PURELY INTERNAL TO ASTROBJ
 // --------------------------
private:
 double potential(double r, double theata) const;

 double intersection(double r) ;

 double transcendental(double xM, double par[4]) const ;

 double bisection_transcendental_neg(double param[4], double r_min, double r_max) const ;

 double bisection_intersection_neg(double r1_min, double r1_max) ;

 double bisection_intersection_pos(double r2_min, double r2_max) ;

 double bessi0(double xx) const;
 double bessi1(double xx) const;
 double bessk0(double xx) const;
 double bessk1(double xx) const;
 double bessk(int nn, double xx) const;

 // Outputs
 // -------
 public:

 /// Display
 friend std::ostream& operator<<(std::ostream& , const PolishDoughnut& ) ;

 public:
 

};

#endif
