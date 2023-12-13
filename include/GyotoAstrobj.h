/**
 * \file GyotoAstrobj.h
 * \brief Astronomical objects (light emitters)
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2011-2020 Thibaut Paumard, Frederic Vincent

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


#ifndef __GyotoAstrobj_H_ 
#define __GyotoAstrobj_H_ 

#include "GyotoConfig.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>
#include <GyotoConverters.h>
#include <GyotoObject.h>
#include <eigen3/Eigen/Dense>

namespace Gyoto{
  class Photon;
  namespace Register { class Entry; }
  namespace Metric { class Generic; }
  class FactoryMessenger;
  namespace Astrobj {
    class Generic;
    class Properties;

    /**
     * This is a more specific version of the
     * SmartPointee::Subcontractor_t type. An Astrobj::Subcontrator_t
     * is called by the Gyoto::Factory to build an instance of the
     * kind of astronomical object specified in an XML file (see
     * Register()). The Factory and Subcontractor_t function
     * communicate through a Gyoto::FactoryMessenger. A template is
     * provided so that you may not have to code anything.
     */
    typedef SmartPointer<Gyoto::Astrobj::Generic>
      Subcontractor_t(Gyoto::FactoryMessenger*, std::vector<std::string> const &);
    ///< A function to build instances of a specific Astrobj::Generic sub-class
 
    /**
     * Instead of reimplementing the wheel, your subcontractor can simply be
     * Gyoto::Astrobj::Subcontractor<MyKind>.
     *
     * If MyKind accepts any XML parameters, it should re-implement
     * Astrobj::Generic::setParameter() or, if low-level access to the
     * FactoryMessenger is needed, Generic::setParameters().
     *
     * \tparam T Gyoto::Astrobj::Generic sub-class
     */
    template<typename T> SmartPointer<Astrobj::Generic> Subcontractor
      (FactoryMessenger* fmp, std::vector<std::string> const &plugin) {
      SmartPointer<T> ao = new T();
      ao -> plugins(plugin) ;
#ifdef GYOTO_USE_XERCES
      if (fmp) ao -> setParameters(fmp);
#endif
      return ao;
    }
    ///< A template for Subcontractor_t functions

    /// Query the Astrobj register
   /**
     * Query the Astrobj register to get the Astrobj::Subcontractor_t
     * corresponding to a given kind name. This function is normally
     * called only from the Factory. If plugin is specified, only a
     * subcontractor matching both name and plugin will be returned,
     * loading the plug-in if necessary. If plugin is the empty
     * string, then the first subcontractor matching name will be
     * returned, and the name of the plug-in it belongs to will be
     * returned in plugin upon output.
     *
     * \param[in] name e.g. "Star"
     * \param[inout] plugin e.g. "stdplug".
     * \param[in] errmode 1 to return NULL in case of failure instead of
     * throwing an Error.
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Astrobj::Subcontractor_t* getSubcontractor(std::string name,
						      std::vector<std::string> &plugin,
						      int errmode = 0);

    /**
     * Use the Astrobj::initRegister() once in your program to
     * initiliaze it, the Astrobj::Register() function to fill it, and
     * the Astrobj::getSubcontractor() function to query it.
     */
    extern Gyoto::Register::Entry * Register_;
    ///< The Astrobj register

     /**
      *  This must be called once.
      */
    void initRegister(); 
    ///< Empty the Astrobj register

    /**
     * Register a new Astrobj::Generic sub-class so that the
     * Gyoto::Factory knows it.
     *
     * \param name The kind name which identifies this object type in
     * an XML file, as in &lt;Astrobj kind="name"&gt;
     *
     * \param scp A pointer to the subcontractor, which will
     * communicate with the Gyoto::Factory to build an instance of
     * the class from its XML description
     */
    void Register(std::string name, Gyoto::Astrobj::Subcontractor_t* scp);
    ///< Make an Astrobj kind known to the Factory
  }
}

/**
 * \namespace Gyoto::Astrobj
 * \brief Access to astronomical objects 
 *
 *  Objects which are supposed to be the target of the ray-tracing
 *  code should inherit from the Gyoto::Astrobj::Generic class.
 *
 *  When implementing a new object, you must:
 *    - make sure the object can be loaded from XML by providing a
 *      Subcontractor_t using the Gyoto::Astrobj::Register(std::string
 *      name, Gyoto::Astrobj::Subcontractor_t* scp) function;
 *    - make sure this subcontractor is registerred in the initialization
 *      routine of your plug-in;
 *    - make sure  Generic::Impact() works (see below).
 *
 *  In addition, you should make sure that your object plays nicely in
 *  the Yorick plug-in, which means:
 *    - implement the copy constructor and the Generic::clone() method;
 *    - implement the fillElement method, used for printing and saving to
 *      XML.
 *
 *  There are basically two ways of making Generic::Impact() work:
 *  either by making the Astrobj a sub-class of the low-level
 *  Gyoto::Astrobj::Generic class ans providing your own
 *  Generic::Impact() function (which, in principle, should rely on
 *  Generic::processHitQuantities()), or by making the Astrobj a
 *  sub-class of the higher-level Gyoto::Astrobj::Standard class and
 *  implementing two lower level, simpler functions which are
 *  used by the Standard::Impact():
 *    - Standard::operator()() yields a distance or potential defining
 *      the interior of the object;
 *    - Standard::getVelocity() yields the velocity field of the fluid.
 *
 *  Generic::processHitQuantities() itself is an intermediate-level
 *  function which you may choose to reimplement. It uses three
 *  low-level, easy to implement functions:
 *    - Generic::emission();
 *    - Generic::integrateEmission();
 *    - Generic::transmission().
 *  Default implementations of these three functions exist, they have
 *  little physical relevance but allow quick 0-th order vizualisation
 *  of your object.
 *
 * To be usable, a Astrobj::Generic (or Astrobj::Standard) sub-classe
 * should register an Astrobj::Subcontractor_t function using the
 * Astrobj::Register() function. See also \ref writing_plugins_page
 * . If your clas implements setParameter() and/or, if necessary,
 * setParameters(), registering it is normally done using the provided
 * template:
 * \code
 * Astrobj::Register("MyKind", &(Astrobj::Subcontractor<Astrobj::MyKind>));
 * \endcode
 */
/**
 * \class Gyoto::Astrobj::Generic
 * \brief Base class for astronomical object
 *
 * See introduction in the Gyoto::Astrobj namespace.
 */
class Gyoto::Astrobj::Generic
: public Gyoto::SmartPointee,
  public Gyoto::Object
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Generic>;


  // Data : 
  // -----
private:

  /**
   * \brief Whether some virtual methods are implemented
   *
   * The default implementations of some methods call
   * one-another. This member is used internally to avoid infinite
   * recursion.
   */
  int __defaultfeatures;

 protected:

  /**
   * \brief The Metric in this end of the Universe
   */
  SmartPointer<Gyoto::Metric::Generic> gg_;


  /**
   * Maximum distance from the center of the coordinate system at
   * which a photon may hit the object.  Child classes may compute a
   * decent value for #rmax_ at any time if #rmax_ is
   * DBL_MAX. External classes (Photons in particular) must use rMax()
   * to access this information.
   *
   * #rmax_ is in geometrical units.
   */                                         
  double rmax_; ///< Maximum distance to the center of the coordinate system [geometrical units]

  /**
   * Maximum step for Photon integation when inside the sphere of
   * radius #rmax_.
   *
   * #deltamaxinsidermax_ is in geometrical units.
   */
  double deltamaxinsidermax_; ///< Maximum Photon integration step inside #rmax_ [geometrical units]

  bool flag_radtransf_; ///< 1 if radiative transfer inside Astrobj, else 0

  int shadow_; ///< 1 to highlight the shadow region in the image
  int noredshift_; ///< 1 to impose redshift factor g = 1
  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

  /**
   *  #kind_ =  "Default", #rmax_ = DBL_MAX
   */
  Generic(); ///< Default constructor.

  /**
   *  #kind_ =  "Default", #rmax_ = radmax
   */
  Generic(double radmax); ///< Set rmax in constructor.

  /**
   *  #kind_ =  kind, #rmax_ = DBL_MAX
   */
  Generic(std::string kind); ///< Set kind in constructor.

  /**
   * Make a deep copy of an Astrobj::Generic instance
   */
  Generic(const Generic& ) ; ///< Copy constructor.

  /**
   * This method must be implemented by the various Astrobj::Generic
   * subclasses in order to support cloning:
   * \code
   * SmartPointer<Astrobj> deep_copy = original->clone();
   * \endcode
   *
   * Cloning is necessary for multi-threading, recommended for
   * interaction with the Yorick plug-in etc.
   *
   * Implementing it is very straightforward, as long as the copy
   * constructor Generic(const Generic& ) has been implemented:
   * \code
   * MyAstrobj* MyAstrobj::clone() const { return new MyAstrobj(*this); }
   * \endcode
   */
  virtual Generic* clone() const = 0 ; ///< Cloner
  
  virtual ~Generic() ; ///< Destructor: does nothing.

  // Accessors
  // ---------
 public:
  /**
   * \brief Get the Metric #gg_
   */
  virtual SmartPointer<Metric::Generic> metric() const;

  /**
   * \brief Set the Metric #gg_
   */
  virtual void metric(SmartPointer<Metric::Generic>) ;

  /**
   *  Get maximal distance from center of coordinate system at which a
   *  Photon may hit the object.
   *  
   *  Child classes may use the #rmax_ member to cache this value, if
   *  its current value is DBL_MAX.
   *
   *  It can also be set using rMax().
   *
   *  \return rmax_ in geometrical units
   */
  virtual double rMax(); ///< Get maximal distance from center of coordinate system
  virtual double rMax() const; ///< Get maximal distance from center of coordinate system

  /**
   *  Call rMax() and convert result to unit.
   *
   *  \param unit string
   *  \return double rmax converted to unit
   */
  virtual double rMax(std::string const &unit); ///< Get rmax_ is specified unit
  virtual double rMax(std::string const &unit) const; ///< Get rmax_ is specified unit

  /// Get max step constraint for adaptive integration
  /**
   * \param[in] coord position
   * \return max step to find this object reliably
   */
  virtual double deltaMax(double coord[8]);

  /**
   *  Set maximal distance from center of coordinate system at which a
   *  Photon may hit the object.
   *  
   *  \param val new #rmax_ in geometrical units.
   */
  virtual void rMax(double val); ///< Set maximal distance from center of coordinate system

  /**
   *  Call Generic::rMax(double val) after converting val from unit
   *  to geometrical units.
   *
   *  \param val #rmax_ expressed in unit "unit";
   *  \param unit string...
   */
  virtual void rMax(double val, std::string const &unit); ///< Set maximal distance from center of coordinate system

  GYOTO_OBJECT_ACCESSORS_UNIT(deltaMaxInsideRMax);

  /**
   * Set flag indicating that radiative transfer should be integrated,
   * i.e. the object is to be considered optically thin.
   * \param flag: 1 if optically thin, 0 if optically thick.
   */
  void opticallyThin(bool flag);
  ///< Set whether the object is optically thin.
  /**
   * See opticallyThin(bool flag).
   */
  bool opticallyThin() const ;
  ///< Query whether object is optically thin.

  void showshadow(bool flag);
  bool showshadow() const ;

  void redshift(bool flag);
  bool redshift() const ;

  /**
   * Return a Gyoto::Quantity_t suitable as input to
   * Gyoto::Scenery::setRequestedQuantities() to set de default
   * quantities to compute for this object. The default of these
   * defaults GYOTO_QUANTITY_INTENSITY.
   */
  virtual Gyoto::Quantity_t getDefaultQuantities();
  ///< Which quantities to compute if know was requested

  //XML I/O
 public:

#ifdef GYOTO_USE_XERCES
  /**
   * \brief Main loop in Subcontractor_t function
   *
   * The Subcontractor_t function for each Astrobj kind should look
   * somewhat like this (templated as
   * Gyoto::Astrobj::Subcontractor<MyKind>):
   * \code
   * SmartPointer<Astrobj::Generic>
   * Gyoto::Astrobj::MyKind::Subcontractor(FactoryMessenger* fmp) {
   *   SmartPointer<MyKind> ao = new MyKind();
   *   ao -> setParameters(fmp);
   *   return ao;
   * }
   * \endcode
   *
   * Each object kind should implement setParameter(string name,
   * string content, string unit) to interpret the individual XML
   * elements. setParameters() can be overloaded in case the specific
   * Astrobj class needs low level access to the FactoryMessenger. See
   * UniformSphere::setParameters().
   */
  virtual void setParameters(FactoryMessenger *fmp);


#endif
  
  // Outputs
  // -------
 public:
  /**
   * Impact() checks whether a Photon impacts the object between two
   * integration steps of the photon's trajectory (those two steps are
   * photon->getCoord(index, coord1) and photon->getCoord(index+1,
   * coord2)). Impact returns 1 if the photon impacts the object
   * between these two steps, else 0. In many cases of geometrically
   * thick obects, the implementation Astrobj::Standard::Impact() will
   * be fine.
   *
   * Impact will call Generic::processHitQuantities() (which is
   * virtual and may be re-implemented) to compute observable
   * properties on demand: if the data pointer is non-NULL, the object
   * will look in it for pointers to properties which apply to its
   * kind. If a pointer to a property known to this object is present,
   * then the property is computed and store at the pointed-to
   * address. For instance, all objects know the "intensity"
   * property. If data->intensity != NULL, the instensity is computed
   * and stored in *data->intensity.
   *
   * If data is non-NULL and only in this case, processHitQuantities()
   * will also call ph->transmit() to update the transmissions of the
   * Photon (see Photon::transmit(size_t, double)). This must not be
   * done if data is NULL (see Astrobj::Complex::Impact() for an
   * explanation).
   *
   * Impact() may not extend the ph Worldline. The only two dates that
   * are guaranteed to be defined are at indices index and index+1.
   *
   * \param ph   Gyoto::Photon aimed at the object;
   * \param index    Index of the last photon step;
   * \param data     Pointer to a structure to hold the observables at impact.
   *
   * \return 1 if impact, 0 if not.
   */
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL) = 0 ;
  ///< Does a photon at these coordinates impact the object?
  
  /**
   * \brief Fills Astrobj::Properties
   *
   * processHitQuantities fills the requested data in Impact. To use
   * it, you need to call it in the Impact() method for your object in
   * case of hit. It will fill Redshift, Intensity, Spectrum,
   * BinSpectrum and update the Photon's transmission by calling
   * Photon::transmit(), only if data==NULL.
   *
   * You can overload it for your Astrobj. The generic implementation
   * calls emission(), integrateEmission() and transmission() below.
   */
  virtual void processHitQuantities(Photon * ph, state_t const &coord_ph_hit,
                                   double const * coord_obj_hit, double dt,
                                   Astrobj::Properties* data) const;
  virtual void processHitQuantities(Photon* ph, double* coord_ph_hit,
                                   double* coord_obj_hit, double dt,
                                   Astrobj::Properties* data) const = delete ;

  /**
   * \brief Specific intensity I<SUB>&nu;</SUB>
   *
   * Called by the default implementation for processHitQuantities().
   *
   * emission() computes the intensity I<SUB>&nu;</SUB> emitted by the
   * small volume of length ds<SUB>em</SUB>, in the emitter's
   * frame. It should take self-absorption along ds<SUB>em</SUB> into
   * account.
   *
   * Reminder :
   *  - intensity = I<SUB>&nu;</SUB> [J m^-2 s^-1 ster^-1 Hz^-1];
   *
   *  - invariant intensity = I<SUB>&nu;</SUB>/&nu;<SUP>3</SUP>, which
   *    has the same value in any frame;
   *
   *  - emission coefficient = j<SUB>&nu;</SUB> [J m^-3 s^-1 ster^-1
   *    Hz^-1] , defined by dI<SUB>&nu;</SUB> = j<SUB>&nu;</SUB>*ds,
   *    where ds is the distance travelled by the photon inside the
   *    object;
   *  - invariant emission coef = j<SUB>&nu;</SUB>/&nu;<SUP>2</SUP>,
   *    which has the same value in any frame.
   *
   * The equation used for radiative transfer (without absorption) is:
   *
   *    d(I<SUB>&nu;</SUB>/&nu;<SUP>3</SUP>)/d&lambda; = (j<SUB>&nu;</SUB>/&nu;<SUP>2</SUP>)  [*]
   *
   * where &lambda; is the integration parameter along the null geodesic.
   *
   * NB: Let us consider a particular observer, with &nu; being the
   * frequency measured by this observer, and ds being the proper
   * distance (as measured by the observer) that the photon travels
   * as it moves from &lambda; to &lambda;+d&lambda; along its
   * geodesic.  Then it can be shown that:
   *
   *    d&lambda; = ds/&nu;
   *
   * This shows that Eq. [*] is homogeneous.
   *
   * The default implementation returns 1. if optically thick and ds<SUB>em</SUB>
   * if optically thin. It allows for a quick implementation of your
   * object for visualization purposes.
   *
   * \param nu_em Frequency at emission [Hz]
   * \param dsem length over which to integrate inside the object
   *        [geometrical units]
   * \param coord_ph Photon coordinate
   * \param coord_obj Emitter coordinate at current photon position
   */
  virtual double emission(double nu_em, double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL)
    const ;
  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL)
    const = delete; ///< Obsolete, update your code;

  /**
   * \brief Specific intensity I<SUB>&nu;</SUB> for several values of &nu;<SUB>em</SUB>
   *
   * Called by the default implementation for processHitQuantities().
   *
   * emission() computes the intensity I<SUB>&nu;</SUB> emitted by the small
   * volume of length dsem. It should take self-absorption along dsem
   * into account.
   *
   * Same as emission(double nu_em, double dsem, double coord_ph[8],
   *		  double coord_obj[8]=NULL) const
   * looping on several values of nu_em.
   *
   * \param Inu[nbnu] Output (must be set to a previously allocated
   *        array of doubles)
   * \param nu_em[nbnu] Frequencies at emission
   * \param nbnu Size of Inu[] and nu_em[] 
   * \param dsem Length over which to integrate inside the object
   * \param coord_ph Photon coordinate
   * \param coord_obj Emitter coordinate at current photon position
   * \return I<SUB>&nu;</SUB> or dI<SUB>&nu;</SUB> [W m-2 sr-2]
   */
  virtual void emission(double Inu[], double const nu_em[], size_t nbnu,
			double dsem, state_t const &coord_ph,
			double const coord_obj[8]=NULL) const ;
  virtual void emission(double Inu[], double nu_em[], size_t nbnu,
			double dsem, double coord_ph[8],
			double coord_obj[8]=NULL) const = delete ; ///< Obsolete, update your code

  /**
   * \brief emission and transmission together
   */
  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;
  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double nu_em[], size_t nbnu,
			  double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL) const = delete ;
  /**
   * \brief Compute the increment of Stokes parameters and transmission matrix. Polarised version of RadiaveQ
   * 
   * First function to be called for radiative quantities. 
   * If exist, i.e. implemented in an Astrobj, return the Stokes parameters emitted by the small
   * volume of length dsem.
   * 
   * Warning : 
   *  - The basis used to determine the Stokes coefficients is different from the observer parallel transported polarisation basis.
   *    One should use getChi function to compute the angle between these two basis.
   * 
   *  - The non polarized case must also be implemented in this function to avoid error.
   * 
   * See exemple in SimplePolarStar.C.
   * 
   * 
   * \param Inu[nbnu] Output increment of intensity (must be set to a previously allocated
   *        array of doubles)
   * \param Qnu[nbnu] Output increment of Stokes parameter Q (must be set to a previously allocated
   *        array of doubles)
   * \param Unu[nbnu] Output increment of Stokes parameter U (must be set to a previously allocated
   *        array of doubles)
   * \param Vnu[nbnu] Output increment of Stokes parameter V (must be set to a previously allocated
   *        array of doubles)
   * \param Onu[nbnu] Output transmission (Eigen) matrix (must be set to a previously allocated
   *        array of Matrix4d)
   * \param nu_em[nbnu] Frequencies at emission
   * \param nbnu Size of nu_em[], Inu[], Qnu[], Unu[], Vnu[], and Onu[]
   * \param dsem Length over which to integrate inside the object
   * \param cph Photon coordinate
   * \param co Emitter coordinate at current photon position
   * \return Increment of the Stokes parameters (I,Q,U,V) and local Transmission matrix (O).
   */
  virtual void radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
			  Eigen::Matrix4d *Onu, double const *nuem , size_t nbnu, double dsem,
			  state_t const &cph, double const *co) const ;

  /**
   * Compute the integral of emission() from &nu;<SUB>1</SUB> to
   * &nu;<SUB>2</SUB>. The default implementation is a numerical
   * integrator which works well enough and is reasonably fast if
   * emission() is a smooth function (i.e. no emission or absorption
   * lines). If possible, it is wise to implement an analytical
   * solution. It is used by processHitQuantities to compute the
   * "BinSpectrum" quantity which is the most physical: it is the only
   * quantity that can be actually measured directly by a real-life
   * instrument.
   */
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   state_t const &c_ph, double const c_obj[8]=NULL) const;
    ///< &int;<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP> I<SUB>&nu;</SUB> d&nu; (or j<SUB>&nu;</SUB>)
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const=delete; ///< Obsolete, update your code

  /**
   * Like double integrateEmission(double nu1, double nu2, double
   * dsem, double c_ph[8], double c_obj[8]) const for each
   * Spectrometer channel.
   */
  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, state_t const &cph, double const *co) const;
    ///< &int;<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP> I<SUB>&nu;</SUB> d&nu; (or j<SUB>&nu;</SUB>)
  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, double *cph, double *co) const = delete;
  ///< Obsolete, update your code

  /**
   * transmission() computes the transmission of this fluid element or
   * 0 if optically thick. The default implementation returns 1. (no
   * attenuation) if optically thin, 0. if optically thick.
   *
   * \param nuem frequency in the fluid's frame
   * \param coord Photon coordinate
   * \param dsem geometrical length in geometrical units
   */
  virtual double transmission(double nuem, double dsem, state_t const &coord_ph, double const coord_obj[8]) const ;
     ///< Transmission: exp( &alpha;<SUB>&nu;</SUB> * ds<SUB>em</SUB> )
  virtual double transmission(double nuem, double dsem, state_t const &coord) const = delete;  ///< Obsolete, update your code
  virtual double transmission(double nuem, double dsem, double coord[8]) const = delete;
  ///< Obsolete, update your code

  /**
   * Omatrix() computes the exponential of the Mueller matrix 
   * which contains the absorption and Faraday coefficients 
   * and is used in the polarized radiative transfer equation.
   * 
   * \param Onu output matrix (must be allocated previously)
   * \param alphanu[4] array containing the 4 absorption coefficients in the Stokes basis (I,Q,U,V)
   * \param rnu[3] array containing the 3 Faraday coefficients in the Stokes basis (Q,U,V)
   * \param Chi angle of rotation between the parallel transported observer polarization basis and the Stokes basis
   * \param dsem geometrical length in geometrical units
   */
  Eigen::Matrix4d Omatrix(double alphanu[4], double rnu[3], double Chi, double dsem) const;
  Eigen::Matrix4d Omatrix(double alphaInu, double alphaQnu, double alphaUnu, double alphaVnu,
        double rQnu, double rUnu, double rVnu, double Chi, double dsem) const;
  Eigen::Matrix4d Omatrix(double alphanu[4], double rnu[3], double sin2Chi, double cos2Chi, double dsem) const;
  Eigen::Matrix4d Omatrix(double alphaInu, double alphaQnu, double alphaUnu, double alphaVnu,
        double rQnu, double rUnu, double rVnu, double sin2Chi, double cos2Chi, double dsem) const;

  Eigen::Matrix4d Pmatrix(double alphaInu, double alphaQnu, double alphaUnu, double alphaVnu,
        double rQnu, double rUnu, double rVnu, double sin2Chi, double cos2Chi, double dsem) const;
  
  /**
   * Apply the rotation matrix with angle Chi to the emission Stokes vector
   * constructed in the fonction from the individual coefficients
   */
  Eigen::Vector4d rotateJs(double jInu, double jQnu, double jUnu, double jVnu, double sin2Chi, double cos2Chi) const;
  Eigen::Vector4d rotateJs(double jInu, double jQnu, double jUnu, double jVnu, double Chi) const;

  /**
   * Get Chi angle.
   * Return the angle between the parallel transported observer polarization basis (Ephi,Etheta)
   * and the Stokes basis in the rest frame of the emitter defined by the 4-vector magnetic/electric field.
   * 
   * \param fourvect 4-vector magnetic/electric field depending on elec (false/true)
   * \param cph Photon coordinate, must contain the Ephi and Etheta vectors i.e. size(cph)==16
   * \param vel Fluid velocity at the photon coordinate
   */
  double getChi(double const fourvect[4], state_t const &cph, double const vel[4], bool elec=false) const;
  /**
   * Get the cosinus and sinus of 2*Chi angle.
   * Chi being the angle between the parallel transported observer polarization basis (Ephi,Etheta)
   * and the Stokes basis in the rest frame of the emitter defined by the 4-vector magnetic/electric field.
   * 
   * \param fourvect 4-vector magnetic/electric field depending on elec (false/true)
   * \param cph Photon coordinate, must contain the Ephi and Etheta vectors i.e. size(cph)==16
   * \param vel Fluid velocity at the photon coordinate
   */
  void getSinCos2Chi(double const fourvect[4], state_t const &cph, double const vel[4], double* sin2Chi, double* cos2Chi, bool elec=false) const;
  /**
   * Compute the magnetic field 4-vector at a given position and a given configuration
   * 
   * \param B4vect the output 4-vector magnetic field
   * \param magneticConfig string which determine the wanted configuration ('Radial', 'Azimuthal', 'Vertical')
   * \param co coordinate of the object which also contain its velocity
   * \param cph coordinate of the photon
   */
  void computeB4vect(double B4vect[4], std::string const magneticConfig, double const co[8], state_t const &cph) const;

  void computeB4vect_ipole(double B4vect[4], std::string const magneticConfig, double const co[8], state_t const &cph, double spin) const;
  /**
   * N-dimensional linear interpolation function
   * 
   * \param N number of dimensions
   * \param array array that should be interpolated. This array has to be a 1D array with the first dimension evolving the slowest
   * \param Xq query position at which the interpolation should be done. 1D array with N elements.
   * \param X 2D array containing the range of all the axes (for exemple : radius, \theta, \varphi) at which the array is knwon. Should be built as a pointer of pointers.
   * \param X_params 1D array that contain the length of each axes.
   * \param cond_limits 1D array of strings that set the boundary conditions that will be applied to each axes
   */
  double interpolate(int const N, double* const array, double* const Xq, double** const X, int* const X_params, std::string const *cond_limits) const;
  /**
   * N-dimensional linear interpolation function
   * 
   * \param N number of dimensions
   * \param array array that should be interpolated. This array has to be a 1D array with the first dimension evolving the slowest
   * \param Xq query position at which the interpolation should be done. 1D array with N elements.
   * \param X_params[N][3] 2D array that contain the minimum value, the maximum value and the length of each axes. The axes will be assumed to be linearly spaced.
   * \param cond_limits 1D array of strings that set the boundary conditions that will be applied to each axes.
   */
  double interpolate(int const N, double* const array, double* const Xq, double** const X_params, std::string const *cond_limits) const;

  private :
  /**
   * Function which returns the 1D linear interpolation.
   * 
   * \param x the x query position along the linear function [0,1].
   * \param y0 the value of the array to be interpolated at the low position.
   * \param y1 the value of the array to be interpolated at the high position.
   */
  double interp1d(double const x, double const y0, double const y1) const;
  /**
   * Function which returns the N- dimension interpolation on an N-dimensional cube
   * 
   * \param N number of dimensions
   * \param Xq query position at which the interpolation should be done. 1D array with N elements.
   * \param X[2**N][N] 2D array that contains for the 2**N points of the N dimensional cube the associated axes values.
   * \param Y[2**N] 1D array that contain the value of the array to be interpolated at each points of the N-dimensional cube.
   * \param cond_limit 1D array of strings that set the boundary conditions that will be applied to each axes.
   */
  double interpNd(int const N, double* const Xq, double** const X, double* const Y, std::string const *cond_limit) const;
  /**
   * Function that return the closest indice of the query value on an axis.
   * 
   * \param xq query value on the axis. Could be modified in case of periodic boundary condition.
   * \param cond_limit boundary conditions that will be applied to the axis.
   * \param X_params[3] array which contain the minimum value, the maximum value and the length of the axe.
   * \param X (optional) array which contains the value of the axis points. Usefull for not linearly spaced axis.
   */
  int getIndice(double &xq, std::string const cond_limit, double const X_params[3], double* const X=NULL) const;

};

/**
 * \class Gyoto::Astrobj::Properties
 * \brief Observable properties of an Astronomical object
 *
 *  The sort of properties one wants to measure on a ray-traced
 *  Gyoto::Photon which hits a Gyoto::Astrobj. Not all Astrobj are
 *  able to fill all of these properties.
 *
 *  An instance of Properties essentially contains a bunch of pointers
 *  to memory areas where the observable quantities (see Quantity_t)
 *  should be stored.
 *
 *  Astrobj::Generic::processHitQuantities() fills the various arrays
 *  upon request.  A quantity is ignored if the corresponding pointer
 *  is NULL.
 *
 *  Scenery::operator()() increments the Properties between each
 *  Photon using Properties::operator++().
 *
 *  The main application (gyoto, the yorick plug-in, or your user
 *  application) is responsible for allocating the various arrays,
 *  filling the various members of Properties, and doing whatever
 *  meaninful with the arrays after they have been filled with values
 *  by the ray-tracing code (e.g. saving them to disk or displaying
 *  them).
 *
 *  Also see Gyoto::Scenery and Gyoto::Quantity_t.
 */
class Gyoto::Astrobj::Properties : public Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Properties>;
 public:
  double *intensity; ///< GYOTO_QUANTITY_INTENSITY   : Intensity
  double *time; ///< GYOTO_QUANTITY_EMISSIONTIME: EmissionTime

  /**
   * Behaves like the square of the closest distance between Photon
   * and Astrobj (but not exactly that). Initialize it to DBL_MAX from
   * float.h.;
   */
  double *distance; ///< GYOTO_QUANTITY_MIN_DISTANCE: MinDistance

  /**
   * First local minimum in distance from object
   */
  double * first_dmin; ///< GYOTO_QUANTITY_FIRST_DMIN  : FirstDmin

  /**
   * Properties::first_dmin will be set to the first local minimum and
   * Properties::first_dmin_found will be set to 1 if a local minimum
   * in distance is found. Initialize it to 0.
   */
  int first_dmin_found; ///< Whether Properties::first_dmin was found

  /**
   * Redshift factor &nu;<SUB>obs</SUB>/&nu;<SUB>em</SUB> (necessary
   * for emission lines computation)
   */
  double *redshift; ///< GYOTO_QUANTITY_REDSHIFT    : RedShift

  /**
   * Numbere of times the geodesic crosses the equatorial plane
   */
  double *nbcrosseqplane; ///< GYOTO_QUANTITY_NBCROSSEQPLANE    : nuùber of equatorial plane crossings

  /**
   * I<SUB>&nu;</SUB> (&nu;) (observed specific intensity)
   */
  double *spectrum; ///< GYOTO_QUANTITY_SPECTRUM    : Spectrum

  double *stokesQ;
  double *stokesU;
  double *stokesV;

  /**
   *  I<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP>, the
   *  integral of I<SUB>&nu;</SUB> over each spectral channel
   *  (i.e. what a spectrometer would measure)
   */
  double *binspectrum; ///< GYOTO_QUANTITY_BINSPECTRUM : BinSpectrum

  /**
   *  Spectra elements are separated by offset doubles in memory. In
   *  other words, the ith spectral element is spectrum[i*offset].
   */
  ptrdiff_t offset; ///< How to jump from one spectral element to the next

  /**
   * Coordinates of the object and photon at impact
   */
  double * impactcoords; ///< GYOTO_QUANTITY_IMPACTCOORDS: ImpactCoords

  /**
   * \brief GYOTO_QUANTITY_USER1       : User1
   * Astrobj-specific quantity
   */
  double *user1;

  /**
   * \brief GYOTO_QUANTITY_USER2       : User2
   * Astrobj-specific quantity
   */
  double *user2;

  /**
   * \brief GYOTO_QUANTITY_USER3       : User3
   * Astrobj-specific quantity
   */
  double *user3;

  /**
   * \brief GYOTO_QUANTITY_USER4       : User4
   * Astrobj-specific quantity
   */
  double *user4;

  /**
   * \brief GYOTO_QUANTITY_USER5       : User5
   * Astrobj-specific quantity
   */
  double *user5;
# ifdef HAVE_UDUNITS
  /**
   * \brief Converter between SI (J.m <SUP> -2</SUP>.s<SUP>-1</SUP>.sr<SUP>-1</SUP>.Hz<SUP>-1</SUP>) and requested Intensity unit
   */
  Gyoto::SmartPointer<Gyoto::Units::Converter> intensity_converter_ ;
  /**
   * \brief Converter between SI (J.m <SUP> -2</SUP>.s<SUP>-1</SUP>.sr<SUP>-1</SUP>.Hz<SUP>-1</SUP>) and requested Spectrum unit
   */
  Gyoto::SmartPointer<Gyoto::Units::Converter> spectrum_converter_ ;
  /**
   * \brief Converter between SI (J.m <SUP> -2</SUP>.s<SUP>-1</SUP>.sr<SUP>-1</SUP>) and requested BinSpectrum unit
   */
  Gyoto::SmartPointer<Gyoto::Units::Converter> binspectrum_converter_ ;
# endif

  /// True if buffers are allocated for entire field (npix*npix)
  bool alloc;

 public:
  Properties(); ///< Default constructor (everything is set to NULL);
  Properties (double*, double*); ///<< Set intensity and time pointers.

  /**
   * \brief Initialize observable quantities
   *
   * The pointed-to values are initialized as follows (if the
   * corresponding pointer is not NULL):
   *
   * - intensity, firt_dmin_found, redshift, userN: 0
   * - time, distance, first_dmin: DBL_MAX
   * - for spectrum and binspectrum, nbnuobs values separated by offset in memory are initialized to 0
   * - for impactcoords, 16 contiguous values are initialized to DBL_MAX
   */
  void init(size_t nbnuobs=0);

  /**
   * \brief Increment pointers
   *
   * All valid pointers are incremented by 1 (sizeof(double)), excepted
   * impactcoords which is incremented by 16.
   */
  Properties& operator++();

  /**
   * \brief Increment pointers by offset
   *
   * All valid pointers are incremented by offset (sizeof(double)), excepted
   * impactcoords which is incremented by 16*offset.
   */
  Properties& operator+=(ptrdiff_t offset);

  operator Gyoto::Quantity_t () const;

# ifdef HAVE_UDUNITS
  void intensityConverter(Gyoto::SmartPointer<Gyoto::Units::Converter>);
  ///< Set Properties::intentity_converter_
  void intensityConverter(std::string);
  ///< Set Properties::intentity_converter_
  void spectrumConverter(Gyoto::SmartPointer<Gyoto::Units::Converter>);
  ///< Set Properties::spectrum_converter_
  void spectrumConverter(std::string);
  ///< Set Properties::spectrum_converter_
  void binSpectrumConverter(Gyoto::SmartPointer<Gyoto::Units::Converter>);
  ///< Set Properties::binspectrum_converter_
  void binSpectrumConverter(std::string);
  ///< Set Properties::binspectrum_converter_
# endif
};

#endif
