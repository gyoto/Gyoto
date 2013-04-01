/**
 * \file GyotoAstrobj.h
 * \brief Astronomical objects (light emitters)
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2011-2013 Thibaut Paumard, Frederic Vincent

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
      Subcontractor_t(Gyoto::FactoryMessenger*);
    ///< A function to build instances of a specific Astrobj::Generic sub-class
 
    /**
     * Instead of reimplementing the wheel, your subcontractor can simply be
     * Gyoto::Astrobj::Subcontractor<MyKind>
     */
    template<typename T> SmartPointer<Astrobj::Generic> Subcontractor
      (FactoryMessenger* fmp) {
      SmartPointer<T> ao = new T();
      ao -> setParameters(fmp);
      return ao;
    }
    ///< A template for Subcontractor_t functions

   /**
     * Query the Astrobj register to get the Astrobj::Subcontractor_t
     * correspondig to a given kind name. This function is normally
     * called only from the Factory.
     *
     * \param name e.g. "Star"
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Astrobj::Subcontractor_t* getSubcontractor(std::string name,
						      int errmode = 1);
    ///< Query the Astrobj register

#if defined GYOTO_USE_XERCES
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
     * communicate whith the Gyoto::Factory to build an instance of
     * the class from its XML description
     */
    void Register(std::string name, Gyoto::Astrobj::Subcontractor_t* scp);
    ///< Make an Astrobj kind known to the Factory
#endif
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
\code
Astrobj::Register("MyKind", &(Astrobj::Subcontractor<Astrobj::MyKind>));
\endcode
 */
/**
 * \class Gyoto::Astrobj::Generic
 * \brief Base class for astronomical object
 *
 * See introduction in the Gyoto::Astrobj namespace.
 */
class Gyoto::Astrobj::Generic : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Generic>;


  // Data : 
  // -----
 protected:

  /**
   * The Metric in this end of the Universe
   */
  SmartPointer<Gyoto::Metric::Generic> gg_;


  /**
   * Maximum distance from the center of the coordinate system at
   * which a photon may hit the object.  Child classes may choose to
   * update rmax at all time or to use it to cache the value, for
   * instance when getRmax() is called. External classes (Photons in
   * particular) must use getRmax() to access this information.
   *
   * rmax_set_==1 means that rmax_ was set using setRmax() or the
   * constructor. In this case, getRmax() must always return this
   * value, not recompute it.
   *
   * rmax_ is in geometrical units.
   */                                         
  double rmax_; ///< Maximum distance to the center of the coordinate system [geometrical units]

  /**
   * rmax_set_==1 means that rmax_ was set using setRmax() or the
   * constructor. In this case, getRmax() must always return this
   * value, not recompute it.
   *
   * Use unsetRmax() to reset rmax_set_ to 0.
   *
   */                                         
  int rmax_set_; ///< Never recompute rmax: it was externally set

  /**
   * The kind should match the name of the class, e.g. "Star" for a
   * Gyoto::Star.
   */
  const std::string kind_; ///< Kind of object (e.g. "Star"...)

  int flag_radtransf_; ///< 1 if radiative transfer inside Astrobj, else 0

  // Constructors - Destructor
  // -------------------------
 public:
  /**
   *  kind_ =  "Default", rmax_ = 0., rmax_set_ = 0.
   */
  Generic(); ///< Default constructor.

  /**
   *  kind_ =  "Default", rmax_ = radmax, rmax_set_ = 1.
   */
  Generic(double radmax); ///< Set rmax in constructor.

  /**
   *  kind_ =  kind, rmax_ = 0., rmax_set_ = 0.
   */
  Generic(std::string kind); ///< Set kind in constructor.

  /**
   * Make a deep copy of an Astrobj::Generic instance
   */
  Generic(const Generic& ) ; ///< Copy constructor.

  /**
   * This method must be implemented by the various Astrobj::Generic
   * subclasses in order to support cloning:
\code
SmartPointer<Astrobj> deep_copy = original->clone();
\endcode
   *
   * Implementing it is very straightforward, as long as the copy
   * constructor Generic(const Generic& ) has been implemented:
\code
MyAstrobj* MyAstrobj::clone() const { return new MyAstrobj(*this); }
\endcode
   */
  virtual Generic* clone() const = 0 ; ///< "Virtual" copy constructor
  
  virtual ~Generic() ; ///< Destructor: does nothing.

  // Accessors
  // ---------
 public:
  /**
   * Get the Metric
   */
  virtual SmartPointer<Metric::Generic> getMetric() const;

  /**
   * Set the Metric
   */
  virtual void setMetric(SmartPointer<Metric::Generic>) ;

  /**
   *  Get maximal distance from center of coordinate system at which a
   *  Photon may hit the object.
   *  
   *  Child classes may use the rmax_ member to cache this value.
   *
   *  It can also be set using setRmax(). If setRmax has been used
   *  to set rmax_, getRmax() must not recompute it.
   *
   *  \return rmax_ in geometrical units
   */
  virtual double getRmax(); ///< Get maximal distance from center of coordinate system

  /**
   *  Call setRmax() and convert result to unit.
   *
   *  \param string unit
   *  \return double rmax converted to unit
   */
  virtual double getRmax(std::string unit); ///< Get rmax_ is specified unit

  const std::string getKind() const; ///< Get the kind of the Astrobj (e.g. "Star")

  /**
   *  Set maximal distance from center of coordinate system at which a
   *  Photon may hit the object.
   *  
   *  Side effect: set rmax_set_ to 1.
   *  \param double val rmax_ in geometrical units.
   */
  virtual void setRmax(double val); ///< Set maximal distance from center of coordinate system

  /**
   *  Call Generic::setRmax(double val) after converting val from unit
   *  to geometrical units.
   *
   *  \param double val rmax_ expressed in unit "unit";
   *  \param string unit unit...
   */
  virtual void setRmax(double val, std::string unit); ///< Set maximal distance from center of coordinate system

  /**
   * getRmax() will then be free to recompute rmax_. Astrobjs
   * which cache rmax_ may need to update it when unsetRmax() is
   * called.
   */
  virtual void unsetRmax() ; ///< Set rmax_set_ to 0.

  /**
   * Set flag indicating that radiative transfer should be integrated,
   * i.e. the object is to be considered optically thin.
   * \param flag: 1 if optically thin, 0 if optically thick.
   */
  void setFlag_radtransf(int flag);
  ///< Set whether the object is optically thin.
  /**
   * See setFlag_radtransf(int flag).
   */
  int getFlag_radtransf() const ;
  ///< Query whether object is optically thin.

  /**
   * Return a Gyoto::Quantity_t suitable as input to
   * Gyoto::Scenery::setRequestedQuantities() to set de default
   * quantities to compute for this object. The default of these
   * defaults GYOTO_QUANTITY_INTENSITY.
   */
  virtual Quantity_t getDefaultQuantities();

  //XML I/O
 public:
  /**
   * \brief Set parameter by name
   *
   * Assume MyKind is a subclass of Astrobj::Generic which has two
   * members (a string StringMember and a double DoubleMember):
   * \code
   * int MyKind::setParameter(std::string name,
   *                          std::string content,
   *                          std::string unit) {
   *   if      (name=="StringMember") setStringMember(content);
   *   else if (name=="DoubleMember") setDoubleMember(atof(content.c_str()),
   *                                                  unit);
   *   else return Generic::setParameter(name, content, unit);
   *   return 0;
   * }
   * \endcode
   * If MyKind is not a direct subclass of Generic but is a subclass
   * of e.g. Standard, UniformSphere of ThinDisk, it should call the
   * corresponding setParameter() implementation instead of
   * Generic::setParameter().
   *
   * \param name XML name of the parameter
   * \param content string representation of the value
   * \param unit string representation of the unit
   * \return 0 if this parameter is known, 1 if it is not.
   */
  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit) ;

#ifdef GYOTO_USE_XERCES
  /**
   * Astrobj implementations should impement fillElement to save their
   * parameters to XML and call the generic implementation to save
   * generic parts such as Flag_radtrans: Generic::fillElement(fmp).
   */

  virtual void fillElement(FactoryMessenger *fmp) const ;
                                             ///< called from Factory

  /**
   * \brief Main loop in Subcontractor_t function
   *
   * The Subcontractor_t function for each Astrobj kind should look
   * somewhat like this (templated as
   * Gyoto::Astrobj::Subcontractor<MyKind>):
\code
SmartPointer<Astrobj::Generic>
Gyoto::Astrobj::MyKind::Subcontractor(FactoryMessenger* fmp) {
  SmartPointer<MyKind> ao = new MyKind();
  ao -> setParameters(fmp);
  return ao;
}
\endcode
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
   * adress. For instance, all objects know the "intensity"
   * property. If data->intensity != NULL, the instensity is computed
   * and stored in *data->intensity.
   *
   * If data is non-NULL and only in this case, processHitQuantities()
   * will also call ph->transmit() to update the transmissions of the
   * Photon (see Photon::transmit(size_t, double)). This must not be
   * done if data is NULL (see Astrobj::Complex::Impact() for an
   * explanation).
   *
   * \param ph   Gyoto::Photon aimed at the object;
   * \param index    Index of the last photon step;
   * \param data     Pointer to a structure to hold the observables at impact.
   *
   * \return 1 if impact, 0 if not.
   */
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL) = 0 ;
  ///< does a photon at these coordinates impact the object?
  
  /**
   * processHitQuantities fills the requested data in Impact. To use
   * it, you need to call it in the Impact() method for your object in
   * case of hit. It will fill Redshift, Intensity, Spectrum,
   * BinSpectrum and update the Photon's transmission by calling
   * Photon::transmi(), only if data==NULL.
   *
   * You can overload it for your Astrobj. The generic implementation
   * calls emission(), integrateEmission() and transmission() below.
   */
  virtual void processHitQuantities(Photon* ph, double* coord_ph_hit,
                                   double* coord_obj_hit, double dt,
                                   Astrobj::Properties* data) const;

  /**
   * Called by the default implementation for processHitQuantities().
   *
   * emission() computes the intensity I<SUB>&nu;</SUB> emitted by the small
   * volume of length dsem. It should take self-absorption along dsem
   * into account.
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
   *   d(I<SUB>&nu;</SUB>/&nu;<SUP>3</SUP>)/d&lambda; = (j<SUB>&nu;</SUB>/&nu;<SUP>2</SUP>)  [*]
   *  where &lambda; is the integration parameter along the null geodesic.
   *
   *      NB: Let us consider a particular observer, 
   *          with &nu; being the frequency measured by this observer,
   *          and ds being the proper distance (as measured by the observer) 
   *          that the photon travels as it moves
   *          from &lambda; to &lambda;+d&lambda; along its geodesic.
   *          Then it can be shown that :
   *                          d&lambda; = ds/&nu;
   *          This shows that Eq. [*] is homogeneous.
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
  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL)
    const ; ///< INVARIANT emission j<SUB>&nu;</SUB>/&nu;<SUP>2</SUP>

  /**
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
  virtual void emission(double Inu[], double nu_em[], size_t nbnu,
			double dsem, double coord_ph[8],
			double coord_obj[8]=NULL) const ; 

  /**
   * Compute the integral of emission() from &nu;<SUB>1</SUB> to &nu;<SUB>2</SUB>. The default
   * implementation is a numerical integrator which works well enough
   * and is reasonably fast if emission() is a smooth function
   * (i.e. no emission or absorption lines). If possible, it is wise
   * to implement an analytical solution. It is used by
   * processHitQuantities to compute the "BinSpectrum" quantity which
   * is the most physical: it is the only quantity that can be
   * actually measured directly by a real-life instrument.
   */
  virtual double integrateEmission(double nu1, double nu2, double dsem,
                                  double c_ph[8], double c_obj[8]=NULL) const;
    ///< &int;<SUB>nu<SUB>1</SUB></SUB><SUP>nu<SUB>2</SUB></SUP> I<SUB>&nu;</SUB> dnu (or j<SUB>&nu;</SUB>)
  /**
   * Like integrateEmission(double nu1, double nu2, double dsem,
   * double c_ph[8], double c_obj[8]) for each Spectrometer channel.
   */
  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, double *cph, double *co) const;


  /**
   * transmission() computes the transmission of this fluid element or
   * 0 if optically thick. The default implementation returns 1. (no
   * attenuation) if optically thin, 0. if optically thick.
   *
   * \param nuem frequency in the fluid's frame
   * \param coord Photon coordinate
   * \param dsem geometrical length in geometrical units
   */
  virtual double transmission(double nuem, double dsem, double coord[8]) const ;
     ///< Transmission: exp( &alpha;<SUB>&nu;</SUB> * ds<SUB>em</SUB> )

};

/**
 * \class Gyoto::Astrobj::Properties
 * \brief Observable properties of an Astronomical object
 *
 *  The sort of properties one wants to measure on a ray-traced Gyoto::Photon which hits a Gyoto::Astrobj. Not all Astrobj are able to fill all of these properties.
 *
 */
class Gyoto::Astrobj::Properties : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Properties>;
 public:
  double *intensity; ///< Apparent intensity (takes beaming into account); 
  double *time; ///< Date of impact (= date of emission of the photon);
  double *distance; ///< Behaves like the square of the closest distance between Photon and Astrobj (but not exactly that). Initialize it to DBL_MAX from float.h.;
  double * first_dmin; ///< first local minimum in distance from object
  int first_dmin_found; ///< first_dmin will be set to the first local minimum and first_dmin_found will be set to 1 if a local minimum in distance is found. Initialize it to 0.
  double *redshift; ///< redshift factor &nu;<SUB>obs</SUB>/&nu;<SUB>em</SUB> (necessary for emission lines computation)
  double *spectrum; ///< I<SUB>&nu;</SUB> (&nu;) (observed specific intensity)
  double *binspectrum; ///< I<SUB>&nu;<SUB>1</SUB></SUB><SUP>&nu;<SUB>2</SUB></SUP>, the integral of I<SUB>&nu;</SUB> over each spectral channel (i.e. what a spectrometer would measure)
  int offset; ///< spectra elements are separated by offset doubles in memory. In other words, the ith spectral element is a spectrum[i*offset].
  double * impactcoords; ///< Coordinates of the object and photon at impact
  double *user1, *user2, *user3, *user4, *user5; ///< Quantities specific to Astrobj
# ifdef HAVE_UDUNITS
  Gyoto::SmartPointer<Gyoto::Units::Converter> intensity_converter_ ;
  Gyoto::SmartPointer<Gyoto::Units::Converter> spectrum_converter_ ;
  Gyoto::SmartPointer<Gyoto::Units::Converter> binspectrum_converter_ ;
# endif
 public:
  Properties(); ///< Default constructor (everything is set to NULL);
  Properties (double*, double*); ///<< set intensity and time pointers.
  void init(size_t nbnuobs=0);
  Properties operator++();
# ifdef HAVE_UDUNITS
  void setIntensityConverter(Gyoto::SmartPointer<Gyoto::Units::Converter>);
  void setIntensityConverter(std::string);
  void setSpectrumConverter(Gyoto::SmartPointer<Gyoto::Units::Converter>);
  void setSpectrumConverter(std::string);
  void setBinSpectrumConverter(Gyoto::SmartPointer<Gyoto::Units::Converter>);
  void setBinSpectrumConverter(std::string);
# endif
};

#endif
