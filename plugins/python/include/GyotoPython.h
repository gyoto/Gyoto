/*
    Copyright 2015-2016 Thibaut Paumard

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

/**
 * \file GyotoPython.h
 * \brief Extending Gyoto using Python
 *
 * The classes provided here allow implementing a Spectrum, an Astrobj
 * or a Metric in Python. Together, they form the "python" Gyoto
 * plug-in.
 *
 * This is complementary to, but distinct from the "gyoto" Python
 * extension. Here, we are embedding Python inside Gyoto so that a few
 * functions are coded in Python. The Python extension does the
 * reverse: it allows calling Gyoto functions from within
 * Python. Beware of the vocabulary: here, we call "plug-in" a shared
 * library that extends Gyoto, and "extension" a shared library that
 * extends Python.
 *
 * The plug-in works within the gyoto command-line utility as well as
 * when Gyoto is used inside Python or inside Yorick. The only caveat
 * is that the python plug-in of Gyoto should not be loaded into a
 * Python interpreter different from the one that was used for
 * building the plug-in.
 *
 * For this reason, the name of this plug-in depends on the Python
 * interpreter that was used when building. It can be simply "python",
 * or a be versionned: for instance "python2.7" or "python3.4". This
 * way, it is possible to keep several copies of the plug-in, one for
 * each version of the Python interpreter that are installed on the
 * machine. Any version can be used in the gyoto command-line utility
 * or in Yorick, but when Gyoto is used inside Python, only the
 * matching version of this plug-in may be used.
 *
 * Implementing a Spectrum, Astrobj or Metric kind in Python is much
 * easier than implementing a new C++ plug-in for Gyoto. This saves in
 * development time. However, there is a cost in terms of computing
 * time. While this cost may not be noticeable for Spectra and is
 * moderate for Astrobjs (at least for simple ones), it is quite
 * significant for Metrics, because the gmunu and christoffel methods
 * are evaluated several times per integration step, for every
 * photon. On one example using the Minkowski Metric, the integration
 * of a full image with the Python implementation took approx. 150-200
 * more time than the same integration with the C++
 * implementation. So, for Metrics, the Python implementation can
 * serve as a prototyping test-bed, but most users will probably still
 * want to re-implement their Metrics in C++ eventually.
 *
 * Note also that multi-threading is not very efficient for the
 * Metric::Python class, because only one thread may interact with the
 * Python interpreter at any time. MPI multi-processing runs much
 * faster. Here again, this limitation is less problematic for Spectra
 * and Astrobjs than it is for Metrics.
 *
 * The principle of these classes is very simple: the plugin embeds a
 * Python interpreter. Each instance of the Gyoto classes
 * Gyoto::Metric::Python, Gyoto::Spectrum::Python,
 * Gyoto::Astrobj:Python::Standard and
 * Gyoto::Astrobj::Python::ThinDisk instantiate a Python class in this
 * interpreter, and delegate certain methods from the Gyoto API to
 * this instance.
 *
 * In simple cases, the Python instance does not even need to know
 * that it is running in Gyoto. It simply exposes an interface that is
 * being called. However, Gyoto sets a few attributes in each
 * method. Most notably, if the "gyoto" python extension is available,
 * Gyoto will set the attribute "this" to the C++ instance that
 * created the Python class instance, so that the Python code can
 * access C++-side information.
 *
 */

#ifndef __GyotoPython_H_ 
#define __GyotoPython_H_ 
#include <GyotoSpectrum.h>
#include <GyotoMetric.h>
#include <GyotoStandardAstrobj.h>
#include <GyotoThinDisk.h>
#include <Python.h>

namespace Gyoto {
  /**
   * \namespace Gyoto::Python
   * \brief Helpers for the classes deriving from Gyoto::Python::Base
   */

  namespace Python {
    class Base;

    /// Return new reference to method, or NULL if method not found.
    PyObject * PyInstance_GetMethod(PyObject* pInstance, const char *name);

    /// Return refernce to the gyoto module, or NULL.
    PyObject * PyImport_Gyoto();

    /// Set "this" attribute in instance
    void PyInstance_SetThis(PyObject * pInstance,
			    PyObject * pNew,
			    void * ptr);

    /// Check whether method accepts the varargs argument
    bool PyCallable_HasVarArg(PyObject * pMethod);

    /// Create module from Python source code in a C string
    PyObject * PyModule_NewFromPythonCode(const char * code);

    /// Get reference to the Spectrum constructor in the gyoto Python extension
    PyObject * pGyotoSpectrum() ;
    /// Get reference to the Metric constructor in the gyoto Python extension
    PyObject * pGyotoMetric() ;
    /// Get reference to the StandardAstrobj constructor in the gyoto Python extension
    PyObject * pGyotoStandardAstrobj() ;
    /// Get reference to the ThinDisk constructor in the gyoto Python extension
    PyObject * pGyotoThinDisk() ;
  }
  namespace Spectrum {
    class Python;
  }
  namespace Metric {
    class Python;
  }
  namespace Astrobj {
  /**
   * \namespace Gyoto::Astrobj::Python
   * \brief Classes that wrap Python classes as Gyoto Astrobj implementations
   */
    namespace Python {
      class Standard;
      class ThinDisk;
    }
  }
}

/**
 * \class Gyoto::Python::Base
 *
 * \brief Base class for classes in the Python plug-in.
 *
 * All classes have those three Properties:
 * - Module (string): the module in which the Python class is
 *   implemented;
 * - Class (string): the name of the Python class, in module Module,
 *   to interface with;
 * - Parameters (vector<double>): list of parameters for this
 *   class. These parameters are passed one by one to the Python
 *   instance using __setitem__ with numerical keys.
 *
 * All the Gyoto instances of the classes descending from
 * Gyoto::Python::Base expose themselves to the Python instance they
 * wrap immediately after instantiation by setting the 'this'
 * attribute. If the 'gyoto' Python extension can be loaded, then
 * 'this' will be an instance of one of the classes gyoto.Metric,
 * gyoto.Spectrum, gyoto.StandardAstrobj or gyoto.ThinDisk pointing to
 * the underlying C++ instance. If the 'gyoto' extension is not
 * available, 'this' will be None.
 *

 */
class Gyoto::Python::Base {
 protected:
  /**
   * \brief Name of the Python module that holds the class
   *
   * For instance, if the class is implemented in toto.py, the module
   * name is "toto". Property name: Module.
   */
  std::string module_;

  /**
   * \brief Python source code for module that holds the class
   *
   */
  std::string inline_module_;

  /**
   * \brief Name of the Python class that we want to expose
   *
   * Property name: Class.
   */
  std::string class_;

  /**
   * \brief Parameters that this class needs
   *
   * A list of parameters (doubles) can be passed in the Property
   * Parameters. They will be sent to the Python instance using
   * __setitem__.
   */
  std::vector<double> parameters_;

  /**
   * \brief Reference to the python module once it has been loaded.
   */
  PyObject * pModule_;

  /**
   * \brief Reference to the python instance once it has been instantiated.
   */
  PyObject * pInstance_;

 public:
  Base();
  Base(const Base&);
  ~Base();

  virtual std::string module() const ; ///< Return module_
  virtual std::string inlineModule() const ; ///< Return inline_module_

  /**
   * \brief Set #module_ and import the Python module
   *
   * Side effects:
   *   - sets #inline_module_ to "";
   *   - calls #klass(#class_) if #class_ is already known, so #module_
   *     can be set (or reset) after #class_.
   */
  virtual void module(const std::string&);

  /**
   * \brief Set #inline_module_ and import the Python module
   *
   * Side effects:
   *   - sets #module_ to "";
   *   - calls #klass(#class_) if #class_ is already known, so #module_
   *     can be set (or reset) after #class_.
   */
  virtual void inlineModule(const std::string&);

  /// Retrieve #class_.
  virtual std::string klass() const ;

  /**
   * \brief Set #class_ and instantiate the Python class.
   *
   * Sets #pInstance_.
   *
   * This generic implementation takes care of the common ground, but
   * does not set 'this' or call #parameters(#parameters_). Therefore,
   * all the derived classes should reimplement this method and at
   * least call Python::Base::klass(c) and
   * #parameters(#parameters_). Between the two is the right moment to
   * check that the Python class implements the required API and to
   * cache PyObject* pointers to class methods.
   */
  virtual void klass(const std::string& c);

  /// Retrieve #parameters_
  virtual std::vector<double> parameters() const;
  /**
   * \brief Set #parameters_ and send them to #pInstance_
   *
   * The parameters are sent to the class instance using the
   * __setitem__ method with numerical keys.
   */
  virtual void parameters(const std::vector<double>&);

};



/**
 * \class Gyoto::Spectrum::Python
 *
 * \brief Loader for Python classes implementing the Spectrum interface
 *
 * It interfaces with a Python class which must implement at least the
 * __call__ method.
 *
 * Sample XML file:
 * \include example-python.xml
 * Sample Python module:
 * \include gyoto_sample_spectra.py
 */
class Gyoto::Spectrum::Python
: public Gyoto::Spectrum::Generic,
  public Gyoto::Python::Base
{
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::Python>;
 protected:

  /**
   * \brief Reference to ___call__
   *
   * __call__ is the method in the underlying Python class that
   * implements Gyoto::Spectrum::Generic::operator()().
   */
  PyObject * pCall_;

  /**
   * \brief Reference to the (optional) integrate method
   */
  PyObject * pIntegrate_;

  /**
   * \brief Whether __call__ is overloaded
   *
   * This is determined automatically by looking at the parameters
   * accepted by __call__:
   * \code
   *   def __call__(self, nu)
   * \endcode
   * In this case call is not overloaded and implements only virtual
   * double operator()(double nu) const;
   * \code
   *   def __call__(self, *args)
   * \endcode
   * In this case __call__ is overloaded and implements both double
   * operator()(double nu) const and virtual double operator()(double
   * nu, double opacity, double ds) const.
   */
  bool pCall_overloaded_;

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

  Python();

  Python(const Python&);

  virtual Python * clone() const;

  ~Python();

  // For some reason we need to implement the bunch although only one
  // is non-trivial
  virtual std::string module() const ;
  virtual void module(const std::string&);
  virtual std::string inlineModule() const ;
  virtual void inlineModule(const std::string&);
  virtual std::string klass() const ;
  virtual void klass(const std::string&);
  virtual std::vector<double> parameters() const;
  virtual void parameters(const std::vector<double>&);

  virtual double operator()(double nu) const;
  virtual double operator()(double nu, double opacity, double ds) const;

  virtual double integrate(double nu1, double nu2) ;

};


/**
 * \class Gyoto::Metric::Python
 * \brief Metric coded in Python
 *
 * Loader for Python Metric classes. It interfaces with a Python class
 * which must implement at least the gmunu and christoffel methods.
 *
 * Other methods are optional: getRmb, getRms, getSpecificAngularMomentum,
 * getPotential.
 *
 * Use &lt;Cartesian&gt; or &lt;/Spherical&gt; to select the coordinate system
 * kind.
 *
 * Sample XML file:
 * \include example-python.xml
 * Sample Python module:
 * \include gyoto_sample_metrics.py
 */
class Gyoto::Metric::Python
: public Gyoto::Metric::Generic,
  public Gyoto::Python::Base
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::Python>;

 private:
  // Variables to cache Python objects:
  /**
   * \brief Reference to the gmunu method
   */
  PyObject * pGmunu_;

  /**
   * \brief Reference to the christoffel method
   */
  PyObject * pChristoffel_;

  /**
   * \brief Reference to the getRmb method
   */
  PyObject * pGetRmb_;

  /**
   * \brief Reference to the getRms method
   */
  PyObject * pGetRms_;

  /**
   * \brief Reference to the  method getSpecificAngularMomentum
   */
  PyObject * pGetSpecificAngularMomentum_;

  /**
   * \brief Reference to the getPotential method
   */
  PyObject * pGetPotential_;

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  Python();
  Python(const Python&);
  ~Python();
  virtual Python* clone() const ;

  // Accessors for the Gyoto::Property members:
  // Those are mere wrappers arround Generic::coordKind(), useful for
  // declaring a boolen property using the macro GYOTO_PROPERTY_BOOL:
  void spherical(bool);
  bool spherical() const;
  virtual std::string module() const ;
  virtual void module(const std::string&);
  virtual std::string inlineModule() const ;
  virtual void inlineModule(const std::string&);
  virtual std::string klass() const ;
  virtual void klass(const std::string&);
  virtual std::vector<double> parameters() const;
  virtual void parameters(const std::vector<double>&);
  using Gyoto::Metric::Generic::mass;
  virtual void mass(double m);

  // The minimal Gyoto::Metric API:
  void gmunu(double g[4][4], const double * x) const ;
  int christoffel(double dst[4][4][4], const double * x) const ;

  // Little more
  double getRmb() const;
  double getRms() const;
  double getSpecificAngularMomentum(double rr) const;
  double getPotential(double const pos[4], double l_cst) const;

};

/**
 * \class Gyoto::Astrobj::Python::Standard
 * \brief Coding a Gyoto::Astrobj::Standard in Python
 *
 * Sample XML file:
 * \include example-python-standard.xml
 * Sample Python module:
 * \include gyoto_sample_standard.py
 */
class Gyoto::Astrobj::Python::Standard
: public Gyoto::Astrobj::Standard,
  public Gyoto::Python::Base
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Python::Standard>;

 private:
  PyObject *pEmission_, *pIntegrateEmission_, *pTransmission_, *pCall_,
    *pGetVelocity_, *pGiveDelta_;
  bool pEmission_overloaded_, pIntegrateEmission_overloaded_;

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

  /* Birth and Death*/
  Standard();
  Standard(const Standard&);
  ~Standard();
  Standard* clone() const;

  /* Astrobj::Generic API */
  virtual double emission(double nu_em, double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;

  virtual void emission(double Inu[], double const nu_em[], size_t nbnu,
			double dsem, state_t const &coord_ph,
			double const coord_obj[8]=NULL) const ;

  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   state_t const &c_ph, double const c_obj[8]=NULL) const;

  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, state_t const &cph, double const *co) const;

  virtual double transmission(double nuem, double dsem, state_t const &cph, double const *co) const ;

  /* Astrobj::Standard API */
  virtual double operator()(double const coord[4]) ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  virtual double giveDelta(double coord[8]);

  /* Python::Base */
  virtual std::string module() const ;
  virtual void module(const std::string&);
  virtual std::string inlineModule() const ;
  virtual void inlineModule(const std::string&);
  virtual std::string klass() const ;
  virtual void klass(const std::string&);
  virtual std::vector<double> parameters() const;
  virtual void parameters(const std::vector<double>&);
  virtual double criticalValue() const ;
  virtual void criticalValue(double) ;

};

/**
 * \class Gyoto::Astrobj::Python::ThinDisk
 * \brief Coding a Gyoto::Astrobj::ThinDisk in Python
 *
 * Sample XML file:
 * \include example-python-thindisk.xml
 * Sample Python module:
 * \include gyoto_sample_thindisks.py
 */
class Gyoto::Astrobj::Python::ThinDisk
: public Gyoto::Astrobj::ThinDisk,
  public Gyoto::Python::Base
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Python::ThinDisk>;

 private:
  PyObject *pEmission_, *pIntegrateEmission_, *pTransmission_, *pCall_,
    *pGetVelocity_, *pGiveDelta_;
  bool pEmission_overloaded_, pIntegrateEmission_overloaded_;

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

  /* Birth and Death*/
  ThinDisk();
  ThinDisk(const ThinDisk&);
  ~ThinDisk();
  ThinDisk* clone() const;

  /* Astrobj::Generic API */
  virtual double emission(double nu_em, double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;

  virtual void emission(double Inu[], double const nu_em[], size_t nbnu,
			double dsem, state_t const &coord_ph,
			double const coord_obj[8]=NULL) const ;

  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   state_t const &c_ph, double const c_obj[8]=NULL) const;

  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, state_t const &cph, double const *co) const;

  virtual double transmission(double nuem, double dsem, state_t const &cph ,double const *co) const ;

  /* Astrobj::ThinDisk API */
  virtual double operator()(double const coord[4]) ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;

  /* Python::Base */
  virtual std::string module() const ;
  virtual void module(const std::string&);
  virtual std::string inlineModule() const ;
  virtual void inlineModule(const std::string&);
  virtual std::string klass() const ;
  virtual void klass(const std::string&);
  virtual std::vector<double> parameters() const;
  virtual void parameters(const std::vector<double>&);

};


#endif
