/*
    Copyright 20015 Thibaut Paumard

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
 * The classes provided here allow implementing a Spectrum or a Metric
 * in Python. Together, they form the a Gyoto plug-in.
 *
 * This is complementary to, but distinct from the Gyoto * Python
 * extension. Here, we are embedding Python inside Gyoto so that a few
 * functions are coded in Python. The Python extension does the
 * reverse: it allows calling Gyoto functions from within
 * Python. Beware of the vocabulary: here, we call "plug-in" a shared
 * library that extends Gyoto, and "extension" a shared library that
 * extends Python.
 *
 * The plug-in works within the gyoto command-line utility as well as
 * when Gyoto is used inside Python or inside Yorick. The only caveat
 * is that the Python plug-in of Gyoto should not be loaded into a
 * Python interpreter different from the one that was used for
 * building the plug-in.
 *
 * For this reason, the name of this plug-in depends on the Python
 * interpreter that was used when building. It can be simply "python",
 * or a be versionned: for instance "python2.7" or "python3.4". This
 * way, it is possible to keep several copies of the plug-in. Any
 * version can be used in the gyoto command-line utility or in Yorick,
 * but when Gyoto is used inside Python, only the matching version of
 * this plug-in may be used.
 *
 * Imlpementing a Spectrum or Metric kind in Python is much easier
 * than implementing a new C++ plug-in for Gyoto. This saves in
 * development time. However, there is a cost in terms of computing
 * time. While this cost may not be noticeable for Spectra, it is
 * quite significant for Metrics. On one example using the Minkowski
 * Metric, the integration of a full image with the Python
 * implementation took approx. 150-200 more time than the same
 * integration with the C++ implementation. So the Python
 * implementation can serve as a prototyping test-bed, but most users
 * will probably still want to re-implement their Metrics in C++
 * eventually.
 *
 * Note also that multi-threading is not very efficient for the
 * Metric::Python class, because only one thread may interact with the
 * Python interpreter at any time. MPI multi-processing runs much
 * faster.
 *
 */

#ifndef __GyotoPython_H_ 
#define __GyotoPython_H_ 
#include <GyotoSpectrum.h>
#include <GyotoMetric.h>
#include <GyotoStandardAstrobj.h>
#include <Python.h>

namespace Gyoto {
  namespace Python {
    class Base;
    PyObject * PyInstance_GetMethod(PyObject* pInstance, const char *name);
    bool PyCallable_HasVarArg(PyObject * pMethod);
  }
  namespace Spectrum {
    class Python;
  }
  namespace Metric {
    class Python;
  }
  namespace Astrobj {
    namespace Python {
      class Standard;
    }
  }
}

/**
 * \class Gyoto::Python::Base
 *
 * \brief Base class for classes in the Python plug-in.
 *
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
   * \brief Name of the Python class that we want to expose
   *
   * Property name: Class.
   */
  std::string class_;

  /**
   * \brief Parameters that this class needs
   *
   * A list of parameters (doubles) can be passed in the Property
   * Parameters.
   */
  std::vector<double> parameters_;

  /**
   * \brief Reference to the python module once it has been loaded.
   */
  PyObject * pModule_;

  /**
   * \brief Reference to the python instance once it has been instanciated.
   */
  PyObject * pInstance_;

 public:
  Base();
  Base(const Base&);
  ~Base();

  virtual std::string module() const ; ///< Return module_

  /**
   * \brief Set #module_ and import the Python module
   *
   * Also calls #klass(#class_) if #class_ is already known, so #module_
   * can be set (or reset) after #class_.
   */
  virtual void module(const std::string&);

  /// Retrieve #class_.
  virtual std::string klass() const ;

  /**
   * \brief Set #class_ and instanciate the Python class.
   *
   * Sets #pInstance_.
   *
   * This generic implementation takes care of the common ground, but
   * does not call #parameters(#parameters_). Therefore, all the derived
   * classes should reimplement this method and at least call
   * Python::Base::klass(c) and #parameters(#parameters_). Between the
   * two is the right moment to check that the Python class implements
   * the required API and to cache PyObject* pointers to class methods.
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
 *  XML stanza:
 *  \code
 *    <Spectrum kind="Python">
 *      <Module>my_python_module</Module>
 *      <Class>my_python_class</Class>
 *      <Parameters> 0. 1. 2. ... </Parameters>
 *    </Spectrum>
 *  \endcode
 *
 * Example:
 *
\code
class PowerLaw:
    def __setitem__(self, key, value):
        '''
        This is how Gyoto sends the <Parameter/> XML entity:
        spectrum[i]=value
        i=0: set constant
        i=1: set exponent
        '''
        if (key==0):
            self.constant = value
        elif (key==1):
            self.exponent = value
        else:
            raise IndexError

    def __getitem__(self, key, value):
        '''
        Implementing this is absolutely not necessary (Gyoto does not
        use it, as of now), but we can: it allows retrieving the
        parameters like __setitem__ sets them:

        spectrum[0] == spectrum.constant
        spectrum[1] == spectrum.exponent
        '''
        if (key==0):
            return self.constant
        elif (key==1):
            return self.exponent
        else:
            raise IndexError

    def __call__(self, *args):
        '''spectrum(frequency_in_Hz) = constant * nu**exponent

        This function implements both
        Spectrum::Python::operator()(double nu).
        and
        Spectrum::Python::operator()(double nu, double opacity, double ds).

        This behavior is obtained by having the varargs *args as
        second argument instead of a normal variable.

	To implement only Spectrum::Python::operator()(double nu), use
	this definition instead: def __call__(self, nu):

        The second overloaded function is here exactly the same as the
        C++ generic implementation and therefore useless. It is here
        to illustrate the API.

        '''
        nu=args[0]
        if (len(args)==1):
            return self.constant * math.pow(nu, self.exponent)
        else:
            opacity=args[1]
            ds=args[2]
            thickness=(opacity*ds)
            if (thickness):
                return self(nu) * (1.-math.exp(-thickness))
            return 0.

    def integrate(self, nu1, nu2):
        '''
        If present, this function implements
        Gyoto::Spectrum::Python::integrate(double nu1, double nu2)

        If absent, the generic integrator is used.
        '''
        if (self.exponent == -1.):
            return self.constant * (math.log(nu2) -math.log(nu1))
        return self.constant * (math.pow(nu2, self.exponent+1)- math.pow(nu1, self.exponent+1)) / (self.exponent+1)
\endcode
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

  Python();

  Python(const Python&);

  virtual Python * clone() const;

  ~Python();

  // For some reason we need to implement the bunch although only one
  // is non-trivial
  virtual std::string module() const ;
  virtual void module(const std::string&);
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
 * which must implement at least the methods detailed below.
 *
 * Use &lt;Cartesian&gt; or &lt;/Spherical&gt; to select the coordinate system
 * kind.
 *
 *  XML stanza:
 *  \code
 *    <Metric kind="Python">
 *      <Module>my_python_module</Module>
 *      <Class>my_python_class</Class>
 *      <Spherical/> or Cartesian
 *      <Mass unit="kg"> 1. </Mass>
 *      <Parameters> 0. 1. 2. ... </Parameters>
 *    </Spectrum>
 *  \endcode
 *
 * Example:
 *
\code
class Minkowski:
    '''Flat space metric

    Implemented for both Cartesian and spherical coordinates.

    Every Gyoto Metric implemented in Python must implement the three
    methods illustrated here.

    '''
    spherical = False

    def __setitem__(self, key, value):
        '''Set parameters.

        Mandatory.

        At least 'spherical' and 'mass' must be supported. If only one
        kind (Spherical or Cartesian) is supported, __setitem__ must
        still accept both True or False as valid, it is not the right
        place to raise an exception. Do this in gmunu or christoffel
        or both.

        Additional parameters, if any, will be sent using integer keys
        like in the SPectrum examples.

        '''
        if key == "spherical":
            spherical = value
        elif key == "mass":
            # C++ may send a mass, we accept it but ignore it.
            pass
        else:
            raise IndexError

    def gmunu(self, g, x):
        ''' Gyoto::Metric::Generic::gmunu(double dst[4][4], const double pos[4])

        Mandatory.

        C++ will send two NumPy arrays.

        '''
        for mu in range(0, 4):
            for nu in range(0, 4):
                g[mu][nu]=g[nu][mu]=0
        g[0][0]=-1;
        if not self.spherical:
            for mu in range(1, 4):
                g[mu][mu]=1.
            return
        r=x[1]
        theta=x[2]
        tmp=r*math.sin(theta)
        g[1][1]=1.
        g[2][2]=r*r
        g[3][3]=tmp*tmp

    def christoffel(self, dst, x):
        '''Gyoto::Metric::Generic::christoffel(double dst[4][4][4], const double pos[4])

        Mandatory.

        C++ will send two NumPy arrays.

        '''
        for alpha in range(0, 4):
            for mu in range(0, 4):
                for nu in range(0, 4):
                    dst[alpha][mu][nu]=0.
        if not self.spherical:
            return 0
        r=x[1]
        theta=x[2]
        sth=math.sin(theta)
        cth=math.cos(theta)
        dst[1][2][2]=-r
        dst[1][3][3]=-r*sth*sth
        dst[2][1][2]=dst[2][2][1]= 1./r
        dst[2][3][3]=-sth*cth
        dst[3][1][3]=dst[3][3][1]= dst[2][1][2]
        dst[3][2][3]=dst[3][3][2]= math.tan(math.pi*0.5 - x[2])
        return 0
\endcode
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

 public:
  GYOTO_OBJECT;
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
  virtual std::string klass() const ;
  virtual void klass(const std::string&);
  virtual std::vector<double> parameters() const;
  virtual void parameters(const std::vector<double>&);
  using Gyoto::Metric::Generic::mass;
  virtual void mass(double m);

  // The minimal Gyoto::Metric API:
  void gmunu(double g[4][4], const double * x) const ;
  int christoffel(double dst[4][4][4], const double * x) const ;

};

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

  /* Birth and Death*/
  Standard();
  Standard(const Standard&);
  ~Standard();
  Standard* clone() const;

  /* Astrobj::Generic API */
  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL) const ;

  virtual void emission(double Inu[], double nu_em[], size_t nbnu,
			double dsem, double coord_ph[8],
			double coord_obj[8]=NULL) const ;

  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const;

  virtual void integrateEmission(double * I, double const * boundaries,
				 size_t const * chaninds, size_t nbnu,
				 double dsem, double *cph, double *co) const;

  virtual double transmission(double nuem, double dsem, double coord[8]) const ;

  /* Astrobj::Standard API */
  virtual double operator()(double const coord[4]) ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  virtual double giveDelta(double coord[8]);

  /* Python::Base */
  virtual std::string module() const ;
  virtual void module(const std::string&);
  virtual std::string klass() const ;
  virtual void klass(const std::string&);
  virtual std::vector<double> parameters() const;
  virtual void parameters(const std::vector<double>&);
  virtual double criticalValue() const ;
  virtual void criticalValue(double) ;

};


#endif
