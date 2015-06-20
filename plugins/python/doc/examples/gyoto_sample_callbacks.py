'''Sample spectra for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Spectrum implementations using Gyoto's "python" plug-in. Note that
   this plug-in can be renamed to whatever matches the particular
   version of Python it has been built against (e.g. python3.4).

   The goal is to be able to instanciate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto
   gyoto.loadPlugin("python") # or python2.7 or python3.4...
   sp=gyoto.Spectrum("Python")
   sp.set("Module", "gyoto_sample_callbacks")
   sp.set("Class", "PowerLaw") # or "BlackBody6000"
   sp.set("Parameters", (1., 2.))
   val=sp(3e8/2e-6)

   Spectra:

   Classes that aim at implementing the Gyoto::Spectrum::Generic
   interface do so by providing the following methods:

   __call__(self, nu): mandatory;
   __setitem__: optional;
   integrate: optional.

   Metrics:

   Classes that aim at implementing the Gyoto::Metric::Generic
   interface do so by providing the following methods:

   gmunu(self, dst, pos): mandatory;
   christoffel(self, dst, pos): mandatory
   __setitem__(self, key, value): mandatory

   Astrobjs:

   Classes that aim at implementing the Gyoto::Metric::Generic
   interface do so by providing the following methods:

   __call__: required
   getVelocity: required
   giveDelta, emission, integrateEmission, transmission, __setitem__:
              optional.
   emission and integrateEmission can be overloaded by using the
   varargs argument.

'''

import math
import numpy
import gyoto

class BlackBody6000:
    '''Black-body spectrum at 6000K

    Parameters: none.

    This example is pretty minimal: it's a black-body with fixed
    temperature. We implement the only mandaroty function:
    __call__(self, nu).
    '''

    def __call__(self, nu):
        '''spectrum(frequency_in_Hz) = black-body distribution for T=6000K.

        This function implements only
        Gyoto::Spectrum::Python::operator()(double nu).

        It does so by not accepting the varargs argument. Contrast
        with the PowerLaw definition of __call__.

        '''
        temperature=6000.
        PLANCK_OVER_BOLTZMANN=4.7992373e-11
        return nu*nu*nu/(math.exp(PLANCK_OVER_BOLTZMANN*nu/temperature)-1.);


class PowerLaw:
    '''Powerlaw spectrum

    Parameters: (constant, exponent)

    This example is pretty complete. It implements everything usefull
    and some eye-candy.

    '''

    constant=0.
    exponent=0.

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


class Minkowski:
    '''Flat space metric

    Implemented for both Cartesian and spherical coordinates.

    Every Gyoto Metric implemented in Python must implement the three
    methods illustrated here.

    '''

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
            pass
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
        spherical=self.this.get('Spherical')
        for mu in range(0, 4):
            for nu in range(0, 4):
                g[mu][nu]=g[nu][mu]=0
        g[0][0]=-1;
        if not spherical:
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
        spherical=self.this.get('Spherical')
        for alpha in range(0, 4):
            for mu in range(0, 4):
                for nu in range(0, 4):
                    dst[alpha][mu][nu]=0.
        if not spherical:
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


class FixedStar:
    ''' Sample class for Astrobj::Python::Standard
    '''
    def __init__(self):
        '''Initialize instance

        Needed here to make a non-static array data member.
        '''
        self.pos = numpy.zeros((4), float)

    def __setitem__(self, key, value):
        '''Set parameters

        Here, the parameters will be the 3 space coordinates of the
        center of the blob.

        '''
        if key in (0, 1, 2):
            self.pos[key+1]=value
        else:
            raise IndexError
        self.coord_st=self.to_cartesian(self.pos)

    def to_cartesian(self, coord):
        '''Helper function, not in the API

        '''
        gg=self.this.metric()
        spherical=False
        if gg is not None:
            spherical = gg.coordKind() == gyoto.GYOTO_COORDKIND_SPHERICAL
        if spherical:
            rs=coord[1]
            ths=coord[2]
            phs=coord[3]
            st=math.sin(ths)
            ct=math.cos(ths)
            sp=math.sin(phs)
            cp=math.cos(phs)
            return numpy.array((coord[0], rs*st*cp, rs*st*sp, rs*ct))
        return coord

    def __call__(self, coord):
        ''' Astrobj::Standard::operator()()

        Required
        '''
        coord_ph=self.to_cartesian(coord)
        coord_st=self.coord_st
        dx = coord_ph[1]-coord_st[1]
        dy = coord_ph[2]-coord_st[2]
        dz = coord_ph[3]-coord_st[3]
        return math.sqrt(dx*dx + dy*dy + dz*dz)

    def getVelocity(self, coord, vel):
        ''' Velocity field

        Required
        '''
        vel[0]=1.
        for i in range(1, 4):
            vel[i]=0.

    def emission(self, nuem, dsem, cph, co):
        ''' emission

        Optional
        '''
        return 1.

class ThinDisk:
    '''A ThinDisk with funny emission
    '''
    def emission(self, nuem, dsem, cph, co):
        return dsem
