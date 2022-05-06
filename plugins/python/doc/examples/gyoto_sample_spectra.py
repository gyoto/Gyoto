'''Sample spectra for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Spectrum implementations using Gyoto's "python" plug-in. Note that
   this plug-in can be renamed to whatever matches the particular
   version of Python it has been built against (e.g. python3.4).

   The goal is to be able to instantiate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto.core
   gyoto.core.requirePlugin("python") # or python2.7 or python3.4...
   sp=gyoto.core.Spectrum("Python")
   sp.set("Module", "gyoto_sample_spectra")
   sp.set("Class", "PowerLaw") # or "BlackBody6000"
   sp.set("Parameters", (1., 2.))
   val=sp(3e8/2e-6)

   Classes that aim at implementing the Gyoto::Spectrum::Generic
   interface do so by providing the following members and methods:

   __call__(self, nu): mandatory;
   __setitem__: optional;
   integrate: optional;
   properties, set, get: optional.

'''

import math
import numpy
import gyoto.core

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

    properties={"Constant": "double", "Exponent": "double"}

    def __setitem__(self, key, value):
        '''
        This is how Gyoto sends the <Parameters/> XML entity:
        spectrum[i]=value
        i=0: set constant
        i=1: set exponent
        '''
        if (key==0 or key == "Constant"):
            self.constant = value
        elif (key==1 or key == "Exponent"):
            self.exponent = value
        else:
            raise IndexError
    set=__setitem__

    def get(self, key):
        '''
        Implementing this is absolutely not necessary (Gyoto does not
        use it, as of now), but we can: it allows retrieving the
        parameters like __setitem__ sets them:

        spectrum[0] == spectrum.constant
        spectrum[1] == spectrum.exponent
        '''
        if (key==0 or key=="Constant"):
            return self.constant
        elif (key==1 or key=="Exponent"):
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
