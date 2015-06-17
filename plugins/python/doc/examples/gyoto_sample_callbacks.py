'''
   Sample spectra for using with Gyoto Python plug-in.

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
   sp.set("Class", "PowerLaw") # or "BlackBody"
   sp.set("Parameters", (1., 2.))
   val=sp(3e8/2e-6)
'''

import math
class BlackBody:
    '''
    Black-body spectrum

    Parameters: temperature(K), scaling.
    '''

    temperature=1.
    scaling=1.
    PLANCK_OVER_BOLTZMANN=4.7992373e-11
    
    def setParameters(self, *args):
        '''
        spectrum.setParameters([temperature[, scaling]])
        '''
        if (len(args)>0):
            self.temperature = args[0]
            if (len(args)>1):
                self.scaling = args[1]
    
    def __call__(self, nu):
        '''
        spectrum(frequency_in_Hz) = black-body distribution
        '''
        return self.scaling*nu*nu*nu/(math.exp(self.PLANCK_OVER_BOLTZMANN*nu/self.temperature)-1.);

class PowerLaw:
    '''
    Powerlaw spectrum

    Parameters: (constant, exponent)
    '''

    constant=0.
    exponent=0.

    def setParameters(self, *args):
        '''
        spectrum.setparameters([constant[, exponent]])
        '''
        if (len(args)>=1):
            self.constant = args[0]
            if (len(args)>=2):
                self.exponent = args[1]

    def __call__(self, nu):
        '''
        spectrum(frequency_in_Hz) = constant * nu**exponent
        '''
        return self.constant * math.pow(nu, self.exponent)

    def integrate(self, nu1, nu2):
        if (self.exponent == -1.):
            return self.constant * (math.log(nu2) -math.log(nu1))
        return self.constant * (math.pow(nu2, self.exponent+1)- math.pow(nu1, self.exponent+1)) / (self.exponent+1)
