'''Sample Astrobj::Standard for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Astrobj::Standard implementations using Gyoto's "python"
   plug-in. Note that this plug-in can be renamed to whatever matches
   the particular version of Python it has been built against
   (e.g. python3.4).

   The goal is to be able to instantiate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto.core
   gyoto.core.requirePlugin("python") # or python2.7 or python3.4...
   sp=gyoto.core.Astrobj("Python::Standard")
   sp.set("Module", "gyoto_sample_standard")
   sp.set("Class", "FixedStar")
   sp.set("Position", (0., 0., 0.))
   sp.set("Radius", 1.)

   Classes that aim at implementing the Gyoto::Astrobj::Standard
   interface do so by providing the following methods and variable:

   __call__: required
   getVelocity: required
   giveDelta, emission, integrateEmission, transmission, set, get,
     properties: optional.
   emission and integrateEmission can be overloaded by using the
   varargs argument.

   Ad-hoc properties are declared in member properties which is a
   dict in the form:
     properties = { Key1: type1, Key2: type2}
   where Keyi and typei are strings.

   They are then handled by set(self, key, val) and get(self, key).

'''

import math
import numpy
import gyoto.core

class FixedStar:
    ''' Sample class for Astrobj::Python::Standard
    '''

    properties={"Position": "vector_double",
                "Radius": "double",
                "Spectrum": "spectrum",
                "Opacity": "spectrum"}
    '''Properties handled by set() and get()
    '''

    position=numpy.zeros(3, dtype=float)
    '''3 coordinates of the center of the star
    '''

    radius=1.
    '''Radius of the Star in geometrical units
    '''

    spectrum=gyoto.core.Spectrum("BlackBody")
    opacity=gyoto.core.Spectrum("PowerLaw")

    def set(self, key, val):
        '''Set parameters as ad hoc entities

        Python:
          instance.set("Position", (x1, x2, x3))
          instance.set("Radius", radius)

        XML:
          <Position> x1 x2 x3 </Position>
          <Radius> radius </Radius>

        Setting Radius actually sets CriticalValue and SafetyValue
        (Properties of Astrobj::Standard).

        '''
        if (key=="Position"):
            self.position = val
        elif key == "Radius":
            self.radius = val
        elif key == "Spectrum":
            self.spectrum = val
        elif key == "Opacity":
            self.opacity = val
        # Note: since we set attributes here, __setattr__(self, key,
        # val) defined below will be called and implements useful side
        # effects.

    def get(self, key):
        '''Retrieve properties

        Optional, but simple to implement and very usefull. Required
        to write back to XML, which is used for MPI parallelization.
        '''
        if key == "Position":
            return self.position
        if key == "Radius":
            return self.radius
        if key == "Spectrum":
            return self.spectrum
        if key == "Opacity":
            return self.opacity

    def __setattr__(self, key, val):
        '''Set attributes

        Optional, but allows reacting when attributes (self.xxx)
        change. In particular, Gyoto will set self.this to the C++
        class instance upon initialization.
        '''
        # First, actually store the attribute. This is what would
        # happen if we did not overload __setattr__.
        self.__dict__[key]=val
        if key == "this":
            # Then, if key is "this", trigger the side effects of
            # setting position and radius. This will reenter
            # __setattr__ with the following cases.
            self.position=self.position
            self.radius=self.radius
        elif key == "position":
            # when position is set, cache Cartesian expression
            pos=numpy.zeros(4, dtype=float)
            pos[1:]=val
            self.coord_st=self.to_cartesian(pos)
        elif key == "radius":
            # Update two properties from the Gyoto::Standard interface
            self.this.set("CriticalValue", val**2)
            self.this.set("SafetyValue", val**2*1.1+0.1)

    def to_cartesian(self, coord):
        '''Helper function, not in the API

        '''
        gg=self.this.metric()
        spherical=False
        if gg is not None:
            spherical = gg.coordKind() == gyoto.core.GYOTO_COORDKIND_SPHERICAL
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
        return dx*dx + dy*dy + dz*dz

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
        if self.this.opticallyThin():
            return self.spectrum(nuem, self.opacity(nuem), dsem)
        return self.spectrum(nuem)

    def transmission(self, nuem, dsem, cph, co):
        if (not self.this.opticallyThin()):
            return 0.
        opac=self.opacity(nuem)
        if opac == 0.:
            return 1.
        return math.exp(-opac*dsem)
