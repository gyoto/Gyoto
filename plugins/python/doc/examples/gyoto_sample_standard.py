'''Sample Astrobj::Standard for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Astrobj::Standard implementations using Gyoto's "python"
   plug-in. Note that this plug-in can be renamed to whatever matches
   the particular version of Python it has been built against
   (e.g. python3.4).

   The goal is to be able to instantiate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto
   gyoto.loadPlugin("python") # or python2.7 or python3.4...
   sp=gyoto.Astrobj("Python::Standard")
   sp.set("Module", "gyoto_sample_standard")
   sp.set("Class", "FixedStar")
   sp.set("Parameters", (0., 0., 0.))

   Classes that aim at implementing the Gyoto::Astrobj::Standard
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
