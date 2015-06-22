'''Sample metrics for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Metric implementations using Gyoto's "python" plug-in. Note that
   this plug-in can be renamed to whatever matches the particular
   version of Python it has been built against (e.g. python3.4).

   The goal is to be able to instanciate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto
   gyoto.loadPlugin("python") # or python2.7 or python3.4...
   gg=gyoto.Metric("Python")
   gg.set("Module", "gyoto_sample_metric")
   gg.set("Class", "Minkowski")

   Classes that aim at implementing the Gyoto::Metric::Generic
   interface do so by providing the following methods:

   gmunu(self, dst, pos): mandatory;
   christoffel(self, dst, pos): mandatory
   __setattr__(self, key, value): optional
   __setitem__(self, key, value): optional

'''

import math
import numpy
import gyoto

class Minkowski:
    '''Flat space metric

    Implemented for both Cartesian and spherical coordinates.

    Every Gyoto Metric implemented in Python must implement the three
    methods illustrated here.

    '''
    def __setattr__(self, key, value):
        '''Set attributes.

        Optional.

        C++ will set several attributes. By overloading __setattr__,
        on can react when that occurs, in particular to make sure this
        knows the coordinate kind as in this example.

        Attributes set by the C++ layer:

          this: if the Python extension "gyoto" can be imported, it
                will be set to a gyoto.Metric instance pointing to the
                C++-side instance. If the "gyoto" extension cannot be
                loaded, this will be set to None.

          spherical: when the spherical(bool t) method is called in
                the C++ layer, it sets the spherical attribute in the
                Python side.

          mass: when the mass(double m) method is called in the C++
                side, it sets the spherical attribute in the Python
                side.

        This example initializes coordKind in the C++ side if it is
        not already set, since this Minkowski class can work in
        either.

        '''
        # First, actually store the attribute. This is what would
        # happen if we did not overload __setattr__.
        self.__dict__[key]=value
        # Then, if key is "this", ensure this knows a valid coordKind.
        if (key is "this"):
            cK=value.coordKind()
            if cK is gyoto.GYOTO_COORDKIND_UNSPECIFIED:
                value.set("Spherical", False)
            # We could do without this, since this will tell us later
            # anyway.
            else:
                self.spherical = (cK is gyoto.GYOTO_COORDKIND_SPHERICAL)

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
