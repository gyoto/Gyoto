'''Sample metrics for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Metric implementations using Gyoto's "python" plug-in. Note that
   this plug-in can be renamed to whatever matches the particular
   version of Python it has been built against (e.g. python3.4).

   The goal is to be able to instantiate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto.core
   gyoto.core.requirePlugin("python") # or python2.7 or python3.4...
   gg=gyoto.core.Metric("Python")
   gg.set("Module", "gyoto_sample_metric")
   gg.set("Class", "Minkowski")

   Classes that aim at implementing the Gyoto::Metric::Generic
   interface do so by providing the following methods:

   gmunu(self, dst, pos): mandatory;
   christoffel(self, dst, pos): mandatory;
   getRms(self): optional
   getRmb(self): optional
   getSpecificAngularMomentum(self, rr): optional
   getPotential(self, pos, l_cst): optional
   __setattr__(self, key, value): optional, useful to react
     when the C++ Metric object sets attributes:
          this: if the Python extension "gyoto.core" can be imported,
                it will be set to a gyoto.core.Metric instance
                pointing to the C++-side instance. If the "gyoto.core"
                extension cannot be loaded, this will be set to None.

          spherical: when the spherical(bool t) method is called in
                the C++ layer, it sets the spherical attribute in the
                Python side.

          mass: when the mass(double m) method is called in the C++
                side, it sets the spherical attribute in the Python
                side.
   __setitem__(self, key, value): optional, mandatory to support the
     "Parameters" Property to set arbitrary parameters for this
     metric.
   set(self, key, val) and get(self, key): optional, both need to be
     implemented in order to support additional "Properties" with
     arbitrary names. In addition, member self.properties must be set
     to a dictionary listing key:datatype pairs, where both key and
     datatype are strings. As of writing, only "double" is supported
     as a datatype.
'''

import math
import numpy
import gyoto.core

class Minkowski:
    '''Flat space metric

    Implemented for both Cartesian and spherical coordinates.
    '''
    def __setattr__(self, key, value):
        '''Set attributes.

        Optional.

        C++ will set several attributes. By overloading __setattr__,
        one can react when that occurs, in particular to make sure `this'
        knows the coordinate kind as in this example.

        Attributes set by the C++ layer:

          this: if the Python extension "gyoto.core" can be imported, it
                will be set to a gyoto.core.Metric instance pointing to the
                C++-side instance. If the "gyoto.core" extension cannot be
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
        if (key == "this"):
            cK=value.coordKind()
            if cK is gyoto.core.GYOTO_COORDKIND_UNSPECIFIED:
                value.set("Spherical", False)
            # We could do without this, since this will tell us later
            # anyway.
            else:
                self.spherical = (cK is gyoto.core.GYOTO_COORDKIND_SPHERICAL)

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

class KerrBL:
    '''A Python implementation of Gyoto::Metric::KerrBL

    Spin and HorizonSecurity may be set either using the Parameters
    Property, as in:
        <Parameters> O.5 0.01 </Parameters>
    or using two distinct ad hoc "properties":
        <Spin> 0.5 </Spin>
        <HorizonSecurity> 0.01 </HorizonSecurity>

    The various methods are essentially cut-and-paste from the C++ class.

    Only gmunu and christoffel absolutely need to be implemented, but
    various Astrobj require more. See below.

    This is for educational and testing purposes. For all practical uses,
    the C++ implementation is much faster.

    '''

    # This is needed to support setting Spin and HorizonSecurity by name
    # in addition to the set and get methods
    properties={"Spin":"double", "HorizonSecurity":"double"}

    spin=0.
    drhor=gyoto.core.GYOTO_KERR_HORIZON_SECURITY
    rsink=2.+gyoto.core.GYOTO_KERR_HORIZON_SECURITY

    def __setitem__(self, key, value):
        '''Set parameters

        Optional, one way to handle parameters

        This is how Gyoto sends the <Parameters/> XML entity:
        metric[key]=value
        key=0: set spin (0.)
        key=1: set drhor, thickness of sink layer around horizon (0.01)
        '''
        if (key==0):
            self.spin = value
        elif (key==1):
            self.drhor = value
        else:
            raise IndexError

    def set(self, key, value):
        '''Set parameters by name

        Optional, one way to handle parameters

        This is how Gyoto sends custom  XML entities, e.g.
          <MyEntity>value</MyEntity>
        metric.set("MyEntity", value)

        Here:
        key="Spin": set spin (0.)
        key="HorizonSecurity": set drhor, thickness of sink layer
             around horizon (0.01)
        '''
        if (key=="Spin"):
            self.spin = value
        elif (key=="HorizonSecurity"):
            self.drhor = value
        else:
            raise IndexError

    def get(self, key):
        '''Get parameters by name

        Optional, mandatory if "properties" and "set" are defined.
        '''
        if (key=="Spin"):
            return self.spin
        elif (key=="HorizonSecurity"):
            return self.drhor
        else:
            raise IndexError

    def __setattr__(self, key, value):
        '''Set attributes.

        Optional.

        C++ will set several attributes. By overloading __setattr__,
        one can react when that occurs, in particular to make sure
        `this' knows the coordinate kind as in this example.

        In addition, update self.rsink each time self.spin or
        self.drhor change.
        '''
        # First, actually store the attribute. This is what would
        # happen if we did not overload __setattr__.
        self.__dict__[key]=value
        # Then, if key is "this", ensure `this' knows a valid coordKind.
        if (key == "this"):
            self.this.set("Spherical", True)
        elif key in ("spin", "drhor"):
            if self.spin > 1:
                self.rsink=drhor
            else:
                self.rsink=1.+math.sqrt(1.-self.spin**2)+self.drhor

    def gmunu(self, g, x):
        ''' Gyoto::Metric::Generic::gmunu(double g[4], const double pos[4])

        Mandatory.

        C++ will send two NumPy arrays.

        Note that the user will not be able to call this method directly but
        through self.this.gmunu which has a different calling sequence:
          g=self.this.gmunu(x)
        '''
        spin_=self.spin
        a2_=spin_**2
        r=x[1]
        sth2=math.sin(x[2])**2
        cth2=math.cos(x[2])**2
        r2=r*r
        sigma=r2+a2_*cth2
        delta=r2-2.*r+a2_

        for mu in range(0, 4):
            for nu in range(0, 4):
                g[mu][nu]=g[nu][mu]=0

        g[0][0] = -1.+2.*r/sigma;
        g[1][1] = sigma/delta;
        g[2][2] = sigma;
        g[3][3] = (r2+a2_+2.*r*a2_*sth2/sigma)*sth2;
        g[0][3] = g[3][0] = -2*spin_*r*sth2/sigma;

    def christoffel(self, dst, x):
        '''Gyoto::Metric::Generic::christoffel(double dst[4][4][4], const double pos[4])

        Mandatory.

        C++ will send two NumPy arrays.

        Like gmunu, the call will actually be:
          dst=metric.gmunu(x)
        where `metric' is self.this.

        '''
        for alpha in range(0, 4):
            for mu in range(0, 4):
                for nu in range(0, 4):
                    dst[alpha][mu][nu]=0.

        spin_=self.spin
        a2_=spin_**2
        r=x[1]
        sth=math.sin(x[2])
        cth=math.cos(x[2])
        sth2 = sth*sth
        cth2 = cth*cth
        sth4=sth2*sth2
        s2th = 2.*sth*cth
        c2th=cth2-sth2
        s4th = 2.*s2th*c2th
        s2th2= s2th*s2th
        ctgth=cth/sth
        r2=r*r
        r4=r2*r2
        r6=r4*r2;
        Sigma=r2+a2_*cth2
        Sigma2=Sigma*Sigma
        Delta=r2-2.*r+a2_
        Deltam1=1./Delta
        Sigmam1=1./Sigma
        Sigmam2=Sigmam1*Sigmam1
        Sigmam3=Sigmam2*Sigmam1
        a2cthsth=a2_*cth*sth
        rSigmam1=r*Sigmam1
        Deltam1Sigmam2=Deltam1*Sigmam2
        r2plusa2 = r2+a2_

        dst[1][1][1]=(1.-r)*Deltam1+rSigmam1;
        dst[1][2][1]=dst[1][1][2]=-a2cthsth*Sigmam1;
        dst[1][2][2]=-Delta*rSigmam1;
        dst[1][3][3]=-Delta*sth2*(r+(a2_*(-2.*r2+Sigma)*sth2)/Sigma2)/Sigma;
        dst[1][3][0]=dst[1][0][3]=spin_*Delta*(-2*r2+Sigma)*sth2*Sigmam3;
        dst[1][0][0]=-Delta*(-2.*r2+Sigma)*Sigmam3;
        dst[2][1][1]=a2cthsth*Deltam1*Sigmam1;
        dst[2][2][1]=dst[2][1][2]=rSigmam1;
        dst[2][2][2]=-a2cthsth*Sigmam1;
        dst[2][3][3]=-sth*cth*Sigmam3 * (Delta*Sigma2 + 2.*r*r2plusa2*r2plusa2);
        dst[2][0][3]=dst[2][3][0]=spin_*r*r2plusa2*s2th*Sigmam3;
        dst[2][0][0]=-2.*a2cthsth*r*Sigmam3;
        dst[3][3][1]=dst[3][1][3]=Deltam1*Sigmam2 * (r*Sigma*(Sigma-2.*r) + a2_*(Sigma-2.*r2)*sth2);
        dst[3][3][2]=dst[3][2][3]=Sigmam2*ctgth * (-(Sigma+Delta)*a2_*sth2 + r2plusa2*r2plusa2);
        dst[3][0][1]=dst[3][1][0]=spin_*(2.*r2-Sigma)*Deltam1Sigmam2;
        dst[3][0][2]=dst[3][2][0]=-2.*spin_*r*ctgth*Sigmam2;
        dst[0][3][1]=dst[0][1][3]=-spin_*sth2*Deltam1Sigmam2 * (2.*r2*r2plusa2 + Sigma*(r2-a2_));
        dst[0][3][2]=dst[0][2][3]=Sigmam2*spin_*a2_*r*sth2*s2th;
        dst[0][0][1]=dst[0][1][0]=(a2_+r2)*(2.*r2-Sigma)*Deltam1Sigmam2;
        dst[0][0][2]=dst[0][2][0]=-a2_*r*s2th*Sigmam2;

        return 0

    def getRms(self):
        aa=self.spin;
        a2_=aa*aa
        z1 = 1. + pow((1. - a2_),1./3.)*(pow((1. + aa),1./3.) + pow((1. - aa),1./3.))
        z2 = pow(3.*a2_ + z1*z1,1./2.);

        return (3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),1./2.));

    def getRmb(self):
        spin_=self.spin
        return 2.-spin_+2.*math.sqrt(1.-spin_);

    def getSpecificAngularMomentum(self, rr):
        aa=self.spin
        sqrtr=math.sqrt(rr);
        return (rr*rr-2.*aa*sqrtr+aa*aa)/(rr**1.5-2.*sqrtr+aa);

    def getPotential(self, pos, l_cst):
        # this is W = -ln(|u_t|) for a circular equatorial 4-velocity
        # Don't try to call directly gmunu from `self',
        # Call it through `this', a pointer to the C++ instance of
        # Gyoto::Metric::Python
        # Note that the calling sequence is not gmunu(self, gg, pos)
        gg=self.this.gmunu(pos)
        gtt = gg[0,0];
        gtp = gg[0,3];
        gpp = gg[3,3];
        Omega = -(gtp + l_cst * gtt)/(gpp + l_cst * gtp) ;
  
        W = (0.5 * math.log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp))
             - math.log(abs(gtt + Omega * gtp))) ;

        return  W ;

    def isStopCondition(self, coord):
        return coord[1] < self.rsink

    def circularVelocity(self, coor, vel, d):
        '''Velocity of a particle in circular motion at point

        Note: this is called only if this->keplerian_ is False.
        '''

        sinth = math.sin(coor[2]);
        coord = [coor[0], coor[1]*sinth, math.pi*0.5, coor[3]];

        vel[1] = vel[2] = 0.;
        #vel[3] = 1./((d*pow(coord[1], 1.5) + spin_)*sinth);
        vel[3] = 1./((d*pow(coord[1], 1.5) + self.spin));

        vel[0] = self.this.SysPrimeToTdot(coor, vel[1:]);
        vel[3] *= vel[0];
