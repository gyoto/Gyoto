'''Utility functions for Gyoto

Coord1dSet  -- easily initialize a gyoto.core.Coord1dset
rayTrace    -- a wrapper around gyoto.core.Scenery.rayTrace which
               hides most of its complexity
readScenery -- short hand for reading a Scenery XML file
'''


from gyoto import core
import numpy
import numbers

def Coord1dSet(k, res, sz):
    '''Easily initialize a gyoto.core.Coord1dset

synopsis:
set = Coord1dSet(k, res, sz)

parameters:
k   -- something to convert to Coord1dSet
res -- Screen resolution
sz  -- Screen size in that direction

caveats:
Gyoto indices start at 1. This function takes care of tranlating from
the more standard Python 0-based indices to Gyoto 1-based indices.

returns:
if k is None:
  a gyoto.core.Range() covering all pixels according to res and sz

if k is a Python range:
  the corresponding gyoto.core.Range

if k is a scalar integer:
  a gyoto.core.Indices instance containing this value

if k is an array-like containing only integers:
  the corresponding gyoto.core.Indices instance

    '''
    data=None
    if (k is None):
        k=core.Range(res//2-sz//2+1, res//2-sz//2+sz, 1)
    elif type(k) is range:
        k=core.Range(k.start+1, k.stop, k.step)
    elif type(k) is slice:
        start=k.start
        if start is None:
            start=0
        elif start < 0:
            start=res+start
        start += 1
        stop=k.stop
        if stop is None:
            stop=res
        elif stop < 0:
            stop=res+stop
        step=k.step
        if step is None:
            step=1
        k=core.Range(start, stop, step)
    elif numpy.isscalar(k):
        if isinstance(k, numbers.Integral):
            data=numpy.array([k+1], numpy.uint64)
            k=core.Indices(data)
        elif isinstance(k, numbers.Real):
            data=numpy.array([k])
            k=core.Angles([k])
        else:
            raise ValueError('unknown scalar type')
    else:
        if all([isinstance(n, numbers.Integral) for n in k]):
            data=numpy.asarray(k, dtype=numpy.uint64)
            data += 1
            k=core.Indices(data)
        else:
            data=numpy.asarray(k)
            k=core.Angles(k)
    return k

def rayTrace(sc,
             j=None, i=None,
             coord2dset=core.Grid,
             prefix='\r j = ',
             height=None, width=None):
    '''Ray-trace scenery

First form:

results=scenery.rayTrace([j, i, [coord2dset [,prefix]]])

optional parameters:
j       -- vertical specification of the part of the field to trace
           (see gyoto.util.Coord1dSet)
i       -- horizontal specification of the part of the field to trace
           (see gyoto.util.Coord1dSet)
coord2dset -- a Coord2dSet subclass. Default: gyoto.core.Grid. The other
           value that makes sense is gyoto.core.Bucket.
prefix  -- prefix to be written in front of the row number for
           progress output
height, width -- vertical and horizontal resolution (overrides what
           is specified in scenery.screen().resolution()

Output:
results -- dict containing the various requested quantities as per
           scenery.requestedQuantitiesString().

CAVEAT:
This high level-wrapper is Pythonic and take the arguments as j, i,
0-based indices whereas most Gyoto functions take them as i, j,
1-based indices.

TODO:
Support impactcoords.

Second form:

    '''
    if isinstance(i, core.AstrobjProperties):
        ij=j
        aop=i
        if not isinstance(coord2dset, type):
            ipct=coord2dset
        else:
            ipct=None
        core._core.Scenery_rayTrace(sc, ij, aop, ipct)
        return

    # If needed, read sc
    if type(sc) is str:
        sc=core.Factory(sc).scenery()

    # Determine resolution, width and height
    res=sc.screen().resolution()
    if width is not None:
        res = width
    else:
        width=res
    if height is not None and height > width:
        res = height
    else:
        height=res
    sc.screen().resolution(res)

    # Prepare coord2dset
    scalars=int(numpy.isscalar(i))+int(numpy.isscalar(j))
    nx=None
    if isinstance(coord2dset, type):
        if not issubclass(coord2dset, core.Coord2dSet):
            raise TypeError("when coord2dset is a type, it must be a subclass of gyoto.core.Coord2dSet")
        if not isinstance(i, core.Coord1dSet):
            i=Coord1dSet(i, res, width)
        if not isinstance(j, core.Coord1dSet):
            j=Coord1dSet(j, res, height)
        try:
            coord2dset=coord2dset(i, j, prefix)
        except TypeError:
            coord2dset=coord2dset(i, j)
        nx=i.size()
        ntot=coord2dset.size()
        ny=ntot//nx
    elif isinstance(coord2dset, core.Coord2dSet):
        nx=ntot=coord2dset.size()
        ny=1
    else:
        raise TypeError('coord2dset must be a gyoto.core.Coord2dSet subclass or instance')

    if isinstance(coord2dset, core.Grid) and scalars is 0 :
        dims=(ny, nx)
        array_double_fromnumpy1or2=core.array_double.fromnumpy2
        array_double_fromnumpy2or3=core.array_double.fromnumpy3
    else:
        dims=(ntot,)
        array_double_fromnumpy1or2=core.array_double.fromnumpy1
        array_double_fromnumpy2or3=core.array_double.fromnumpy2

    # Prepare arrays to store results
    res = dict()
    aop=core.AstrobjProperties()
    aop.offset=ntot

    if sc.getSpectralQuantitiesCount():
        nsamples=sc.screen().spectrometer().nSamples()

    if 'Intensity' in sc.requestedQuantitiesString():
        intensity=numpy.zeros(dims)
        pintensity=array_double_fromnumpy1or2(intensity)
        aop.intensity=pintensity
        res['Intensity']=intensity

    if 'EmissionTime' in sc.requestedQuantitiesString():
        time=numpy.zeros(dims)
        ptime=array_double_fromnumpy1or2(time)
        aop.time=ptime
        res['EmissionTime'] = time

    if 'MinDistance' in sc.requestedQuantitiesString():
        distance=numpy.zeros(dims)
        pdistance=array_double_fromnumpy1or2(distance)
        aop.distance=pdistance
        res['MinDistance'] = distance

    if 'FirstDistMin' in sc.requestedQuantitiesString():
        first_dmin=numpy.zeros(dims)
        pfirst_dmin=array_double_fromnumpy1or2(first_dmin)
        aop.first_dmin=pfirst_dmin
        res['FirstDistMin'] = first_dmin

    if 'Redshift' in sc.requestedQuantitiesString():
        redshift=numpy.zeros(dims)
        predshift=array_double_fromnumpy1or2(redshift)
        aop.redshift=predshift
        res['Redshift'] = redshift

    if 'NbCrossEqPlane' in sc.requestedQuantitiesString():
        nbcrosseqplane=numpy.zeros(dims)
        pnbcrosseqplane=array_double_fromnumpy1or2(nbcrosseqplane)
        aop.nbcrosseqplane=pnbcrosseqplane
        res['NbCrossEqPlane'] = nbcrosseqplane

    if 'ImpactCoords' in sc.requestedQuantitiesString():
        impactcoords=numpy.zeros(dims+(16,))
        pimpactcoords=array_double_fromnumpy2or3(impactcoords)
        aop.impactcoords=pimpactcoords
        res['ImpactCoords'] = impactcoords

    if 'User1' in sc.requestedQuantitiesString():
        user1=numpy.zeros(dims)
        puser1=array_double_fromnumpy1or2(user1)
        aop.user1=puser1
        res['User1'] = user1
    
    if 'User2' in sc.requestedQuantitiesString():
        user2=numpy.zeros(dims)
        puser2=array_double_fromnumpy1or2(user2)
        aop.user2=puser2
        res['User2'] = user2

    if 'User3' in sc.requestedQuantitiesString():
        user3=numpy.zeros(dims)
        puser3=array_double_fromnumpy1or2(user3)
        aop.user3=puser3
        res['User3'] = user3

    if 'User4' in sc.requestedQuantitiesString():
        user4=numpy.zeros(dims)
        puser4=array_double_fromnumpy1or2(user4)
        aop.user4=puser4
        res['User4'] = user4

    if 'User5' in sc.requestedQuantitiesString():
        user5=numpy.zeros(dims)
        puser5=array_double_fromnumpy1or2(user5)
        aop.user5=puser5
        res['User5'] = user5

    if 'Spectrum' in sc.requestedQuantitiesString():
        spectrum=numpy.zeros((nsamples,)+dims)
        pspectrum=array_double_fromnumpy2or3(spectrum)
        aop.spectrum=pspectrum
        res['Spectrum'] = spectrum

    if 'SpectrumStokesQ' in sc.requestedQuantitiesString():
        stokesQ=numpy.zeros((nsamples,)+dims)
        pstokesQ=array_double_fromnumpy2or3(stokesQ)
        aop.stokesQ=pstokesQ
        res['SpectrumStokesQ'] = stokesQ

    if 'SpectrumStokesU' in sc.requestedQuantitiesString():
        stokesU=numpy.zeros((nsamples,)+dims)
        pstokesU=array_double_fromnumpy2or3(stokesU)
        aop.stokesU=pstokesU
        res['SpectrumStokesU'] = stokesU

    if 'SpectrumStokesV' in sc.requestedQuantitiesString():
        stokesV=numpy.zeros((nsamples,)+dims)
        pstokesV=array_double_fromnumpy2or3(stokesV)
        aop.stokesV=pstokesV
        res['SpectrumStokesV'] = stokesV

    if 'BinSpectrum' in sc.requestedQuantitiesString():
        binspectrum=numpy.zeros((nsamples,)+dims)
        pbinspectrum=array_double_fromnumpy2or3(binspectrum)
        aop.binspectrum=pbinspectrum
        res['BinSpectrum'] = binspectrum

    # Perform the actual ray-tracing
    sc.rayTrace(coord2dset, aop)

    if scalars is 2:
        for key in res:
            res[key]=res[key][0]

    return res

def Scenery_getitem(self, args):
    '''Shortcut for Scenery.rayTrace(i, j)
'''
    return self.rayTrace(*args)

def readScenery(filename):
    '''Read Scenery from XML file'''
    return core.Factory(filename).scenery()

def writeObject(obj, filename):
    '''Write Gyoto object (e.g. Scenery) to XML file'''
    core.Factory(obj).write(filename)

### Pythonic extension of methods
# Worldline
# This version of  getCartesian accepts numpy arrays
#   wl.getCartesian(t, x, y, z, [xprime, yprime, zprime])
# in addition to to the raw level arrays and dimension
def _Worldline_getCartesian(self, t=None, *arrays):
    # if no time provided, return all computed coordinates
    if t is None:
        n=self.get_nelements()
        t=numpy.empty(n)
        self.get_t(t)
        x=numpy.empty(n)
        y=numpy.empty(n)
        z=numpy.empty(n)
        xprime=numpy.empty(n)
        yprime=numpy.empty(n)
        zprime=numpy.empty(n)
        _Worldline_getCartesian(self, t, x, y, z, xprime, yprime, zprime)
        return t, x, y, z, xprime, yprime, zprime
    # if t is provided, return 8-coord for this time. t can be index or time.
    if not len(arrays):
        scalarcase=False
        if numpy.isscalar(t):
            scalarcase=True
            if numpy.issubdtype(type(t), numpy.integer):
                n=self.get_nelements()
                ta=numpy.empty(n)
                self.get_t(ta)
                if t >= 0:
                    t=ta[t]
                else:
                    t=ta[n+t]
            t=numpy.full(1, t)
        n=t.size
        x=numpy.empty(n)
        y=numpy.empty(n)
        z=numpy.empty(n)
        xprime=numpy.empty(n)
        yprime=numpy.empty(n)
        zprime=numpy.empty(n)
        _Worldline_getCartesian(self, t, x, y, z, xprime, yprime, zprime)
        if scalarcase:
            return numpy.asarray([t, x[0], y[0], z[0], xprime[0], yprime[0], zprime[0]])
        else:
            return t, x, y, z, xp, yp, zp
    # if t is a gyoto.core.array_double, call low-level C++ function
    if isinstance(t, core.array_double):
        core._core.Worldline_getCartesian(self, t, *arrays)
    # else fill pre-allocated arrays
    else:
        sizes=numpy.asarray([v.size for v in arrays])
        if not (sizes == t.size).all():
            raise ValueError('all arrays must be the same size')
        core._core.Worldline_getCartesian(
            self,
            core.array_double.fromnumpy1(t),
            t.size,
            *[core.array_double.fromnumpy1(v) for v in arrays]
            )

# Same for getCoord
def _Worldline_getCoord(self, t, *arrays):
    # if no time provided, return all computed coordinates
    if t is None:
        n=self.get_nelements()
        t=numpy.empty(n)
        self.get_t(t)
        x=numpy.empty(n)
        y=numpy.empty(n)
        z=numpy.empty(n)
        tdot=numpy.empty(n)
        xdot=numpy.empty(n)
        ydot=numpy.empty(n)
        zdot=numpy.empty(n)
        _Worldline_getCoord(self, t, x, y, z, tdot, xdot, ydot, zdot)
        return t, x, y, z, tdot, xdot, ydot, zdot
    # if t is provided, return 8-coord for this time. t can be index or time.
    if numpy.isscalar(t) and not len(arrays):
        if numpy.issubdtype(type(t), numpy.integer):
            n=self.get_nelements()
            ta=numpy.empty(n)
            self.get_t(ta)
            if t >= 0:
                t=ta[t]
            else:
                t=ta[n+t]
        ta=numpy.full(1, t)
        x=numpy.empty(1)
        y=numpy.empty(1)
        z=numpy.empty(1)
        tdot=numpy.empty(1)
        xdot=numpy.empty(1)
        ydot=numpy.empty(1)
        zdot=numpy.empty(1)
        _Worldline_getCoord(self, ta, x, y, z, tdot, xdot, ydot, zdot)
        return numpy.asarray([t, x[0], y[0], z[0], tdot[0], xdot[0], ydot[0], zdot[0]])
    # if t is a gyoto.core.array_double, call low-level C++ function
    if isinstance(t, core.array_double) or numpy.isscalar(t) :
        core._core.Worldline_getCoord(self, t, *arrays)
    # else fill pre-allocated arrays
    else:
        sizes=numpy.asarray([v.size for v in arrays])
        if not (sizes == t.size).all():
            raise ValueError('all arrays must be the same size')
        core._core.Worldline_getCoord(
            self,
            core.array_double.fromnumpy1(t),
            t.size,
            *[core.array_double.fromnumpy1(v) for v in arrays]
            )
