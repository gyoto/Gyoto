'''Utility functions for Gyoto

Coord1dSet  -- easily initialize a gyoto.core.Coord1dset
rayTrace    -- a wrapper around gyoto.core.Scenery.rayTrace which
               hides most of its complexity
readScenery -- short hand for reading a Scenery XML file
'''


from gyoto import core
import numpy
import numbers
import keyword
import re

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
             height=None, width=None,
             quantities=None):
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
height, width -- vertical and horizontal resolution (override
           scenery.Screen.Resolution)
quantities -- quantities to compute as a space-separated string or
           iterable (tuplet, list...) of strings (overrides
           scenery.Quantities). Possible quantities:
               Intensity, EmissionTime, MinDistance, FirstDistMin,
               Redshift, NbCrossEqPlane, ImpactCoords,
               SpectrumStokesQ, SpectrumStokesU, SpectrumStokesV,
               Spectrum, SpectrumStokes, BinSpectrum, User1, User2,
               User3, User4, User5.

Output:
results -- dict containing the various requested quantities as per
           scenery.Quantities.

CAVEAT:
This high level-wrapper is Pythonic and take the arguments as j, i,
0-based indices whereas most Gyoto functions take them as i, j,
1-based indices.

Second form:

    '''
    ## Seconform: the C++ wrapper.
    # redistributes the parameters
    if isinstance(i, core.AstrobjProperties):
        ij=j
        aop=i
        if not isinstance(coord2dset, type):
            ipct=coord2dset
        else:
            ipct=None
        core._core.Scenery_rayTrace(sc, ij, aop, ipct)
        return

    ## First form: Pythonic wrapper.
    # If needed, read sc
    if type(sc) is str:
        sc=core.Factory(sc).scenery()

    # Determine resolution, width and height
    resbefore = sc.Screen.Resolution
    if width is None and height is None:
        res = resbefore
        height = res
        width = res
    elif height is None:
        height = width
        res = width
    elif width is None:
        width = height
        res = height
    elif width > height:
        res = width
    else:
        res = height

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

    if isinstance(coord2dset, core.Grid) and scalars == 0 :
        dims=(ny, nx)
    else:
        dims=(ntot,)

    # Set requested quantities if needed
    quantbefore = sc.Quantities
    if quantities is not None:
        if isinstance(quantities, str):
            sc.Quantities = quantities
        else:
            sc.Quantities = " ".join(quantities)

    # Prepare AstrobjProperties and dict of arrays to store results
    if sc.getSpectralQuantitiesCount():
        nsamples=sc.screen().spectrometer().nSamples()
    else:
        nsamples = None

    aop = core.AstrobjProperties(sc.Quantities, shape=dims, nsamples=nsamples)
    res = aop.data

    # Perform the actual ray-tracing
    sc.rayTrace(coord2dset, aop)

    if scalars == 2:
        for key in res:
            res[key]=res[key][...,0]

    # Reset Resolution and Quantities
    sc.Screen.Resolution = resbefore
    sc.Quantities = quantbefore
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

def valid_identifier(s):
    '''to_valid_identifier(str) -> identifier

    Bulds a valid identifier from any string.
    '''
    # If s is empty, replace with '_'
    if not s:
        s = '_'
    # Replace invalid characters with underscores.
    s = re.sub(r'[^a-zA-Z0-9_]', '_', s)
    # Add leading underscore if s starts with a digit.
    if s and s[0].isdigit():
        s = '_' + s
    # Check that s is not a reserved keyword:
    while keyword.iskeyword(s):
        return s + '_'
    return s

def make_class(namespace, classname, plugin=None, identifier=None, module=None,
               base=None):
    '''Dynamically wrap a Python class around a C++ class from a Gyoto plugin.

klass = make_class(namespace, classname, plugin=None, identifier=None, module=None])
obj = constructor()

Parameters:
namespace  -- one of gyoto.metric, gyoto.astrobj, gyoto.spectrum,
              gyoto.spectrometer.
classname  -- actual name of the class (e.g."KerrBL"), str.
plugin     -- optional name of the plugin to look in, e.g. "stdplug".
identifier -- optional identifier (for the docstring). Defaults to
              valid_identifier(classname).
module     -- optional module name (for the docstring).

Output:
klass      -- a Python class to instanciate the C++ class. 
'''
    if identifier is None:
        identifier=valid_identifier(classname)

    if module is None:
        if plugin is not None and plugin == valid_identifier(plugin):
            module = f'gyoto.{plugin}'
        else:
            module = 'gyoto.utils'

    if base is None:
        base=namespace.Generic

    plgname = plugin if plugin is not None else '(unknown)'

    doc = __doc__='''obj = '''+identifier+'''()
    Expose class \''''+classname+'''\' from Gyoto plugin \''''+plgname+'''\'.'''

    class klass(base):
        _classname  = classname
        _plugin     = plugin
        def __init__(self, *args, **kwargs):
            if self._plugin is None:
                plugins=()
            else:
                plugins=(self._plugin,)
            if args:
                super().__init__(*args, **kwargs)
            else:
                super().__init__(self._classname, plugins, **kwargs)

    klass.__name__   = identifier
    klass.__doc__    = doc
    klass.__module__ = module

    return klass

def convert(val, from_unit, to_unit, metric=None):
    ''' Convert values
    '''
    if from_unit==to_unit:
        return val

    meter=core.Unit('m')
    second=core.Unit('s')
    converter = None

    if to_unit=='geometrical' or to_unit=='geometrical_time':
        if from_unit=='geometrical_time' or from_unit == 'geometrical':
            return val;
        if metric is None:
            raise ValueError('metric is needed to convert to geometrical units')
        from_u = core.Unit(from_unit);
        from_factor=1.
        to_factor = 1./metric.unitLength()
        if core.areConvertible(from_u, meter):
            to_u = meter
        elif core.areConvertible(from_u, second):
            to_u = second
            to_factor *= core.GYOTO_C
        else: raise ValueError(f'cannot convert from {from_unit} to geometrical units')
        converter = lambda v: core.Converter(from_u, to_u)(v * from_factor) * to_factor;
    elif from_unit=='geometrical_time' or from_unit == 'geometrical':
        if metric is None:
            raise ValueError('metric is needed to convert to geometrical units')
        to_u = core.Unit(to_unit);
        to_factor=1.
        from_factor = metric.unitLength()
        if core.areConvertible(to_u, meter):
            from_u = meter
        elif core.areConvertible(to_u, second):
            from_u = second
            from_factor /= core.GYOTO_C
        else: raise ValueError(f'cannot convert from geometrical units to {to_unit}')
        converter = lambda v: core.Converter(from_u, to_u)(v * from_factor) * to_factor;
    else:
        from_u=core.Unit(from_unit)
        to_u=core.Unit(to_unit)
        if core.areConvertible(from_u, to_u):
            converter = core.Converter(from_u, to_u)
        else:
            if core.areConvertible(from_u, meter) and core.areConvertible(to_u, second):
                converter = lambda v: core.Converter(second, to_u)(core.Converter(from_u, meter)(v) / core.GYOTO_C)
            elif core.areConvertible(from_u, second) and core.areConvertible(to_u, meter):
                converter = lambda v: core.Converter(meter, to_u)(core.Converter(from_u, second)(v) * core.GYOTO_C)
            else:
                raise ValueError (f'cannot convert from {from_unit} to {to_unit}')

    try:
        return list([converter(val[i]) for i in range(len(val))])
    except TypeError:
        return converter(val)
