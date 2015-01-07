#/bin/env python
# -*- coding: utf-8 -*-
# Example file for gyoto
# Before loading the gyoto_std extension:
#  - gyoto must be imported;
#  - the stdplug Gyoto plug-in must be loaded. This normally happens
#    automatically, unless the user has set her environment otherwise.
import numpy
import matplotlib as ml
import matplotlib.pyplot as plt
import gyoto
gyoto.loadPlugin("stdplug")
import gyoto_std

a=gyoto.Factory("../doc/examples/example-moving-star.xml")
sc=a.getScenery()
sc.nThreads(8)
sc.astrobj().opticallyThin(False)

# Trace and plot NULL geodesic:

ph=gyoto.Photon()
ph.setInitialCondition(sc.metric(), sc.astrobj(), sc.screen(), 0., 0.)
ph.hit()
n=ph.get_nelements()
t=gyoto.array_double(n)
r=gyoto.array_double(n)
theta=gyoto.array_double(n)
phi=gyoto.array_double(n)
ph.getCoord(t, r, theta, phi)

t2=numpy.zeros(n)
r2=numpy.zeros(n)
theta2=numpy.zeros(n)
phi2=numpy.zeros(n)

for i in range(0, n):
    t2[i]=t[i]
    r2[i]=r[i]
    theta2[i]=theta[i]
    phi2[i]=phi[i]

plt.plot(t2, r2)
plt.show()

# Trace and plot timelike geodesic
# We need to cast the object to a gyoto_std.Star:

wl=gyoto_std.Star(sc.astrobj())
wl.xFill(1000)

n=wl.get_nelements()
t=gyoto.array_double(n)
r=gyoto.array_double(n)
theta=gyoto.array_double(n)
phi=gyoto.array_double(n)
wl.getCoord(t, r, theta, phi)

t2=numpy.zeros(n)
r2=numpy.zeros(n)
theta2=numpy.zeros(n)
phi2=numpy.zeros(n)

for i in range(0, n):
    t2[i]=t[i]
    r2[i]=r[i]
    theta2[i]=theta[i]
    phi2[i]=phi[i]

plt.plot(r2*numpy.cos(phi2), r2*numpy.sin(phi2))
plt.show()

# Ray-trace scenery

res=sc.screen().resolution()
intensity=numpy.zeros((res, res), dtype=float)
time=numpy.zeros((res, res), dtype=float)
distance=numpy.zeros((res, res), dtype=float)
aop=gyoto.AstrobjProperties()
aop.Intensity(intensity)
aop.EmissionTime(time)
aop.MinDistance(distance)

ii=gyoto.Range(1, res, 1)
jj=gyoto.Range(1, res, 1)
grid=gyoto.Grid(ii, jj, "\rj = ")

sc.rayTrace(grid, aop)

plt.imshow(intensity)
plt.show()
plt.imshow(time)
plt.show()
plt.imshow(distance)
plt.show()

# Another Scenery, with spectrum

sc=gyoto.Factory("../doc/examples/example-polish-doughnut.xml").getScenery()
sc.screen().resolution(32)
res=sc.screen().resolution()
ns=sc.screen().spectrometer().nSamples()
spectrum=numpy.zeros((ns, res, res), dtype=float)

ii=gyoto.Range(1, res, 1)
jj=gyoto.Range(1, res, 1)
grid=gyoto.Grid(ii, jj, "\rj = ")

aop=gyoto.AstrobjProperties()
aop.Spectrum(spectrum)
aop.offset=res*res

sc.rayTrace(grid, aop)

plt.imshow(spectrum[1,:,:])
plt.show()

# Another Scenery, with impact coords, created from within Python

met=gyoto.Metric("KerrBL")
met.mass(4e6, "sunmass")
ao=gyoto.Astrobj("PageThorneDisk")
ao.metric(met)
ao.opticallyThin(False)
ao.rMax(100)
screen=gyoto.Screen()
screen.distance(8, "kpc")
screen.time(8, "kpc")
screen.resolution(64)
screen.inclination(numpy.pi/4)
screen.PALN(numpy.pi)
screen.time(8, "kpc")
screen.fieldOfView(100, "µas")
sc=gyoto.Scenery()
sc.metric(met)
sc.astrobj(ao)
sc.screen(screen)
sc.delta(1, "kpc")
sc.adaptive(True)
sc.nThreads(8)

res=sc.screen().resolution()

ii=gyoto.Range(1, res, 1)
jj=gyoto.Range(1, res, 1)
grid=gyoto.Grid(ii, jj, "\rj = ")

ipct=numpy.zeros((res, res, 16), dtype=float)

aop=gyoto.AstrobjProperties()
aop.ImpactCoords(ipct)
aop.offset=res*res

sc.rayTrace(grid, aop)

plt.imshow(ipct[:,:,0], interpolation="nearest", vmin=-100, vmax=0)
plt.show()

# Trace one line of the above using alpha and delta

N=10

buf=numpy.linspace(screen.fieldOfView()*-0.5, screen.fieldOfView()*0.5, N)
a=gyoto.Angles(buf)
d=gyoto.RepeatAngle(screen.fieldOfView()*-0.5, N)
bucket=gyoto.Bucket(a, d)

ipct=numpy.zeros((N, 1, 16), dtype=float)

aop=gyoto.AstrobjProperties()
aop.ImpactCoords(ipct)
aop.offset=N

sc.rayTrace(bucket, aop)
plt.plot(buf, ipct[:,0,0])
plt.show()

# Trace the diagonal of the above using i and j. The Range and Indices
# definitions below are equivalent.  Range is more efficient for a
# range, Indices can hold arbitrary indices.

ind=numpy.arange(1, res+1, dtype=numpy.uintp) # on 64bit arch...

ii=gyoto.Indices(ind)
jj=gyoto.Range(1, res, 1)
bucket=gyoto.Bucket(ii, jj)

ipct=numpy.zeros((res, 1, 16), dtype=float)

aop=gyoto.AstrobjProperties()
aop.ImpactCoords(ipct)
aop.offset=res

sc.rayTrace(bucket, aop)

t=numpy.clip(ipct[:,0,0], a_min=-200, a_max=0)
plt.plot(t)
plt.show()

# Any derived class can be instanciated from its name, as soon as the
# corresponding plug-in has been loaded into Gyoto. The standard
# plug-in is normally loaded automatically, but this can also be
# forced with gyoto.loadPlugin():
gyoto.loadPlugin('stdplug')
tt=gyoto.Astrobj('Torus')
kerr=gyoto.Metric('KerrBL')

# Most properties that can be set in an XML file can also be accessed
# from Python using the Property/Value mechanism:
p=tt.property("SmallRadius")
tt.set(p, gyoto.Value(0.2))
tt.get(p).toDouble() == 0.2

p=kerr.property("Spin")
kerr.set(p, gyoto.Value(0.95))
kerr.get(p).toDouble() == 0.95

# However, we also have Python extensions around the standard Gyoto
# plug-ins. Beware that the plug-in must be loaded into Gyoto before
# importing the corresponding Python extension:
gyoto.loadPlugin('stdplug')
import gyoto_std
# And if the lorene plug-in has been compiled:
# gyoto.loadPlugin('lorene')
# import gyoto_lorene

# It then becomes possible to access the methods specific to derived
# classes. They can be instanciated directly from the gyoto_* extension:
tr2=gyoto_std.Torus()
# and we can cast a generic pointer (from the gyoto extension) to a
# derived class:
tr=gyoto_std.Torus(tt)
p=tt.property("SmallRadius")
tt.get(p).toDouble() == tr.smallRadius()

# Another example: using a complex (i.e. compound) Astrobj:
cplx=gyoto_std.ComplexAstrobj()
cplx.append(tr)
cplx.append(sc.astrobj())
sc.astrobj(cplx)

cplx.rMax(100)

print("All done, exiting")
