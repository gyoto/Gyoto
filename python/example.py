import numpy
import matplotlib as ml
import matplotlib.pyplot as plt
import gyoto

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

wl=gyoto.castToWorldline(sc.astrobj());
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
aop=gyoto.Properties()
aop.Intensity(intensity)
aop.EmissionTime(time)
aop.MinDistance(distance)

ii=gyoto.Range(1, res, 1)
jj=gyoto.Range(1, res, 1)
grid=gyoto.Grid(ii, jj)

sc.rayTrace(grid, aop)

plt.imshow(intensity)
plt.show()
plt.imshow(time)
plt.show()
plt.imshow(distance)
plt.show()
