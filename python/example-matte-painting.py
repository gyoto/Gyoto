import numpy
from matplotlib import pyplot as plt
import urllib
import PIL

import gyoto
from gyoto.matte_painting import *

### First, define an appropriate Scenery:

# Choose a metric, set its mass:
earth_mass = 5.9722e24 # mass of the Earth, in kg
kerrbl = gyoto.std.KerrBL()
kerrbl.Spin = 0.
kerrbl.Mass = 1*earth_mass, "kg"
kerrks = gyoto.std.KerrKS()
kerrks.Mass = kerrbl.Mass
kerrks.Spin = kerrbl.Spin
minkowski = gyoto.std.Minkowski()
minkowski.Mass=kerrbl.Mass
minkowski.Spherical = True
metric=kerrks

# Define the screen (=camera):
res=512 # resolution of the output image
scr=gyoto.core.Screen()
scr.Inclination = numpy.pi/2.
scr.PALN = numpy.pi
scr.Argument = numpy.pi/2.
scr.Metric = metric
scr.Resolution = res
scr.AngleKind = "Rectilinear"
scr.Distance = 2, "m" # camera is at 2m from the BH
# scr.FieldOfView = numpy.pi/2.
# or: ~40°, 36mm film behind a 50mm lense
scr.FieldOfView = 2.*numpy.atan(36./100.) 

# Use an empty Astrobj, just set rmax to 10 meters:
ao = gyoto.std.ComplexAstrobj()
ao.Metric = metric
ao.RMax= 10., "m"

# Build a Scenery with this Metric, Screen and Astrobj:
sc=gyoto.core.Scenery()
sc.Metric = metric
sc.Astrobj=ao
sc.Screen = scr
sc.NThreads = 16

##/ Second, choose a painter to paint the sky:
# To use a jpeg file containing a full-sky, 360°x180° panorama:
url = ("http://farm1.staticflickr.com/192/"
       "456185667_adde9d2f8e_o_d.jpg")
fp = urllib.request.urlopen(url)
img = numpy.array(PIL.Image.open(fp))
painter=Panorama(img=img)
nx=res
ny=res*10//16; # use a 16/10 aspect ratio
# Or a more conventional picture:
# url = ("https://i.pinimg.com/236x/19/53/11/"
#        "195311bdb5b029cb72214b95833e6dcf.jpg")
# fp = urllib.request.urlopen(url)
# img = numpy.array(PIL.Image.open(fp))
# painter = Picture(img=img, fov=scr.FieldOfView)
# ny = res
# nx = 3*res//4
# To use a p-mode-like pattern:
# painter=PMode(ntheta=80, nphi=80);
# nx=res
# ny=res*10//16; # use a 16/10 aspect ratio


##/ Third, simply call the appropriate function:
# matte_paint() can be called directly on sc, but to specify only a
# region or to call matte_paint repeatedly (e.g. on distinct
# painters, of to fiddle with yaw, pitch & roll), we can precompute
# ipct:
x1=(res-nx)//2+1;
x2=res+1-x1;
y1=(res-ny)//2+1;
y2=res+1-y1;
ii=gyoto.core.Range(x1, x2, 1)
jj=gyoto.core.Range(y1, y2, 1)
grid=gyoto.core.Grid(ii, jj, "\rj = ")
ipct=numpy.zeros((ny, nx, 16), dtype=float)
aop=gyoto.core.AstrobjProperties()
aop.impactcoords=gyoto.core.array_double.fromnumpy3(ipct)
aop.offset=nx*ny
sc.rayTrace(grid, aop)

#ipct=sc(x1:x2,y1:y2,impactcoords=);
#ipct=fits_read("ipct.fits");
#bg=matte_paint(ipct, painter, kind=metric.coordkind());
#
bg=matte_paint(ipct, painter, coordkind=metric.coordKind(), yaw=0., pitch=0., roll=0);

##/ Fourth, display the image:
plt.imshow(bg)
plt.show()
# limits, square=1;
