#/bin/env python
# -*- coding: utf-8 -*-
# Example file for gyoto
#
# Copyright 2026 Thibaut Paumard
#
# This file is part of Gyoto.
#
# Gyoto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gyoto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import numpy
import matplotlib as ml
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import urllib.request
import PIL

import gyoto
from gyoto.matte_painting import *
# Parse command line and optionally switch to PDF output

pdfname=None
dir_path = os.path.dirname(os.path.realpath(__file__))
examples_dir=dir_path+"/../doc/examples/"
painters=['PMode', 'Picture', 'Panorama']
res=512

for param in sys.argv:
    sparam=param.split("=")
    if os.path.basename(sparam[0])==os.path.basename(__file__):
        pass
    elif sparam[0]=="--pdf":
        if len(sparam)==2:
            pdfname=sparam[1]
        else:
            raise ValueError('--pdf argument expects a filename, e.g. --pdf=output.pdf')
    elif sparam[0]=="--examples-dir":
        if len(sparam)==2:
            examples_dir=sparam[1]
        else:
            raise ValueError('--examples_dir argument expects a directory, e.g. --examples-dir=../doc/examples')
    elif sparam[0]=="--painters":
        if len(sparam)==2:
            painters=sparam[1].split(",")
        else:
            raise ValueError('--painters argument expects a '
            'comma-separated list of Painter subclasses to '
            'illustrate, e.g. --painters=PMode,Picture')
    elif sparam[0]=="--resolution":
        if len(sparam)==2:
            res=int(sparam[1])
        else:
            raise ValueError('--resolution argument expects an integer value')
    else:
        raise ValueError(f'unknown argument: {sparam[0]}')

pdf=None if pdfname is None else PdfPages(pdfname)
if len(examples_dir) > 0 and examples_dir[-1] != "/":
    examples_dir += "/"


### First, define an appropriate Scenery:

# Let's define a bunch of metrics and set their mass:
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
scr=gyoto.core.Screen()
scr.Inclination = numpy.pi/2.
scr.PALN = numpy.pi
scr.Argument = numpy.pi/2.
scr.Metric = metric
scr.Resolution = res
scr.AngleKind = "Rectilinear"
scr.Distance = 2, "m" # camera is at 2m from the BH
# ~40°, 36mm film behind a 50mm lense
scr.FieldOfView = 2.*numpy.atan(36./100.) 

# Use an empty Astrobj, just set rmax to 10 meters:
ao = gyoto.std.ComplexAstrobj()
ao.Metric = metric
ao.RMax= 10., "m"

# Build a Scenery with this Screen and Astrobj
# (the Metric will be set later):
sc=gyoto.core.Scenery()
sc.Metric = metric
sc.Astrobj=ao
sc.Screen = scr
sc.NThreads = 16

### We will illustrate the use of three different
### matte_painting.Painter subclasses.

## A p-mode-like pattern:
if 'PMode' in painters:
    painter=PMode(ntheta=80, nphi=80);
    nx = res
    ny = (res*10)//16; # use a 16/10 aspect ratio

    print("Rendering a p-mode packground in Minkowski metric")
    sc.Metric = minkowski
    bg=matte_paint(sc, painter, width=nx, height=ny)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    print("Rendering a p-mode packground in KerrBL metric")
    sc.Metric = kerrbl
    bg=matte_paint(sc, painter, width=nx, height=ny)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    print("Rendering a p-mode packground in KerrKS metric")
    sc.Metric = kerrks
    bg=matte_paint(sc, painter, width=nx, height=ny)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


## A simple picture
#
# Here we load it from the web but fp could be the name of a local
# file
if 'Picture' in painters:
    url = ("https://i.pinimg.com/236x/19/53/11/"
           "195311bdb5b029cb72214b95833e6dcf.jpg")
    fp = urllib.request.urlopen(url)
    img = numpy.array(PIL.Image.open(fp))
    painter = Picture(img=img, fov=scr.FieldOfView)
    ny = res
    nx = 3*res//4

    print("Rendering a picture background in Minkowski metric")
    sc.Metric = minkowski
    bg=matte_paint(sc, painter, width=nx, height=ny)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    print("Rendering a picture background in KerrBL metric")
    sc.Metric = kerrbl
    bg=matte_paint(sc, painter, width=nx, height=ny)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


## A full-sky, 360°x180° panorama:
#
# again from he web but fp could be the anme of a local file
if 'Panorama' in painters:
    url = ("http://farm1.staticflickr.com/192/"
           "456185667_adde9d2f8e_o_d.jpg")
    fp = urllib.request.urlopen(url)
    img = numpy.array(PIL.Image.open(fp))
    painter=Panorama(img=img)
    nx=res
    ny=(res*10)//16; # use a 16/10 aspect ratio

    print("Rendering a panoramic background in Minkowski metric" )
    sc.Metric = minkowski
    bg=matte_paint(sc, painter, width=nx, height=ny)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    # This time we pre-compute to play with yaw, pitch an roll:
    print("Pre-computing in KellBL coordinates")
    sc.Metric = kerrbl
    data = sc.rayTrace(width=nx, height=ny, quantities="ImpactCoords")["ImpactCoords"]

    print("Displaying it")
    bg=matte_paint(data, painter, coordkind=sc.Metric.coordKind(),
                   yaw=0., pitch=0., roll=0.)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    print("Displaying it with yaw = 0.1")
    bg=matte_paint(data, painter, coordkind=sc.Metric.coordKind(),
                   yaw=0.1, pitch=0., roll=0.)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    print("Displaying it with yaw and pitch = 0.1")
    bg=matte_paint(data, painter, coordkind=sc.Metric.coordKind(),
                   yaw=0.1, pitch=0.1, roll=0.)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()


    print("Displaying it with yaw, pitch and roll = 0.1")
    bg=matte_paint(data, painter, coordkind=sc.Metric.coordKind(),
                   yaw=0.1, pitch=0.1, roll=0.1)

    plt.imshow(bg)
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()

print("All done, exiting")
if pdf is not None:
    pdf.close()
