#!/usr/bin/env yorick
/*
   This script examplifies how to make a movie using Gyoto and LibAV.

   As of writing there is a bug in Gyoto which makes this script
   occasionaly fail with a segmentation fault or message such as
yorick(20187,0x7fff706cdcc0) malloc: *** error for object 0x10038d0c8: incorrect checksum for freed object - object was probably modified after being freed.

   mk-movie-fork-raytracing.i tries to work around this bug by
   insulating each call to gyoto_Scenery_rayTrace in its own process.

   Required packages:
    gyoto        (gyoto_Scenery)
    yorick-av    (av_*)
    yutils       (tic, tac)
 */

#include "gyoto.i"
#include "pnm.i"
#include "libav.i"

ifile=GYOTO_EXAMPLES_DIR+"example-moving-star.xml";
ofile="test_single.ogg";
resolution=64;
// use jmin and jmax to fix aspect ratio
jcrop = resolution/8; // aspect ratio 4:3
jmin=1+jcrop;
jmax=resolution-jcrop
nthreads=8;
gyoto_verbose,0;

tic;

t0=1000.;
nframes=2000;
dt=1.;

sc=gyoto_Scenery(ifile);
ut = sc(metric=)(unitlength=)/GYOTO_C;
sc,quantities="Intensity";
sc,nthreads=nthreads;
noop, sc(screen=)(time=t0*ut);
noop, sc(screen=)(resolution=resolution);

// Specific to star astrobj: precompute orbit
ao = sc(astrobj=);
if (ao(kind=)=="Star") noop, ao(xfill=0.);

// Change star size to something reasonably physical, yet nice
ao, radius=1, opticallythin=0;

// Precompute mask, if possible for this astrobj
stt = ao(startrace= t0, t0+nframes*dt)(opticallythin=0, radius=2.*ao(radius=), delta=0.5*ao(radius=));
sc, astrobj=stt;
mask=sc(,,"Intensity");
sc, astrobj=ao;
noop, sc.screen(mask=mask);
// Check the mask at least once!
//pli, mask(,jmin:jmax); limits, square=1; pause, 10000; winkill;

encoder=av_create(ofile);

window;
pli, array(double, 2,2);
palette, query=1, r,g,b;
rgb= transpose([r,g,b]);
top=numberof(r)-1;
winkill;

for (n=1, t=t0; n<=nframes; ++n, t+=dt) {
  write, format="\nRay-tracing frame nr %d of %d, observing time=%e\n",
    n, nframes, t;
  noop, sc(screen=)(time=t*ut);
  im = sc(,jmin:jmax,"Intensity");
  if (n==1) cmax=max(im);
  frame = rgb(,bytscl(im(,0:1:-1), cmax=cmax, top=top)+1);
  av_write, encoder, frame;
 }

av_close, encoder;
write, format="\nRay-traced %d frames in %g seconds\n", nframes, tac();
