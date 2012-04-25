#!/usr/bin/env yorick
/*
   This script examplifies how to make a movie using Gyoto and LibAV.

   There is an unconspicuous bug in Gyoto by which mutlithreaded
   ray-tracing corrupts the process memory (presumably some block is
   freed in one thread whereas it is still used by the main
   thread). We seem to work around this bug by performing the
   ray-tracing of each image in a separate (multi-threaded) yorick
   process.

   Required packages:
    gyoto        (gyoto_Scenery)
    yorick-av    (av_*)
    yorick-svipc (ftok, fork, shm_*, sem_*)
    yutils       (tic, tac)
   
 */
#include "svipc.i"
#include "gyoto.i"
#include "libav.i"

ifile="../doc/examples/example-moving-star.xml";
ofile="test.ogg";
vcodec="libtheora";
resolution=64;
//gyoto_verbose,0;

tic;

t0=1000.;
nframes=200;
dt=1.;

file_i=current_include();
ppid=getpid();
shmid = ftok(file_i, proj=ppid);
semid = ftok(file_i, proj=ppid+1);

sem_init, semid, nums=3;
shm_init, shmid, slots=1;

window;
pli, array(double, 2,2);
palette, query=1, r,g,b;
rgb= transpose([r,g,b]);
top=numberof(r)-1;
winkill;

sc=gyoto_Scenery(ifile);
ut = sc(metric=)(unitlength=)/GYOTO_C;
sc,quantities="Intensity";
sc,nthreads=2;
noop, sc(screen=)(time=t0*ut);
noop, sc(screen=)(resolution=resolution);

sem_give, semid, 1; // tell gyoto process we are ready
encoder=av_create(ofile, vcodec=vcodec);

for (n=1, t=t0; n<=nframes; ++n, t+=dt) {
  if (fork()) {
    sem_take, semid, 0; // wait for gyoto thread to write image
    im = shm_read(shmid, "image");
    if (n==1) cmax=max(im);
    frame = rgb(,bytscl(im(,0:1:-1), cmax=cmax, top=top)+1);
    sem_give, semid, 1; // tell gyoto we are ready
    av_write, encoder, frame;
  } else {
    write, format="frame nr %d, time=%e\n", n, t;
    noop, sc(screen=)(time=t*ut);
    im = gyoto_Scenery_rayTrace(sc);
    sem_take, semid, 1; // the FFmpeg process gives when ready to read
    shm_write, shmid, "image", &im;
    sem_give, semid, 0; // tell FFmpeg process to read
    quit;
  }
 }
  
av_close, encoder;

sem_cleanup, semid;
shm_cleanup, shmid;

write, format="\nRay-traced %d frames in %g seconds\n", nframes, tac();
