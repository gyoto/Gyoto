#!/usr/bin/env yorick
/*
   This script examplifies how to make a movie using Gyoto and LibAV.

   There is an unconspicuous bug in Gyoto by which mutlithreaded
   ray-tracing corrupts the process memory (presumably some block is
   freed in one thread whereas it is still used by the main
   thread). We seem to work around this bug by performing the
   ray-tracing of each image in a separate (multi-threaded) yorick
   process. Therefore, this script is mutli-process (main process
   writes movie file, subprocesses compute images), and the Gyoto
   subprocesses are multi-threaded.

   Required packages:
    gyoto        (gyoto_Scenery)     https://github.com/gyoto/Gyoto
    yorick-av    (av_*)              https://github.com/paumard/yorick-av
    yorick-svipc (ftok, fork, s?m_*) https://github.com/mdcb/yp-svipc
    yutils       (tic, tac)          https://github.com/frigaut/yorick-yutils
   
 */
#include "svipc.i"
#include "gyoto.i"
#include "libav.i"


//// PARAMETERS
ifile="../doc/examples/example-moving-star.xml";
ofile="test.ogg";
vcodec="libtheora";
resolution=64;
t0=1000.;
nframes=1000;
dt=1.;
//gyoto_verbose,0;
//// END PARAMETERS

tic; // Initialize time counter

file_i=current_include();
ppid=getpid();
shmid = ftok(file_i, proj=ppid);
semid = ftok(file_i, proj=ppid+1);

sem_init, semid, nums=2;
shm_init, shmid, slots=1;

// Semaphores are used in this script to lock the shared memory
// segment.  Gyoto write to shared memory and gives 0.  FFmpeg reads
// from shared memory and gives 1.  Initially, FFmpeg gives 1 (nothing
// to read) and takes 0 (waiting for frame).

// Reading palette using plain Yorick requires an X11 window
window;
pli, array(double, 2,2);
palette, query=1, r,g,b;
rgb= transpose([r,g,b]);
top=numberof(r)-1;
winkill;

// Reading and initializing scenery 
sc=gyoto_Scenery(ifile);
ut = sc(metric=)(unitlength=)/GYOTO_C;
sc,quantities="Intensity";
sc,nthreads=2;
noop, sc(screen=)(time=t0*ut);
noop, sc(screen=)(resolution=resolution);
// Specific to star astrobj: precompute orbit
ao = sc(astrobj=);
if (ao(kind=)=="Star") noop, ao(xfill=0.);

sem_give, semid, 1; // Give 1: Gyoto may write to shared memory.
encoder=av_create(ofile, vcodec=vcodec);

for (n=1, t=t0; n<=nframes; ++n, t+=dt) {
  if (fork()) {
    sem_take, semid, 0; // Is there something to read now?
    im = shm_read(shmid, "image");
    if (n==1) cmax=max(im);
    frame = rgb(,bytscl(im(,0:1:-1), cmax=cmax, top=top)+1);
    sem_give, semid, 1; // Reading done, Gyoto may write again.
    write, format="Writting frame nr %d to movie file\n", n;
    av_write, encoder, frame;
  } else {
    write, format="Ray-tracing frame nr %d, time=%e\n", n, t;
    noop, sc(screen=)(time=t*ut);
    im = gyoto_Scenery_rayTrace(sc);
    sem_take, semid, 1; // May we write now?
    shm_write, shmid, "image", &im;
    sem_give, semid, 0; // Image ready to be read.
    quit;
  }
 }
  
av_close, encoder;

sem_cleanup, semid;
shm_cleanup, shmid;

write, format="\nRay-traced %d frames in %g seconds\n", nframes, tac();
