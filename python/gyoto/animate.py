#!/usr/bin/env python3
'''Create movies using Gyoto

This module provides facilities to animate Gyoto Sceneries and render
them as movie files. The module can also be called as a script on the
command line.

rayTrace() takes as argument a user-defined callable func and a
sceneray (among other parameters). func is responsible for mutating
the scenery for each frame. raTrace then ray-traces the scenery. Two
callables suitable as the func argument are provided.

mk_movie() is the top-level function. It opens a video file, calls
rayTrace() in a loop and closes the video.

Classes:
VideoWriter -- abstract interface with video library
OpenCVVideoWriter -- implements VideoWriter using OpenCV-Python
PyAVVideoWriter -- implements VideoWriter using PyAV
NullVideoWriter -- a VideoWriter that does not write anything
orbiting_screen -- a callable that can be used by rayTrace
orbiting_screen_forward -- idem
accelerating_tangential_screen -- idem
growing_mass -- idem

Function:
rayTraceFrame -- mutates Gyoto scenery and raytrace it
defaultScenery -- returns a default Gyoto scenery
defaultTrajectory -- returns a default screen trajectory
mk_video -- make a video

'''
try:
    from . import core, std
except ImportError:
    from gyoto import core, std
import numpy
import os
import matplotlib.pyplot as plt
import matplotlib
import argparse
import warnings

### Helper functions

## An API for video-writing facilities

class VideoWriter:
    '''Generic interface to write videos
    '''

    norm=None
    vmin=0.
    vmax=0.
    cmap=plt.cm.get_cmap('hot')
    
    def __init__(self, filename, fps, width, height):
        '''Initialize video

        Derived classes should open the video in there initializer and
        call the base class __init__.

        '''
        self.fps=fps
        self.width=width
        self.height=height

    def colorize(self, im_float):
        '''Colorize image

        Given a Gyoto frame (numpy array of doubles), returns and RGB
        image (3 planes of uint8 type).

        '''
        if self.norm is None:
            if self.vmax==0:
                self.vmax = numpy.nanmax(im_float)
            if self.vmax != 0.:
                self.norm=matplotlib.colors.Normalize(vmin=0., vmax=self.vmax, clip=True)
                im_float = self.norm(im_float)
        else:
            im_float=self.norm(im_float)
        return numpy.uint8(self.cmap(im_float)*255.)[:, :, :3]

    def write(self, frame):
        '''Write frame to video

        Frame is a numpy RGB image.
        '''
        raise NotImplementedError

    def close(self):
        '''Close video
        '''
        raise NotImplementedError

    def __del__(self):
        self.close()

class NullVideoWriter(VideoWriter):
    '''A VideoWriter that does not write anything

    Mostly useful during preparatory work in cunjunction with
    plot=True.
    '''
    def write(self, frame):
        '''Write frame to video

        Frame is a numpy RGB image.
        '''
        pass

    def close(self):
        '''Close video
        '''
        pass

class OpenCVVideoWriter(VideoWriter):
    '''An implementation of VideoWriter that uses OpenCV-python
    '''
    video=None
    fourcc=0

    def __init__(self, filename, fps, width, height):
        import cv2
        VideoWriter.__init__(self, filename, fps, width, height)
        self.video=cv2.VideoWriter(filename, self.fourcc, fps, (width, height), True)

    def write(self, frame):
        self.video.write(frame[::-1, :, ::-1])

    def close(self):
        if self.video is not None:
            self.video.release()
            self.video=None

class PyAVVideoWriter(VideoWriter):
    '''An implementation of VideoWriter that uses PyAV
    '''
    container=None
    stream=None
    fourcc=0

    def __init__(self, filename, fps, width, height,
                 codec_name='mpeg4', pix_fmt='yuv420p'):
        import av
        self.av=av
        VideoWriter.__init__(self, filename, fps, width, height)
        self.container = av.open(filename, mode='w')
        self.stream = self.container.add_stream(codec_name, rate=fps)
        self.stream.width = width
        self.stream.height = height
        self.stream.pix_fmt = pix_fmt

    def write(self, frame):
        avframe=self.av.VideoFrame.from_ndarray(frame[::-1,:,:], format='rgb24')
        for packet in self.stream.encode(avframe):
            self.container.mux(packet)

    def close(self):
        if self.container is not None:
            for packet in self.stream.encode():
                self.container.mux(packet)
            self.container.close()
            self.stream=None
            self.container=None

## Two types of changes for the screen
#
# The rayTrace function below takes a callable as argument to mutate
# the screen between exposures. We define two such callables (as
# classes) for two kinds of videos below.
#

class accelerating_tangential_screen:
    '''The screen does not move but has increasing velocity
    
    Members:
    maxvel -- norm of velocity for last frame. Should never reach 1!
    '''

    maxvel=0.99

    def __init__(self, **args):
        for key in args:
            setattr(self, key, args[key])

    def __call__(self, sc, k, nframes):
        scr=sc.screen()
        metric=sc.metric()
        pos=scr.getObserverPos()
        ucirc=metric.circularVelocity(pos)
        uzamo=metric.zamoVelocity(pos)
        Gamma=-metric.ScalarProd(pos, ucirc, uzamo)
        Vzamo=ucirc/Gamma-uzamo           # ucirc as seen by zamo
        norm_circ=metric.ScalarProd(pos, Vzamo, Vzamo)
        norm_wanted=self.maxvel*k/(nframes-1)
        Vzamo *= norm_wanted/numpy.sqrt(norm_circ)    # rescale velocity
        Gamma2 = 1./(1.-metric.ScalarProd(pos, Vzamo, Vzamo))
        assert Gamma2 >= 0, 'Gamma2 < 0! VzamoxVzamo='+str(metric.ScalarProd(pos, Vzamo, Vzamo))+', norm_wanted='+str(norm_wanted)+', norm_circ='+str(norm_circ)
        Gamma=numpy.sqrt(Gamma2)
        fourvel=Gamma*(uzamo+Vzamo)
        scr.fourVel(fourvel)

class static_screen:
    '''The screen does not move

    This is meaningful when the astrobj itself is changing. The camera
    will typically be far away, with ObserverKind set to
    ObserverAtInfinity.

    members:
    tref -- reference time, read from sc.screen() the first time
          __call__ is called
    t0 -- coordinate time at movie start
    t1 -- coordinate time at movie end
    unit -- unit in which t0 and t1 are expressed
    '''
    tref=None
    t0=0.
    t1=1000.
    unit='geometric'

    def __init__(self, **args):
        for key in args:
            setattr(self, key, args[key])

    def __call__(self, sc, k, nframes):
        '''update time in screen
        '''
        if self.tref is None:
            self.tref=sc.screen().time(self.unit)
        t = self.t0+k*(self.t1-self.t0)/(nframes-1)
        sc.screen().time(self.tref+t, self.unit)

class orbiting_screen:
    '''The screen follows an orbit, camera looks along -er

    members:
    t0 -- proper time at movie start
    t1 -- proper time at movie end
    trajectory -- a gyoto.std.Star (or anything else that provides
          getCoord(t, coord, proper).
    '''
    t0=0.
    t1=1000.
    trajectory=None

    def __init__(self, **args):
        for key in args:
            setattr(self, key, args[key])

    def __call__(self, sc, k, nframes):
        t = self.t0+k*(self.t1-self.t0)/(nframes-1)
        coord=core.vector_double(8)
        self.trajectory.getCoord(t, coord, True)
        pos=[coord[i] for i in range(4)]
        vel=[coord[i] for i in range(4, 8)]
        screen=sc.screen()
        screen.setObserverPos(pos)
        screen.fourVel(vel)

class orbiting_screen_forward:
    '''The screen follows an orbit, camera looks forward

    members:
    t0 -- proper time at movie start
    t1 -- proper time at movie end
    roll -- roll angle: left is rotated roll degrees towards up
    trajectory -- a gyoto.std.Star (or anything else that provides
          getCoord(t, coord, proper).
    '''
    t0=0.
    t1=1000.
    roll=0.
    trajectory=None

    def __init__(self, **args):
        for key in args:
            setattr(self, key, args[key])
        self.metric=self.trajectory.metric()

    def __call__(self, sc, k, nframes):
        metric=sc.metric()
        tau = self.t0+k*(self.t1-self.t0)/(nframes-1)
        coord=core.vector_double(8)
        self.trajectory.getCoord(tau, coord, True)
        pos=numpy.asarray([coord[i] for i in range(4)])
        vel=numpy.asarray([coord[i] for i in range(4, 8)])
        screen=sc.screen()
        screen.setObserverPos(pos)
        screen.fourVel(vel)

        # The 3 other elements of the tetrad are initialized as
        # follows: front is along the 3-velocity; up is er projected
        # on the plane orthogonal to front; left is the external
        # product up x front, thus (left, up, front) is direct.  Then
        # we orthonormalise this tetrad. Obviously, this will fail if
        # vel3c is along er.
        #
        # We convert everything to cartesian to compute the external
        # product.
        front=numpy.zeros(4)
        metric.cartesianVelocity(coord, front[1:])
        front /= numpy.sqrt((front*front).sum())
        if metric.coordKind()==core.GYOTO_COORDKIND_SPHERICAL:
            t=pos[0]
            r=pos[1]
            theta=pos[2]
            phi=pos[3]
            st=numpy.sin(theta)
            ct=numpy.cos(theta)
            sp=numpy.sin(phi)
            cp=numpy.cos(phi)
            posr=[pos[0], r*st*cp, r*st*sp, r*ct]
        else:
            posr=pos
            r=numpy.sqrt(pos[1:]**2).sum()
        up0=numpy.concatenate(([0.], posr[1:]/r))
        up0 -= (front*up0).sum()*front
        up0 /= numpy.sqrt((up0*up0).sum())
        left0=numpy.zeros(4)
        left0[1]=up0[2]*front[3]-up0[3]*front[2]
        left0[2]=up0[3]*front[1]-up0[1]*front[3]
        left0[3]=up0[1]*front[2]-up0[2]*front[1]
        rollr=self.roll*numpy.pi/180
        cr=numpy.cos(rollr)
        sr=numpy.sin(rollr)
        up=cr*up0-sr*left0
        left=cr*left0+sr*up0

        if metric.coordKind()==core.GYOTO_COORDKIND_SPHERICAL:
            er=posr/r
            ephi=numpy.asarray([0., -sp, cp, 0.])
            etheta=numpy.asarray([0., ct*cp, ct*sp, -st])
            front=numpy.asarray([
                0.,
                (er*front).sum(),
                (etheta*front).sum(),
                (ephi*front).sum()
            ])
            up=numpy.asarray([
                0.,
                (er*up).sum(),
                (etheta*up).sum(),
                (ephi*up).sum()
            ])
            left=numpy.asarray([
                0.,
                (er*left).sum(),
                (etheta*left).sum(),
                (ephi*left).sum()
            ])
        metric.GramSchmidt(pos, vel, up, front, left);
        screen.screenVector1(left)
        screen.screenVector2(up)
        screen.screenVector3(front)

class growing_mass:
    '''The mass of the central object changes

    The Astrobj needs to be PatternDisk-like (it needs to have the
    properties InnerRadius and OuterRadius).

    '''
    delta0=0.01
    deltaMax0=1.
    rin0=0.
    rout0=28.
    rmax0=50.
    d0=28.
    factor_first=100.
    factor_last=2.01/28.

    def __init__(self, scenery=None, **args):
        if scenery is not None:
            self.delta0=scenery.delta()
            scr=scenery.screen()
            self.d0=scr.distance('geometrical')
            ao=scenery.astrobj()
            self.deltaMax0=ao.deltaMaxInsideRMax()
            self.rin0=ao.get('InnerRadius', 'geometrical')
            self.rout0=ao.get('OuterRadius', 'geometrical')
            self.rmax0=ao.rMax()
            self.factor_last=2.01/self.d0
            scenery.screen().observerKind('ZAMO')
        for key in args:
            setattr(self, key, args[key])

    def __call__(self, sc, k, nframes):
        # factor=(((self.factor_last*k)
        #          +(self.factor_first*(self.nframes-1-k)))
        #         /(self.nframes-1))
        # print(factor)
        log2_first=numpy.log2(self.factor_first)
        log2_last=numpy.log2(self.factor_last)
        log2_k=(((log2_last*k)
                 +(log2_first*(nframes-1-k)))
                /(nframes-1))
        factor=2.**log2_k
        sc.delta(factor*self.delta0)
        scr=sc.screen()
        scr.distance(factor*self.d0, 'geometrical')
        print(factor, scr.distance('geometrical'))
        ao=sc.astrobj()
        ao.deltaMaxInsideRMax(factor*self.deltaMax0)
        ao.set('InnerRadius', factor*self.rin0, 'geometrical')
        ao.set('OuterRadius', factor*self.rout0, 'geometrical')
        ao.rMax(factor*self.rmax0)

## Helper function for tracing one frame

def rayTraceFrame(sc, func, k, nframes, width, height):
    '''Ray-trace one frame of a Gyoto video

    Parameters:
      sc:     Scenery to ray-trace
      func:   callable to mutate the Screen
      k:      number of the frame
      width:  width of the video
      height: height of the video

    Returns:
      The raytraced intensity as a NumPy array
    '''
    intensity=numpy.zeros((height, width))
    pintensity=core.array_double.fromnumpy2(intensity)
    func(sc, k, nframes)
    res=max(width, height)
    sc.screen().resolution(res)
    ii=core.Range(res//2-width//2+1, res//2-width//2+width, 1)
    jj=core.Range(res//2-height//2+1, res//2-height//2+height, 1)
    grid=core.Grid(ii, jj)
    aop=core.AstrobjProperties()
    aop.intensity=pintensity
    sc.rayTrace(grid, aop)
    # print(newpos)
    # print(newvel)
    # plt.imshow(intensity)
    # plt.show()
    return intensity
    
## Build default scenery and trajectory

def defaultScenery():
    '''Create a default scenery
    The astrobj is a PatternDisk.
    '''
    metric = core.Metric("KerrBL")
    metric.mass(4e6, "sunmass");
    gridshape=numpy.asarray( (1, 3, 11) , numpy.uint64)
    pgridshape=core.array_size_t.fromnumpy1(gridshape)
    opacity=numpy.zeros(gridshape)
    popacity=core.array_double.fromnumpy3(opacity)
    opacity[:, 0::2, 0::2]=100.
    opacity[:, 1::2, 1::2]=100.
    intensity=opacity*0.+1.;
    pintensity=core.array_double.fromnumpy3(intensity)
    pd=std.PatternDisk()
    pd.velocityKind('ZAMO')
    pd.copyIntensity(pintensity, pgridshape)
    pd.copyOpacity  (popacity, pgridshape)
    pd.innerRadius(0)
    pd.outerRadius(28)
    pd.repeatPhi(8)
    pd.metric(metric)
    pd.rMax(50)
    screen=core.Screen()
    screen.metric(metric)
    screen.resolution(64)
    screen.time(1000., "geometrical_time")
    screen.distance(28., "geometrical")
    # Standard 24x36 field of view after a 55mm objective
    screen.fieldOfView(2.*numpy.arctan(18./55), 'radians')
    screen.anglekind('Rectilinear')
    screen.inclination(95., "degree")
    screen.PALN(180., "degree")
    sc=core.Scenery()
    sc.metric(metric)
    sc.screen(screen)
    sc.astrobj(pd)
    sc.nThreads(8)
    return sc

def defaultTrajectory(screen):
    '''Get default trajectory and adapt screen (field-of-view etc.)
    '''
    screen.observerKind('VelocitySpecified')
    screen.fieldOfView(90, 'degree')
    traj=std.Star()
    traj.metric(screen.metric())
    traj.setInitCoord((0., 28., 0.8, 0.), (0., 0., 0.007))
    return traj

### The main function

def mk_video(scenery=None,
             func="orbiting_screen",
             orbit_trajectory=None, orbit_t0=0., orbit_t1=1000.,
             static_t0=0., static_t1=1000., static_unit='geometrical',
             acceleration_maxvel=0.99,
             growth_factor_first=None,
             growth_factor_last=None,
             duration=10, fps=3, width=128, height=72,
             dangle1=None, dangle2=None, fov=None,
             verbose=0, debug=False, nthreads=8,
             
             output=None,
             backend=OpenCVVideoWriter,
             cmap=None,
             observerkind=None,
             plot=False,
             frame_first=0,
             frame_last=None
             ):
    '''Make a video from a Gyoto Scenery

    Keyword arguments:

    scenery  -- the Scenery to animate, a string (XML file name),
                Gyoto.core.Scenery instance or None in which case a
                default Scenery is used
    func     -- a callable which will be called as func(scenery, k)
                where k is the frame number. func() is responsible to
                mutate the Scenery for each frame, for instance by
                moving the camera, changing its field-of view etc..
                func may also be a string, the name of one of the
                built-in callables: orbiting_screen (default),
                orbiting_screen_forward, static_screen,
                accelerating_tangential_screen or growing_mass. Those
                options take parameters below
    orbit_trajectory -- 
                only if func above starts with 'orbiting_screen', a
                gyoto.std.Star instance or a string (XML file name) or
                None (in which case a default trajectory is used)
                specifying the trajectory of the Screen.
    orbit_t0 -- only if func above starts with 'orbiting_screen',
                proper time along the trajectory at the beginning of
                the movie
    orbit_t1 -- only if func above starts with 'orbiting_screen',
                proper time along the trajectory at the end of the
                movie
    static_t0-- only if func above is 'static_screen', coordinate
                observing time offset at first frame
    static_t0-- only if func above is 'static_screen', coordinate
                observing time offset at last frame
    static_unit --
                only if func above is 'static_screen', unit in which
                static_t0 and static_t1 are expressed ('geometrical').
    acceleration_maxvel --
                only if func is 'accelerating_tangential_screen',
                maximum velocity of the observer, in terms of c
                (default: 0.99)
    growth_factor_first --
                only if func is 'growing_mass', scale factor for first
                frame.
    growth_factor_last --
                only if func is 'growing_mass', scale factor for last
                frame.
    duration -- duration pf the movie in seconds (default: 10.)
    fps      -- frames per second (default: 3)
    width    -- width of the movie in pixels (default: 128)
    height   -- height of the movie in pixels (default: 72)
    dangle1  -- rotate the camera dangle1 degrees horizontally
    dangle2  -- rotate the camera dangle2 degrees vertically
    fov      -- screen field-of-view
    verbose  -- verbosity level (default: 0)
    debug    -- debug mode (default: False)
    nthreads -- number of parallel threads to use (default: 8)
    output   -- output file name
    backend  -- class implementing VideoWriter
    observerkind --
                change observer kind (default: VelocitySpecified)
    plot     -- show each frame (default: False)

    '''

    # Set a few variables
    nframes=int(duration*fps)

    # Read or create Scenery
    if scenery is None:
        sc=defaultScenery()
    elif type(scenery) is core.Scenery:
        sc=scenery
    else:
        sc=core.Factory(scenery).scenery()
    screen=sc.screen()
    metric=sc.metric()

    # Select video type, init func callable
    if func == 'orbiting_screen':
        # Read or create trajectory
        if orbit_trajectory is None:
            traj=defaultTrajectory(screen)
        elif type(orbit_trajectory) is std.Star:
            traj=orbit_trajectory
        else:
            traj=std.Star(core.Factory(orbit_trajectory).astrobj())
            traj.metric(metric)
        screen.observerKind('VelocitySpecified')
        func=orbiting_screen(trajectory=traj, t0=orbit_t0, t1=orbit_t1)
    elif func == 'static_screen':
        func=static_screen(unit=static_unit, t0=static_t0, t1=static_t1)
    elif func == 'orbiting_screen_forward':
        # Read or create trajectory
        if orbit_trajectory is None:
            traj=defaultTrajectory(screen)
        elif type(orbit_trajectory) is std.Star:
            traj=orbit_trajectory
        else:
            traj=std.Star(core.Factory(orbit_trajectory).astrobj())
            traj.metric(metric)
        screen.observerKind('FullySpecified')
        func=orbiting_screen_forward(trajectory=traj, t0=orbit_t0, t1=orbit_t1)
    elif func == 'accelerating_tangential_screen':
        screen.observerKind('VelocitySpecified')
        screen.dangle1(-45, 'degree')
        screen.fieldOfView(90, 'degree')
        func=accelerating_tangential_screen(maxvel=acceleration_maxvel)
    elif func == 'growing_mass':
        func=growing_mass(sc)
        if growth_factor_first is not None:
            func.factor_first=growth_factor_first
        if growth_factor_last is not None:
            func.factor_last=growth_factor_last
    # else assume func is a callable

    # Override some values set on command-line
    if (observerkind is not None):
        screen.observerKind(observerkind)
    if dangle1 is not None:
        screen.dangle1(dangle1, 'degree')
    if dangle2 is not None:
        screen.dangle2(dangle2, 'degree')
    if fov is not None:
        screen.fieldOfView(fov, 'degree')

    # Prepare for ray-tracing
    nframes=int(duration*fps)
    screen.resolution(max(height, width))
    core.verbose(verbose)
    core.debug(debug)
    sc.nThreads(nthreads)

    # Open video
    if backend == 'OpenCV':
        backend=OpenCVVideoWriter
    elif backend == 'PyAV':
        backend=PyAVVideoWriter
    elif backend is None or backend == 'Null':
        backend=NullVideoWriter

    if type(output) == str or backend is NullVideoWriter:
        video=backend(output, fps, width, height)
    elif isinstance(output, VideoWriter):
        video=output
    else:
        raise ValueError('output needs to be a string or VideoWriter')
    if type(cmap) == str:
        video.cmap=plt.cm.get_cmap(cmap)
    elif cmap is not None:
        video.cmap=cmap

    # Loop on frame number
    if frame_last is None:
        frame_last=nframes-1
    for k in range(frame_first, frame_last+1):
        print(k, "/", nframes)
        intensity=rayTraceFrame(sc, func, k, nframes, width, height)
        frame=video.colorize(intensity)
        if plot:
            plt.imshow(frame, origin='lower')
            plt.show()
        video.write(frame)

    # Close video
    del video

### Define arguments

parser = argparse.ArgumentParser(description=
'''Make a video from a Gyoto Scenery

Synopsys:
gyoto mk-video --output=<out.avi> [--options...]

Several types of videos are available and can be selected with
--func. Currently, only the KerrBL Metric implements the methods
needed for those various videos (zamoVelocity, circularVelocity and
observerTetrad are needed). All parameters except --output have
reasonable defaults. For even more advanced usage, see the module
'gyoto.animate' in Python. This script requires one of the modules
OpenCV-python or PyAV to be installed (see --backend option).

gyoto mk-video --output=<out.avi> \\
               --func=static_screen \\
              [--scenery=<scenery.xml>] \\
              [--static-unit=<time unit>] \\
              [--static-t0=<t0>] \\
              [--static-t1=<t1>]

will produce a video of the Scenery <scenery.xml> with the Screen
(a.k.a. camera) fixed, with observing (coordinate) time varying from
Time+<t0> to Time+<t1> in <time unit>, where Time is the observing
date specified in<scenery.xml>. This kind of video makes sense only if
the astrobj in <scenery.xml> is time-variable.

gyoto mk-video --output=<out.avi> \\
               --func=orbiting_screen \\
              [--scenery=<scenery.xml>] \\
              [--orbit-trajectory=<orbit.xml>] \\
              [--orbit-t0=<t0>] \\
              [--orbit-t1=<t1>]

will produce a video of the Scenery <scenery.xml> with the Screen
(a.k.a. camera) orbiting the central object along the trajectory
described in <trajectory.xml>. <trajectory.xml> may be a scenery file
or a file produced by gyotoy and must contain a Star astrobj. The
geodesic of this Star will be used as the camera trajectory. The
Screen is pointing along -er, i.e. in the direction of the central object.


gyoto mk-video --output=<out.avi> \\
               --func=orbiting_screen_forward \\
              [--scenery=<scenery.xml>] \\
              [--orbit-trajectory=<orbit.xml>] \\
              [--orbit-t0=<t0>] \\
              [--orbit-t1=<t1>]

will produce a video of the Scenery <scenery.xml> with the Screen
(a.k.a. camera) orbiting the central object along the trajectory
described in <trajectory.xml>. <trajectory.xml> may be a scenery file
or a file produced by gyotoy and must contain a Star astrobj. The
geodesic of this Star will be used as the camera trajectory. The
Screen is pointing in the direction of motion, i.e. looking forward.


gyoto mk-video --output=<out.avi> \\
               --func=accelerating_tangential_screen \\
              [--scenery=<scenery.xml>] \\
              [--acceleration-maxvel=<maxvel>]

will produce a video where the Screen is fixed but changes velocity
from 0 to maxvel*c. This is not an actual time sequence but helps
vizualizing effects such as light aberration.

gyoto mk-video --output=<out.avi> \\
               --func=growing_mass \\
              [--scenery=<scenery.xml>] \\
              [--growth_factor_first=<ff>] \\
              [--growth_factor_last=<fl>]

will produce a video where the scale of the Scenery changes for each
frame, from ff to fl. For this particular type of video, the Astrobj
in the Scenery must be some kind of PatternDisk.


''',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 prefix_chars='-+')
parser.add_argument('-s', '--scenery', type=str, default=None,
                    help='name of a Gyoto Scenery XML file. Default: use built-in Scenery.')
parser.add_argument('-t', '--orbit-trajectory', type=str, default=None,
                    dest='orbit_trajectory',
                    help='name of a Gyoto Scenery or Astrobj XML file containing a star '
                    'describing the screen motion. Default: use built-in trajectory.')
parser.add_argument('-o', '--output', type=str, default=None,
                    help='name of video file to save the movie in')
parser.add_argument('-B', '--backend', type=str, default='OpenCV',
                    choices=['OpenCV', 'PyAV', 'Null'],
                    help='name of backend to create video')
parser.add_argument('-c', '--cmap', type=str, default='hot',
                    help='name of pyplot color map')
parser.add_argument("-V", "--func", help="type of video to produce.",
                    type=str, default='orbiting_screen',
                    choices=['orbiting_screen', 'orbiting_screen_forward',
                             'accelerating_tangential_screen',
                             'growing_mass', 'static_screen'])
parser.add_argument("-D", "--duration", help="movie duration in seconds",
                    type=float, default=10.)
parser.add_argument("-f", "--fps", help="number of frames per second",
                    type=int, default=3)
parser.add_argument("-W", "--width", help="image width in pixels",
                    type=int, default=128)
parser.add_argument("-H", "--height", help="image height in pixels",
                    type=int, default=72)
parser.add_argument("-a", "--dangle1", help="camera azimuth offset in degrees",
                    type=float, default=None)
parser.add_argument("-b", "--dangle2", help="camera elevation offset in degrees",
                    type=float, default=None)
parser.add_argument("-F", "--fov", help="camera field-of-view in degrees",
                    type=float, default=None)
parser.add_argument("-O", "--observerkind", help="observer kind",
                    type=str, default=None)
parser.add_argument("-v", "--verbose", help="verbosity level",
                    type=int, default=0)
parser.add_argument("-d", "--debug", help="debug mode",
                    dest='debug', action='store_true', default=False)
parser.add_argument("-p", "--plot", help="plot each frame",
                    dest='plot', action='store_true', default=False)
parser.add_argument("-T", "--nthreads", help="number of threads to use",
                    type=int, default=8)
parser.add_argument("--orbit-t0", help="for orbit video type, initial time in geometrical units",
                    dest='orbit_t0', type=float, default=0)
parser.add_argument("--orbit-t1", help="for orbit video type, final time in geometrical units",
                    dest='orbit_t1', type=float, default=1000)
parser.add_argument("--static-t0", help="for static video type, initial time offset",
                    dest='static_t0', type=float, default=0)
parser.add_argument("--static-t1", help="for static video type, final time offset",
                    dest='static_t1', type=float, default=1000)
parser.add_argument("--static-unit", dest="static_unit", help="unit for static-t0 and static-t1",
                    type=str, default="geometrical")
parser.add_argument("--acceleration-maxvel", help="for acceleration video type, max velocity in terms of light velocity",
                    dest='acceleration_maxvel', type=float, default=0.99)
parser.add_argument("--growth_factor_first", help="for growth video type, scale factor on first frame",
                    dest='growth_factor_first', type=float, default=None)
parser.add_argument("--growth-factor-last", help="for growth video type, scale factor on last frame",
                    dest='growth_factor_last', type=float, default=None)

def main():
    args = parser.parse_args()
    mk_video(**args.__dict__)

# If called as script, process command line and produce video
if (__name__ == "__main__"):
    main()
