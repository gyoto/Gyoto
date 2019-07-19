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
orbiting_screen -- a callable that can be used by rayTrace
accelerating_tangential_screen -- idem

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
                self.vmax = im_float.max()
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

## Two types of changes for the screen
#
# The rayTrace function below takes a callable as argument to mutate
# the screen between exposures. We define two such callables (as
# classes) for two kinds of videos below.
#

class accelerating_tangential_screen:
    '''The screen does not move but has increasing velocity
    
    Members:
    dvel -- norm of velocity will be k*dvel. Should never reach 1!
    '''
    dvel=1e-3

    def __init__(self, **args):
        for key in args:
            setattr(self, key, args[key])

    def __call__(self, sc, k):
        pos=numpy.zeros(4)
        scr=sc.screen()
        metric=sc.metric()
        scr.getObserverPos(pos)
        ucirc=metric.circularVelocity(pos)
        uzamo=metric.zamoVelocity(pos)
        Gamma=-metric.ScalarProd(pos, ucirc, uzamo)
        Vzamo=ucirc/Gamma-uzamo           # ucirc as seen by zamo
        norm_circ=metric.ScalarProd(pos, Vzamo, Vzamo)
        norm_wanted=k*self.scale
        Vzamo *= norm_wanted/numpy.sqrt(norm_circ)    # rescale velocity
        Gamma2 = 1./(1.-metric.ScalarProd(pos, Vzamo, Vzamo))
        assert Gamma2 >= 0, 'Gamma2 < 0! VzamoxVzamo='+str(metric.ScalarProd(pos, Vzamo, Vzamo))+', norm_wanted='+str(norm_wanted)+', norm_circ='+str(norm_circ)
        Gamma=numpy.sqrt(Gamma2)
        fourvel=Gamma*(uzamo+Vzamo)
        scr.fourVel(fourvel)

class orbiting_screen:
    '''The screen follows an orbit
    '''
    t0=0.
    dt=1.
    trajectory=None

    def __init__(self, **args):
        for key in args:
            setattr(self, key, args[key])

    def __call__(self, sc, k):
        t = self.t0+k*self.dt
        coord=core.vector_double(8)
        self.trajectory.getCoord(t, coord, True)
        pos=[coord[i] for i in range(4)]
        vel=[coord[i] for i in range(4, 8)]
        screen=sc.screen()
        screen.setObserverPos(pos)
        screen.fourVel(vel)

## Helper function for tracing one frame

def rayTraceFrame(sc, func, k, width, height):
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
    pintensity=core.array_double_fromnumpy2(intensity)
    func(sc, k)
    res=sc.screen().resolution()
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
    pgridshape=core.array_size_t_fromnumpy1(gridshape)
    opacity=numpy.zeros(gridshape)
    popacity=core.array_double_fromnumpy3(opacity)
    opacity[:, 0::2, 0::2]=100.
    opacity[:, 1::2, 1::2]=100.
    intensity=opacity*0.+1.;
    pintensity=core.array_double_fromnumpy3(intensity)
    pd=std.PatternDisk()
    pd.copyIntensity(pintensity, pgridshape)
    pd.copyOpacity  (popacity, pgridshape)
    pd.innerRadius(3)
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
    screen.fieldOfView(2.*numpy.arctan(18./55), 'degree')
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
             acceleration_maxvel=0.99,
             duration=10, fps=3, width=128, height=72,
             dangle1=None, dangle2=None, fov=None,
             verbose=0, debug=False, nthreads=8,
             
             output=None,
             backend=OpenCVVideoWriter,
             observerkind=None,
             plot=False,
             ):
    '''Make video from a Gyoto Scenery

    Keyword arguments:

    scenery  -- the Scenery to animate, a string (XML file name),
                Gyoto.core.Scenery instance or None in which case a
                default Scenery is used
    func     -- a callable which will be called as func(scenery, k)
                where k is the frame number. func() is responsible to
                mutate the Scenery for each frame, for instance by
                moving the camera, changing its field-of view etc..
                func may also be a string, the name of one of the
                built-in callables: orbiting_screen (default) or
                accelerating_tangential_screen. Those options take
                parameters below
    orbit_trajectory -- 
                only if func above is 'orbiting_screen', a
                gyoto.std.Star instance or a string (XML file name) or
                None (in which case a default trajectory is used)
                specifying the trajectory of the Screen.
    orbit_t0 -- only if func above is 'orbiting_screen', time
                coordinate of along the trajectory at the beginning of
                the movie
    orbit_t1 -- only if func above is 'orbiting_screen', time
                coordinate of along the trajectory at the end of the
                movie
    acceleration_maxvel --
                only if func is 'accelerating_tangential_screen',
                maximum velocity of the observer, in terms of c
                (default: 0.99)
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
        screen.observerKind('VelocitySpecified')
        # Read or create trajectory
        if orbit_trajectory is None:
            traj=defaultTrajectory(screen)
        elif type(orbit_trajectory) is std.Star:
            traj=orbit_trajectory
        else:
            traj=std.Star(core.Factory(orbit_trajectory).astrobj())
            traj.metric(metric)
        func=orbiting_screen(trajectory=traj, t0=orbit_t0, dt=orbit_t1/(nframes-1))
    elif func == 'accelerating_tangential_screen':
        screen.observerKind('VelocitySpecified')
        screen.dangle1(-45, 'degree')
        screen.fieldOfView(90, 'degree')
        func=accelerating_tangential_screen(scale=acceleration_maxvel/(nframes-1))
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
    video=backend(output, fps, width, height)

    # Loop on frame number
    for k in range(nframes):
        print(k, "/", nframes)
        intensity=rayTraceFrame(sc, func, k, width, height)
        frame=video.colorize(intensity)
        if plot:
            plt.imshow(frame)
            plt.show()
        video.write(frame)

    # Close video
    del video

### Define arguments

parser = argparse.ArgumentParser(description=__doc__,
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
parser.add_argument("-V", "--func", help="type of video to produce.",
                    type=str, default='orbiting_screen',
                    choices=['orbiting_screen', 'accelerating_tangential_screen'])
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
parser.add_argument("--acceleration-maxvel", help="for acceleration video type, max velocity in terms of light velocity",
                    dest='acceleration_maxvel', type=float, default=0.99)
parser.add_argument("--orbit-t0", help="for orbit video type, initial time in geometrical units",
                    dest='orbit_t0', type=float, default=0)
parser.add_argument("--orbit-t1", help="for orbit video type, final time in geometrical units",
                    dest='orbit_t1', type=float, default=1000)

# If called as script, process command line and produce video
if (__name__ == "__main__"):
    args = parser.parse_args()
    mk_video(**args.__dict__)
