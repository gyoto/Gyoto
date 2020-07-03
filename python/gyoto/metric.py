'''Gyoto::Metric namespace

In order to emulate the C++ Gyoto::Metric namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Metrics
in here.

'''

import gyoto._namespaces as _namespaces
from gyoto.core import Metric as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces
Complex=ComplexMetric

import gyoto.core
import numpy

def jacobian_numerical(metric, pos, epsilon=1e-6):
    '''Estimate the Jacobian matrix of a metric numerically

    This function is intended for debugging using. For production,
    using the method gyoto.core.Metric.jacobian is preferred.

    If `metric' is an instance of a subclass of gyoto.core.Metric,
    jacobian_numerical(metric, pos) should yield the same as
    metric.jacobian(pos) within numerical errors.

    Keyword arguments:
    metric -- the gyoto.core.Metric instance to work on
    pos -- the coordinates at which to estimate the Christoffel symbols
    epsilon -- the step for estimating of the derivatives (default 1e-6)

    '''
    delta=numpy.empty((4, 4, 4))
    posa = numpy.asarray(pos)
    posb = posa.copy()
    ga = metric.gmunu(posa)
    for alpha in range(4):
        posb[alpha] += epsilon
        gb =  metric.gmunu(posb)
        delta[alpha, :, :] = (gb-ga)/epsilon
        posb[alpha]=posa[alpha]
    return delta

def christoffel_numerical(metric, pos, epsilon=1e-6):
    '''Estimate the Christoffel symbols of a metric numerically

    This function is intended for debugging using. It is called by
    gyoto.metric.check_christoffel for this purpose. For production,
    using the method gyoto.core.Metric.christoffel is preferred.

    If `metric' is an instance of a subclass of gyoto.core.Metric,
    christoffel_numerical(metric, pos) should yield the same as
    metric.christoffel(pos) within numerical errors.

    This function estimates the Christoffel symbols by estimating
    numerically the partial derivatives of the metric coefficients
    (given by metric.gmunu(pos)) and inverting (also numerically) the
    covariant metric coefficients matrix as pos.

    Keyword arguments:
    metric -- the gyoto.core.Metric instance to work on
    pos -- the coordinates at which to estimate the Christoffel symbols
    epsilon -- the step for estimating of the derivatives (default 1e-6)

    '''
    Gamma=numpy.empty((4, 4, 4))
    delta=jacobian_numerical(metric, pos, epsilon=epsilon)
    gup=metric.gmunu_up(pos)
    for i in range(4):
        for k in range(4):
            for l in range(4):
                Gamma[i, k, l] = (
                    0.5*gup[i, :]*
                    [
                        delta[l, :, k]
                        +delta[k, :, l]
                        -delta[:, k, l]
                    ]).sum()
    return Gamma


def check_christoffel(metric, poslist=None, epsilon=1e-6, abstol=1e-6, reltol=1e-6):
    '''Check the christoffel method of a gyoto.core.Metric subclass

    This method compares the Christoffel symbols of metric as given by
    metric.christoffel(pos) to those numerically estimated by
    christoffel_numerical(metric, pos). It raises an error if the
    difference is too large.

    The difference between the two estimates should always be smaller
    than abstol.

    In addition, if the value of the symbol is larger than abstol, the
    relative error should be smaller than reltol.

    Keyword arguments:
    metric -- one of:
              1- the gyoto.core.Metric instance to work on
                 (e.g. gyoto.std.KerrBL());
              2- a gyoto.core.Metric subclass
                 (e.g. gyoto.std.KerrBL);
              3- the name of a gyoto.core.metric subclass
                 (e.g. 'KerrBL').

    poslist -- a Python list of 4-coordinates at which to check the
              Christoffel symbols. By default, a small number of
              arbitrary positions that depend on whether the
              coordinate system is spherical or Cartesian and work
              well fr the Kerr metric are used.

    epsilon -- the step for estimating of the derivatives (default 1e-6)

    abstol -- the absolute tolerance

    retol -- the relative tolerance

    '''
    if isinstance(metric, str):
        metric=gyoto.core.Metric(metric)
    elif isinstance(metric, type):
        metric=metric()

    if poslist is None:
        if metric.coordKind()==gyoto.core.GYOTO_COORDKIND_SPHERICAL:
            poslist=[
                (0., 6., numpy.pi/2, 0.),
                (100., 50, numpy.pi/4, numpy.pi/6.)
                ]
        elif metric.coordKind()==gyoto.core.GYOTO_COORDKIND_CARTESIAN:
            poslist=[
                (0., 6., 0., 0.),
                (100., 50., 30., 50),
                (1000., 0., 0., 40.)
                ]
        else:
            raise ValueError('Unknown coordinate kind')

    for pos in poslist:
        G=metric.christoffel(pos)
        Gn=christoffel_numerical(metric, pos, epsilon)
        for a in range(4):
            for m in range(4):
                for n in range(4):
                    e=numpy.abs(G[a, m, n]-Gn[a, m, n])
                    assert e <= abstol, "absolute error {} larger than {} at {} for kind={}, alpha={}, mu={}, nu={}, val={}".format(e, abstol, pos, metric.kind(), a, m, n, G[a, m, n])
                    avg=numpy.abs(0.5*(G[a, m, n]-Gn[a, m, n]))
                    if avg > abstol:
                        assert e/avg <= reltol, "relative error {} larger than {} at {} for kind={}, alpha={}, mu={}, nu={}".format(e/avg, reltol, pos, metric.kind(), a, m, n)
