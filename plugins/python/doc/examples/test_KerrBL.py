# Compare the C++ and Python implementations of the KerrBL metric
#
# This script is for demonstration and test purposes, in all practical case
# the C++ implementation should be preferred.

# Import gyoto, require the Python plug-in. Note that the Python
# plug-in name depends on the running interpreter. Default to basename
# of Python executable, then python3, then python.
import gyoto, gyoto.python, gyoto_sample_metrics, numpy

# Instanciate an instance of the C++ KerrBL and an instance of the Python equivalent
kerrc=gyoto.metric.KerrBL()
kerrp=gyoto.metric.PythonMetric()
kerrp.Module = "gyoto_sample_metrics"
kerrp.Class  = "KerrBL"

# Set spin
kerrc.Spin       =  0.5
kerrp.Parameters = (0.5,)

# Test getRms
Rmsp=kerrp.getRms()
Rmsc=kerrc.getRms()
assert (Rmsp==Rmsc), "Rms is different"

# Test getRmb
Rmbp=kerrp.getRmb()
Rmbc=kerrc.getRmb()
assert (Rmbp==Rmbc), "Rmb is different"

# Test getSpecificAngularMomentum
samp=kerrp.getSpecificAngularMomentum(10.)
samc=kerrc.getSpecificAngularMomentum(10.)
assert (samp==samc), "Specific angular momentum is different"

# Test gmunu
pos=[0, 10, 1, 1.]
ggp=kerrp.gmunu(pos)
ggc=kerrc.gmunu(pos)
assert (ggp==ggc).all(), "gmunu is different"

# Test christoffel
pos=[0, 10, 1, 1.]
chp=kerrp.christoffel(pos)
chc=kerrp.christoffel(pos)
assert (chp==chc).all(), "christoffel is different"

# Test getPotential
Wp=kerrp.getPotential([0, 10, 1, 1.], 10.)
Wc=kerrc.getPotential([0, 10, 1, 1.], 10.)
assert (Wp==Wc), "getPotential is different"


# Test isStopCondition
Wp=kerrp.isStopCondition([0, 10, 1, 1., 0., 0., 0., 0])
Wc=kerrc.isStopCondition([0, 10, 1, 1., 0., 0., 0., 0])
assert (Wp==Wc), "isStopCondition is different"

# Test circularVelocity
kerrp.keplerian(False)
kerrc.keplerian(False)
vp=kerrp.circularVelocity([0, 10, 1, 1.])
vc=kerrc.circularVelocity([0, 10, 1, 1.])
assert (vp==vc).all(), "isStopCondition is different"

vp=numpy.ndarray(4, dtype=float)
vc=numpy.ndarray(4, dtype=float)
kerrp.circularVelocity([0, 10, 1, 1.], vp)
kerrc.circularVelocity([0, 10, 1, 1.], vc)
assert (vp==vc).all(), "isStopCondition is different"

vp=numpy.ndarray(4, dtype=float)
vc=numpy.ndarray(4, dtype=float)
kerrp.circularVelocity([0, 10, 1, 1.], vp, -1)
kerrc.circularVelocity([0, 10, 1, 1.], vc, -1)
assert (vp==vc).all(), "isStopCondition is different"

print ("All tests passed successfully")
