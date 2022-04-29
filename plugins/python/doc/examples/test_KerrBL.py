# Compare the C++ and Python implementations of the KerrBL metric
#
# This script is for demonstration and test purposes, in all practical case
# the C++ implementation should be preferred.

# Import gyoto, require the Python plug-in. Note that the Python
# plug-in name depends on the running interpreter. Default to basename
# of Python executable, then python3, then python.
import gyoto.core, gyoto.std
import sys
import os.path

try:
    gyoto.core.requirePlugin(os.path.basename(sys.executable))
except gyoto.core.Error:
    try:
        gyoto.core.requirePlugin("python3")
    except gyoto.core.Error:
        try:
            gyoto.core.requirePlugin("python")
        except gyoto.core.Error:
            raise gyoto.core.Error("Could not load Python plugin, tried: "+os.path.basename(sys.executable)+", python3 and python")

# Instanciate an instance of the C++ KerrBL and an instance of the Python equivalent
kerrc=gyoto.std.KerrBL()
kerrp=gyoto.core.Metric("Python")
kerrp.set("Module", "gyoto_sample_metrics")
kerrp.set("Class", "KerrBL")

# Set spin
kerrc.spin(0.5)
kerrp.set("Parameters", (0.5,))

# Test a few methods
Rmsp=kerrp.getRms()
Rmsc=kerrc.getRms()
assert (Rmsp==Rmsc), "Rms is different"

Rmbp=kerrp.getRmb()
Rmbc=kerrc.getRmb()
assert (Rmbp==Rmbc), "Rmb is different"

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
assert (Wp==Wc), "Potential is different"

print ("All tests passed successfully")
