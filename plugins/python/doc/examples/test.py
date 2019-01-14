import gyoto.core
import sys
import os.path

# The name of the Gyoto plug-in that can be loaded in a given Python
# session is the same as the name of the Python executable
python_plugin = os.path.basename(sys.executable)
gyoto.core.requirePlugin(python_plugin)

sp=gyoto.core.Spectrum("Python")
sp.set("Module", "gyoto_sample_spectra")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (1., 10.))
print('value: {}'.format(sp(3e8/2e-6)))
sp.set("Parameters", (2., 20.))
print('value: {}'.format(sp(3e8/2e-6)))

#BlackBody6000 does not accept any parameter:
sp=gyoto.core.Spectrum("Python")
sp.set("Module", "gyoto_sample_spectra")
sp.set("Class", "BlackBody6000")
print('value: {}'.format(sp(3e8/2e-6)))

sp.set("Module", "gyoto_sample_spectra")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (2., 0.))
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2=gyoto.core.Spectrum("PowerLaw")
sp2.set("Exponent", 0.)
sp2.set("Constant", 2.)
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

sp.set("Module", "gyoto_sample_spectra")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (2., 1.))
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2=gyoto.core.Spectrum("PowerLaw")
sp2.set("Exponent", 1.)
sp2.set("Constant", 2.)
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

sp.set("Module", "gyoto_sample_spectra")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (5., -1.))
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2=gyoto.core.Spectrum("PowerLaw")
sp2.set("Exponent", -1.)
sp2.set("Constant", 5.)
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

gg=gyoto.core.Metric("Python")
gg.set("Module", "gyoto_sample_metrics")
gg.set("Class", "Minkowski")
