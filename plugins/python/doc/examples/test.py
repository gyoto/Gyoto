import gyoto

sp = gyoto.python.PythonSpectrum()
sp.Module     = "gyoto_sample_spectra"
sp.Class      = "PowerLaw"
sp.Parameters = (1., 10.)
print('value: {}'.format(sp(3e8/2e-6)))
sp.Parameters = (2., 20.)
print('value: {}'.format(sp(3e8/2e-6)))

#BlackBody6000 does not accept any parameter:
sp = gyoto.python.PythonSpectrum()
sp.Module = "gyoto_sample_spectra"
sp.Class  = "BlackBody6000"
print('value: {}'.format(sp(3e8/2e-6)))

sp.Module     = "gyoto_sample_spectra"
sp.Class      = "PowerLaw"
sp.Parameters = (2., 0.)
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2 = gyoto.spectrum.PowerLaw()
sp2.Exponent = 0.
sp2.Constant = 2.
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

sp.Module     = "gyoto_sample_spectra"
sp.Class      = "PowerLaw"
sp.Parameters = (2., 1.)
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2 = gyoto.spectrum.PowerLaw()
sp2.Exponent = 1.
sp2.Constant = 2.
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

sp.Module     = "gyoto_sample_spectra"
sp.Class      = "PowerLaw"
sp.Parameters = (5., -1.)
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2 = gyoto.spectrum.PowerLaw()
sp2.Exponent = -1.
sp2.Constant =  5.
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

gg=gyoto.python.PythonMetric()
gg.Module = "gyoto_sample_metrics"
gg.Class  = "Minkowski"
