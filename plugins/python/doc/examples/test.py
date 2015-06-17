import gyoto
gyoto.loadPlugin("python3.4")
sp=gyoto.Spectrum("Python")
sp.set("Module", "gyoto_sample_callbacks")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (1., 10.))
print('value: {}'.format(sp(3e8/2e-6)))
sp.set("Parameters", (2., 20.))
print('value: {}'.format(sp(3e8/2e-6)))

#BlackBody6000 does not accept any parameter:
sp=gyoto.Spectrum("Python")
sp.set("Module", "gyoto_sample_callbacks")
sp.set("Class", "BlackBody6000")
print('value: {}'.format(sp(3e8/2e-6)))

sp.set("Module", "gyoto_sample_callbacks")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (2., 0.))
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2=gyoto.Spectrum("PowerLaw")
sp2.set("Exponent", 0.)
sp2.set("Constant", 2.)
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

sp.set("Module", "gyoto_sample_callbacks")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (2., 1.))
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2=gyoto.Spectrum("PowerLaw")
sp2.set("Exponent", 1.)
sp2.set("Constant", 2.)
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))

sp.set("Module", "gyoto_sample_callbacks")
sp.set("Class", "PowerLaw")
sp.set("Parameters", (5., -1.))
print('using Python integrate: {}'.format(sp.integrate(1., 2.)))

sp2=gyoto.Spectrum("PowerLaw")
sp2.set("Exponent", -1.)
sp2.set("Constant", 5.)
print('using Gyoto generic integrator: {}'.format(sp2.integrate(1., 2.)))
