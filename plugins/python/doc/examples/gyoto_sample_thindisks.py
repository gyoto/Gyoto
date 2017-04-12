'''Sample Astrobj::ThinDisk for using with Gyoto Python plug-in

   Those classes demonstrate how to use Python classes as Gyoto
   Astrobj::ThibDisk implementations using Gyoto's "python"
   plug-in. Note that this plug-in can be renamed to whatever matches
   the particular version of Python it has been built against
   (e.g. python3.4).

   The goal is to be able to instantiate these from XML, from
   Yorick... and even from Python using the gyoto extension...

   Synopsis:

   import gyoto.core
   gyoto.core.requirePlugin("python") # or python2.7 or python3.4...
   td=gyoto.core.Astrobj("Python::ThinDisk")
   td.set("Module", "gyoto_sample_thindisks")
   td.set("Class", "ThinDisk")

   Classes that aim at implementing the Gyoto::Astrobj::ThinDisk
   interface do so by providing the following methods:

   getVelocity, giveDelta, emission, integrateEmission, transmission,
   __setitem__:
              optional.
   emission and integrateEmission can be overloaded by using the
   varargs argument.

'''

class ThinDisk:
    '''A ThinDisk with trivial emission
    '''
    def emission(self, nuem, dsem, cph, co):
        return dsem
