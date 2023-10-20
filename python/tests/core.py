#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import gyoto.core
import numpy
import gyoto.metric, gyoto.astrobj, gyoto.spectrum, gyoto.spectrometer
import inspect

gyoto.core.requirePlugin('stdplug')

class TestSmartPointer(unittest.TestCase):
    def test_simple_classes(self):
        for classname in ('Scenery', 'Screen', 'Photon'):
            cls=getattr(gyoto.core, classname)
            obj=cls()
            self.assertEqual(obj.getRefCount(), 1)
            clone=obj.clone()
            self.assertEqual(obj.getRefCount(), 1)
            self.assertEqual(clone.getRefCount(), 1)
            rep=obj.__str__()
            self.assertEqual(obj.getRefCount(), 1)

    def test_base_classes(self):
        '''Test reference counting

        All constructors, cloners and destructors must implement and
        decrement the reference counter correctly.

        '''
        for gnspace in ('Metric', 'Astrobj', 'Spectrum', 'Spectrometer'):
            pnspace=gnspace.lower()
            nspace=getattr(gyoto, pnspace)
            generic=getattr(nspace, 'Generic')
            for classname, cls in inspect.getmembers(nspace):
                # Skip abstract classes
                if (classname in ('Generic',
                                  gnspace,
                                  'StandardAstrobj',
                                  'UniformSphere')
                    or not inspect.isclass(cls)):
                    continue
                # The XML name of ComplexAstrobj et al. is 'Complex'
                if 'Complex' in classname:
                    classname='Complex'
                # The XML name of UniformSpectrometer is 'wave'
                if classname in ('UniformSpectrometer', 'Uniform'):
                    classname='wave'
                # Construct instance from default constructor
                obj=cls()
                self.assertEqual(obj.getRefCount(), 1)
                # Cast to base class
                gen=generic(obj)
                self.assertEqual(obj.getRefCount(), 2)
                # Destroy one reference
                del gen
                self.assertEqual(obj.getRefCount(), 1)
                # Clone
                clone=obj.clone()
                self.assertEqual(obj.getRefCount(), 1)
                self.assertEqual(clone.getRefCount(), 1)
                # Print
                rep=obj.__str__()
                self.assertEqual(obj.getRefCount(), 1)
                # Clean
                del rep
                del clone
                del obj
                # Construct instance from XML name
                gen=generic(classname)
                self.assertEqual(gen.getRefCount(), 1)
                # Cast to derived class
                obj=cls(gen)
                self.assertEqual(gen.getRefCount(), 2)
                # Destroy one isntance
                del obj
                self.assertEqual(gen.getRefCount(), 1)
                # Clean
                del gen
                # Construct instance from XML name, setting plugin list
                gen=generic(classname, [])
                self.assertEqual(gen.getRefCount(), 1)
                # Clean
                del gen

    def test_complex_classes(self):
        '''Test that adding, retrieving, deleting members updates refCount
        '''
        for cls in (gyoto.astrobj.Complex, gyoto.spectrometer.Complex):
            cplx=cls()
            cplx1=cls()
            self.assertEqual(cplx.getRefCount(), 1)
            self.assertEqual(cplx1.getRefCount(), 1)
            cplx.append(cplx1)
            self.assertEqual(cplx1.getRefCount(), 2)
            cplx2=cplx[0]
            self.assertEqual(cplx1.getRefCount(), 3)
            del cplx2
            self.assertEqual(cplx1.getRefCount(), 2)
            cplx.remove(0)
            self.assertEqual(cplx1.getRefCount(), 1)
            del cplx1
            del cplx

class TestUnit(unittest.TestCase):

    def test___str__(self):
        self.assertEqual(str(gyoto.core.Unit("km")), "km")

    def test_throw(self):
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Unit("deadbeaf"))

    def test_to(self):
        self .assertEqual(gyoto.core.Unit("m").To(5, gyoto.core.Unit("km")), 5000.)

    def test_from(self):
        self .assertEqual(gyoto.core.Unit("m").From(5000, gyoto.core.Unit("km")), 5.)

class TestConverter(unittest.TestCase):

    def test_trivial(self):
        a=gyoto.core.Converter()
        self.assertEqual(a(2.), 2.)

    def test_converter(self):
        a=gyoto.core.Converter(gyoto.core.Unit("m"), gyoto.core.Unit("km"))
        self.assertEqual(5., a(5000))

    def test_reset_trivial(self):
        a=gyoto.core.Converter(gyoto.core.Unit("m"), gyoto.core.Unit("km"))
        a.reset()
        self.assertEqual(a(2.), 2.)

    def test_reset(self):
        a=gyoto.core.Converter()
        a.reset(gyoto.core.Unit("m"), gyoto.core.Unit("km"))
        self.assertEqual(5., a(5000))

class TestError(unittest.TestCase):

    def test_get_message(self):
        a=gyoto.core.Error("toto")
        self.assertEqual(a.get_message(), "toto")

    def test_throwError(self):
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.throwError("msg"))

    def test_getErrcode(self):
        a=gyoto.core.Error("toto")
        self.assertEqual(a.getErrcode(), 1)

class TestMetric(unittest.TestCase):

    def test_circularVelocity(self):
        gg=gyoto.core.Metric('KerrBL')
        vel=numpy.zeros(4, float)
        vel2=numpy.zeros(4, float)
        pos_list=[0, 6, numpy.pi*0.5, 0]
        pos_numpy=numpy.asarray(pos_list)
        gg.circularVelocity(pos_numpy, vel)
        gg.circularVelocity(pos_list, vel2)
        vel3=gg.circularVelocity(pos_list)
        self.assertTrue((vel == vel2).all())
        self.assertTrue((vel == vel3).all())

    def test_SysPrimeToTdot(self):
        gg=gyoto.core.Metric('KerrBL')
        self.assertAlmostEqual(gg.SysPrimeToTdot([0, 6, numpy.pi/2, 0], (0, 0, 0.1)), 1.8057877962865378)

    def test_isStopCondition(self):
        gg=gyoto.core.Metric('KerrBL')
        self.assertFalse(gg.isStopCondition((0, 6, 3.14, 0, 0, 0, 0, 0)))
        self.assertTrue(gg.isStopCondition((0, 0, 3.14, 0, 0, 0, 0, 0)))

    def test_gmunu(self):
        gg=gyoto.core.Metric('KerrBL')
        tt=gg.gmunu((0, 6, 3.14, 0), 0, 0)
        self.assertAlmostEqual(tt, -0.6666666666666667)
        dst=gg.gmunu((0, 6, 3.14, 0))
        self.assertEqual(tt, dst[0, 0])

    def test_christoffel(self):
        gg=gyoto.core.Metric('KerrBL')
        tt=gg.christoffel((0, 6, 3.14, 0), 0, 0, 0)
        self.assertAlmostEqual(tt, 0)
        dst=numpy.zeros((4, 4, 4), float)
        retval=gg.christoffel(dst, (0, 6, 3.14, 0))
        self.assertEqual(tt, dst[0, 0, 0])
        self.assertEqual(retval, 0)
        dst2=gg.christoffel((0, 6, 3.14, 0))
        self.assertEqual(tt, dst2[0, 0, 0])

    def test_norm(self):
        gg=gyoto.core.Metric('KerrBL')
        gg.set('Spin', 0.95)
        pos=numpy.asarray([0., 6., numpy.pi/2., 0.])
        vel=gg.circularVelocity(pos)
        self.assertAlmostEqual(gg.norm(pos, vel), -1.)

    def test_GramSchmidt_BL(self):
        gg=gyoto.core.Metric('KerrBL')
        gg.set('Spin', 0.95)
        pos=numpy.asarray([0., 6., numpy.pi/2., 0.])
        fourvel=gg.circularVelocity(pos)
        screen1=numpy.zeros(4)
        screen2=numpy.zeros(4)
        screen3=numpy.zeros(4)
        screen1[0]=0.;
        screen1[1]=0.;
        screen1[2]=0.;
        screen1[3]=-1.;
        screen2[0]=0.;
        screen2[1]=0.;
        screen2[2]=-1.;
        screen2[3]=0.;
        screen3[0]=0.;
        screen3[1]=-1.;
        screen3[2]=0.;
        screen3[3]=0.;
        gg.GramSchmidt(pos, fourvel, screen2, screen3, screen1);
        self.assertAlmostEqual(gg.norm(pos, fourvel), -1.)
        self.assertAlmostEqual(gg.norm(pos, screen1), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen2), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen3), 1.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen1), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen2, screen3), 0.)

    def test_GramSchmidt_KS(self):
        gg=gyoto.core.Metric('KerrKS')
        gg.set('Spin', 0.95)
        pos=numpy.asarray([0., 6., 0. , 0.])
        fourvel=gg.circularVelocity(pos)
        screen1=numpy.zeros(4)
        screen2=numpy.zeros(4)
        screen3=numpy.zeros(4)
        rp=numpy.sqrt(pos[1]*pos[1]+pos[2]*pos[2])
        theta=numpy.arctan2(rp, pos[3])
        phi=numpy.arctan2(pos[2], pos[1])
        sp=numpy.sin(phi)
        cp=numpy.cos(phi)
        st=numpy.sin(theta)
        ct=numpy.cos(theta)
        screen1[0]=0.;
        screen1[1]=sp;
        screen1[2]=-cp;
        screen1[3]=0.;
        screen2[0]=0.;
        screen2[1]=-ct*cp;
        screen2[2]=-ct*sp;
        screen2[3]=st;
        screen3[0]=0.;
        screen3[1]=-pos[1];
        screen3[2]=-pos[2];
        screen3[3]=-pos[3];
        gg.GramSchmidt(pos, fourvel, screen2, screen3, screen1);
        self.assertAlmostEqual(gg.norm(pos, fourvel), -1.)
        self.assertAlmostEqual(gg.norm(pos, screen1), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen2), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen3), 1.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen1), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen2, screen3), 0.)

    def test_obsereverTetrad_BL(self):
        gg=gyoto.core.Metric('KerrBL')
        gg.set('Spin', 0.95)
        pos=numpy.asarray([0., 6., numpy.pi/2., 0.])
        fourvel=gg.circularVelocity(pos)
        screen1=numpy.zeros(4)
        screen2=numpy.zeros(4)
        screen3=numpy.zeros(4)
        gg.observerTetrad(pos, fourvel, screen1, screen2, screen3)
        self.assertAlmostEqual(gg.norm(pos, fourvel), -1.)
        self.assertAlmostEqual(gg.norm(pos, screen1), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen2), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen3), 1.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen1), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen2, screen3), 0.)

    def test_GramSchmidt_KS(self):
        gg=gyoto.core.Metric('KerrKS')
        gg.set('Spin', 0.95)
        pos=numpy.asarray([0., 6., 0. , 0.])
        fourvel=gg.circularVelocity(pos)
        screen1=numpy.zeros(4)
        screen2=numpy.zeros(4)
        screen3=numpy.zeros(4)
        gg.observerTetrad(pos, fourvel, screen1, screen2, screen3)
        self.assertAlmostEqual(gg.norm(pos, fourvel), -1.)
        self.assertAlmostEqual(gg.norm(pos, screen1), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen2), 1.)
        self.assertAlmostEqual(gg.norm(pos, screen3), 1.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen1), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, fourvel, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen2), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen1, screen3), 0.)
        self.assertAlmostEqual(gg.ScalarProd(pos, screen2, screen3), 0.)

class TestValue(unittest.TestCase):

    def test_good(self):
        a=gyoto.core.Value(5)
        self.assertEqual(a.toLong(), 5)
        self.assertEqual(a.toULong(), 5)
        self.assertEqual(a.toSizeT(), 5)
        a=gyoto.core.Value((1, 2, 3.))
        self.assertEqual(a.toVDouble(), (1., 2., 3.))
        a=gyoto.core.Value((1, 2, 3))
        self.assertEqual(a.toVULong(), (1, 2, 3))
        a=gyoto.core.Value('a')
        self.assertEqual(a.toString(), 'a')
        s=gyoto.core.Screen()
        a=gyoto.core.Value(s)
        self.assertEqual(a.toScreen().this, s.this)

    def test_assign(self):
        a=gyoto.core.Value(0)
        a.assign(gyoto.core.Value(5))
        self.assertEqual(a.toLong(), 5)
        a.assign(gyoto.core.Value(5.))
        self.assertEqual(a.toDouble(), 5.)

    def test_bad(self):
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value(5).toDouble())
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value(5.).toLong())
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value((1,)).toVDouble())
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value('a').toVULong())

class TestProperty(unittest.TestCase):

    def test_get(self):
        s=gyoto.core.Screen()
        self.assertEqual(s.get('Distance'), 1.0)
        self.assertRaises(gyoto.core.Error, lambda: s.get('NonExistentProperty'))

    def test_set(self):
        s=gyoto.core.Screen()
        s.set('Distance', 8., 'kpc')
        self.assertAlmostEqual(s.get('Distance'), 8.*gyoto.core.GYOTO_KPC, -15)

    def test_describe(self):
        s=gyoto.core.Screen()
        p=s.property('Distance')
        self.assertIn('Distance: double with unit', s.describeProperty(p))

class TestPolar(unittest.TestCase):
    def test_triad(self):
        met=gyoto.core.Metric("KerrBL")
        met.set("Spin",0.99)
        met.set("Mass",4e6, "sunmass")
        s=gyoto.core.Screen()
        s.metric(met)
        s.distance(8., 'kpc')
        s.inclination(95,"°")
        s.PALN(numpy.pi)
        s.argument(-numpy.pi/2)
        s.time(8., 'kpc')
        coord=numpy.zeros(8, float)
        Ephi=numpy.zeros(4, float)
        Etheta=numpy.zeros(4, float)

        fov=4.85e-10 # fov in rad for typical Sgr conditions
        #s.getRayCoord(0.,0., coord)
        compute_polar_basis=True
        s.getRayTriad(fov/2.,fov/2.,
                      coord,
                      compute_polar_basis,
                      Ephi, Etheta)
        #s.getRayTriad(coord, Ephi, Etheta)
        x = coord[0:4]
        k = coord[4:8]
        #print("")
        #print("Photon pos= ",x)
        #print("Photon tangent= ",k)
        #print("Ephi= ",Ephi)
        #print("Etheta= ",Etheta)
        self.assertAlmostEqual(met.ScalarProd(x, k, Ephi), 0., 15)
        self.assertAlmostEqual(met.ScalarProd(x, k, Etheta), 0., 15)
        self.assertAlmostEqual(met.ScalarProd(x, Etheta, Ephi), 0., 15)
        self.assertAlmostEqual(met.ScalarProd(x, Etheta, Etheta), 1., 15)
        self.assertAlmostEqual(met.ScalarProd(x, Ephi, Ephi), 1., 15)
        ao = gyoto.core.Astrobj("Complex")
        ao.rMax(0.)
        ph = gyoto.core.Photon()
        ph.parallelTransport(True)
        ph.setInitialCondition(met, ao, coord, Ephi, Etheta);
        cph=gyoto.core.vector_double()
        ph.getCoord(10., cph)
        coord=numpy.asarray(cph)
        x = coord[0:4]
        k = coord[4:8]
        Ephi = coord[8:12]
        Etheta = coord[12:16]
        self.assertAlmostEqual(met.ScalarProd(x, k, Ephi), 0., 6)
        self.assertAlmostEqual(met.ScalarProd(x, k, Etheta), 0., 6)
        self.assertAlmostEqual(met.ScalarProd(x, Etheta, Ephi), 0., 6)
        self.assertAlmostEqual(met.ScalarProd(x, Etheta, Etheta), 1., 6)
        self.assertAlmostEqual(met.ScalarProd(x, Ephi, Ephi), 1., 6)
