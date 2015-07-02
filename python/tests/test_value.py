import unittest
import gyoto
gyoto.loadPlugin("stdplug")

class TestValue(unittest.TestCase):

    def test_good(self):
        a=gyoto.Value(5)
        self.assertEqual(a.toLong(), 5)
        self.assertEqual(a.toULong(), 5)
        self.assertEqual(a.toSizeT(), 5)
        a=gyoto.Value((1, 2, 3.))
        self.assertEqual(a.toVDouble(), (1., 2., 3.))
        a=gyoto.Value((1, 2, 3))
        self.assertEqual(a.toVULong(), (1, 2, 3))
        a=gyoto.Value('a')
        self.assertEqual(a.toString(), 'a')
        s=gyoto.Screen()
        a=gyoto.Value(s)
        self.assertEqual(a.toScreen().this, s.this)

    def test_assign(self):
        a=gyoto.Value(0)
        a.assign(gyoto.Value(5))
        self.assertEqual(a.toLong(), 5)
        a.assign(gyoto.Value(5.))
        self.assertEqual(a.toDouble(), 5.)

    def test_bad(self):
        self.assertRaises(gyoto.Error, lambda: gyoto.Value(5).toDouble())
        self.assertRaises(gyoto.Error, lambda: gyoto.Value(5.).toLong())
        self.assertRaises(gyoto.Error, lambda: gyoto.Value((1,)).toVDouble())
        self.assertRaises(gyoto.Error, lambda: gyoto.Value('a').toVULong())
