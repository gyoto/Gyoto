import unittest
import gyoto

class TestError(unittest.TestCase):

    def test_get_message(self):
        a=gyoto.Error("toto")
        self.assertEqual(a.get_message(), "toto")

    def test_throwError(self):
        self.assertRaises(gyoto.Error, lambda: gyoto.throwError("msg"))
