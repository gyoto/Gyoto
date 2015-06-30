import unittest
import gyoto

class TestError(unittest.TestCase):

    def test_get_message(self):
        a=gyoto.Error("toto")
        self.assertEqual(a.get_message(), "toto")

if __name__ == '__main__':
    unittest.main()
