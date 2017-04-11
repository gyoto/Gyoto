import numpy
import unittest
import gyoto
import gyoto_std

class TestStar(unittest.TestCase):

    def test_setInitCoord(self):
        st=gyoto_std.Star()
        gg=gyoto_std.KerrBL()
        st.metric(gg)
        pos_list=(600., 9., 1.5707999999999999741, 0)
        vel_list=(0., 0., 0.037037)
        st.setInitCoord(pos_list, vel_list)
        dst=numpy.zeros(8, float)
        dst2=numpy.zeros(8, float)
        st.getInitialCoord(dst)
        st.setPosition(pos_list)
        st.setVelocity(vel_list)
        st.getInitialCoord(dst2)
        self.assertTrue((dst == dst2).all())
