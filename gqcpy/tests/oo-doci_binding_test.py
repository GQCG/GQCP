import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class HubbardQCM(unittest.TestCase):

    def setUp(self):
        """ Iniates variables to be used by tests """
        self.oo_doci_module = gqcpy.DOCINewtonOrbitalOptimizer("data/h2_cristina.xyz", "STO-3G", False, False)
        self.ref_energy = -1.13726333769813
        self.oo_doci_module.solve()

    def tearDown(self):
        pass

    def test_energies(self):
        """ Compare the energy with a reference value """
        self.assertAlmostEqual(self.oo_doci_module.get_energy(), self.ref_energy)

    
if __name__ == '__main__':
    unittest.main()
