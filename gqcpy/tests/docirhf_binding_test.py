import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class DOCIRHFQCM(unittest.TestCase):

    def setUp(self):
        """ Iniates variables to be used by tests """
        self.doci_module = gqcpy.DOCIRHF("data/CO_mulliken.xyz", "STO-3G", False)
        self.ref_energy = -111.29189323578304
        self.doci_module.solve()

    def tearDown(self):
        pass

    def test_energies(self):
        """ Compare the energy with a reference value """
        self.assertAlmostEqual(self.doci_module.get_energy(), self.ref_energy)

    
if __name__ == '__main__':
    unittest.main()
