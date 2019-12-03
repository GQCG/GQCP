import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class DOCIRHFQCM(unittest.TestCase):

    def setUp(self):
        """ Iniates variables to be used by tests """
        F = gqcpy.Nucleus(9, 0, 0, 0)
        self.doci_module = gqcpy.DOCIRHF([F, +1], "STO-3G", False)
        self.ref_energy = -98.32891803226985
        self.doci_module.solve()

    def tearDown(self):
        pass

    def test_energies(self):
        """ Compare the energy with a reference value """
        print(doci_module.get_energy())
        self.assertAlmostEqual(self.doci_module.get_energy(), self.ref_energy)

    
if __name__ == '__main__':
    unittest.main()
