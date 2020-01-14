import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class DOCIRHFTest(unittest.TestCase):

    def setUp(self):
        """Initiate variables to be used by the tests"""

        # Construct the F+ ion
        F = gqcpy.Nucleus(9, 0, 0, 0)
        self.doci_module = gqcpy.DOCIRHF(gqcpy.Molecule([F], +1), "STO-3G", False)

        self.ref_energy = -97.41595197133647
        self.doci_module.solve()


    def tearDown(self):
        pass


    def testEnergies(self):
        """Compare the energy with a reference value"""
        self.assertAlmostEqual(self.doci_module.get_energy(), self.ref_energy)

    
if __name__ == '__main__':
    unittest.main()
