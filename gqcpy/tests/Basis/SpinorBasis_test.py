# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy

import unittest
import numpy as np


class SpinorBasisTest(unittest.TestCase):

    def setUp(self):
        """Initiate variables to be used by the tests"""

        self.molecule = gqcpy.Molecule.ReadXYZ("data/h2_szabo.xyz" , 0)  # create a neutral molecule
        self.spinor_basis = gqcpy.SpinorBasis(self.molecule, "STO-3G")

        self.ref_S = np.array([[1.0,    0.6593],
                               [0.6593, 1.0]])

        self.ref_H_core = np.array([[-1.1204, -0.9584],
                                    [-0.9584, -1.1204]])


    def tearDown(self):
        pass


    def test_one_electron_integrals(self):
        """Compare the one-electron integrals with the reference values"""

        S = self.spinor_basis.quantizeOverlapOperator().parameters()
        
        T = self.spinor_basis.quantizeKineticOperator().parameters()
        V = self.spinor_basis.quantizeNuclearAttractionOperator(self.molecule).parameters()
        H_core = T + V
        
        self.assertTrue(np.allclose(S, self.ref_S, 1.0e-04))
        self.assertTrue(np.allclose(H_core, self.ref_H_core, 1.0e-04))


    def test_two_electron_integrals(self):
        """Compare the two-electron integrals with the reference values"""

        g = self.spinor_basis.quantizeCoulombRepulsionOperator().parameters()

        self.assertAlmostEqual(g[(0,0,0,0)], 0.7746, 4)
        self.assertAlmostEqual(g[(0,0,0,0)], g[(1,1,1,1)], 12)
        self.assertAlmostEqual(g[(0,0,1,1)], 0.5697, 4)
        self.assertAlmostEqual(g[(1,0,0,0)], 0.4441, 4)
        self.assertAlmostEqual(g[(1,0,0,0)], g[(1,1,1,0)], 12)
        self.assertAlmostEqual(g[(1,0,1,0)], 0.2970, 4)


if __name__ == '__main__':
    unittest.main()
