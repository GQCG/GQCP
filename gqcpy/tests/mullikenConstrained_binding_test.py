import unittest
import numpy as np

# force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

# import our `pybind11`-based extension module
import gqcpy


class MullikenConstrainedQCM(unittest.TestCase):

    ''' iniates variables to be used by tests '''
    def setUp(self):
        N = gqcpy.Nuclei(7, 0, 0, 0)
        O = gqcpy.Nuclei(8, 7, 0, 0)
        NO = gqcpy.Molecule([N,O], +1)
        self.constrained_module = gqcpy.MullikenConstrainedFCI()
        "-124.450111059414	-2.27812975203	4.37000013448507	0.294476160440132	-51.4005213413705	-51.0749111177636	-73.0495897180439	-72.7239794944369	-0.651220447214087"


    def tearDown(self):
        pass

        ''' compare energies with reference '''
    def test_energies(self):
        self.assertAlmostEqual(self.energy1, self.ref_energy1)
        self.assertAlmostEqual(self.energy2, self.ref_energy2)

        ''' compare RDMs with reference '''
    def test_1rdms(self):
        self.assertTrue(np.allclose(self.rdm1, self.ref_rdm1))
        self.assertTrue(np.allclose(self.rdm2, self.ref_rdm2))

    
if __name__ == '__main__':
    unittest.main()
