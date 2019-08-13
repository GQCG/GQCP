import unittest
import numpy as np

# force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

# import our `pybind11`-based extension module
import gqcpy


class MullikenConstrainedQCM(unittest.TestCase):

    def setUp(self):
        """ Iniates variables to be used by tests  """
        N = gqcpy.Nuclei(7, 0, 0, 0)
        O = gqcpy.Nuclei(8, 7, 0, 0)
        NO = gqcpy.Molecule([N,O], +1) # The NO+ molecule with an intramolecular distance of 7 bohr
        self.constrained_module = gqcpy.MullikenConstrainedFCI()

        # Reference data obtained from gqcp NO+ calculations
        self.reference_energy = -124.450111059414
        self.reference_population = 4.37000013448507
        self.reference_entropy = 0.294476160440132
        self.reference_N-fragment_energy = -51.4005213413705
        self.reference_self_N-fragment_energy = -51.0749111177636
        self.reference_self_N-fragment_energy = -51.0749111177636
        self.reference_O-fragment_energy = -73.0495897180439
        self.reference_self_O-fragment_energy = -72.7239794944369
        self.reference_interaction_energy = -0.651220447214087

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
