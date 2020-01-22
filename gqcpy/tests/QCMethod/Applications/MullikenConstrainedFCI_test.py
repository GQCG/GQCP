import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class MullikenConstrainedFCITest(unittest.TestCase):

    def setUp(self):
        """Initiate variables to be used by tests"""
        # Create NO+ with an intramolecular distance of 7 bohr
        N = gqcpy.Nucleus(7, 0, 0, 0)
        O = gqcpy.Nucleus(8, 7, 0, 0)
        NO_plus_far = gqcpy.Molecule([N,O], +1)
        basis_targets = [0,1,2,3,4]  # the basis functions that belong to N
        basis_set = "STO-3G"

        self.constrained_module_mulliken = gqcpy.MullikenConstrainedFCI(NO_plus_far, basis_set, basis_targets)

        # Reference data obtained from previous GQCP NO+ calculations
        self.lambda_input = -2.27812975203

        self.reference_energy = -124.450111059414
        self.reference_population = 4.37000013448507
        self.reference_N_fragment_energy = -51.4005213413705
        self.reference_self_N_fragment_energy = -51.0749111177636
        self.reference_O_fragment_energy = -73.0495897180439
        self.reference_self_O_fragment_energy = -72.7239794944369
        self.reference_interaction_energy = -0.651220447214087


        # Create NO+ with an intramolecular distance of 2 bohr
        O = gqcpy.Nucleus(8, 2, 0, 0)
        NO_plus_close = gqcpy.Molecule([N,O], +1)

        self.constrained_module_spin = gqcpy.MullikenConstrainedFCI(NO_plus_close, basis_set, basis_targets)
        self.atomic_sz = 0.938336864443


    def tearDown(self):
        pass


    def test_properties(self):
        """Compare some properties with values previously calculated with GQCP. Since the first test system is NO+ at 7 bohr apart, we aren't too strict on the equality."""

        self.constrained_module_mulliken.solveMullikenDavidson(mulliken_multiplier = self.lambda_input)  # lambda is the Mulliken-constraint multiplier

        # Compare the solutions to reference values
        self.assertAlmostEqual(self.constrained_module_mulliken.get_energy(), self.reference_energy, places=4)
        self.assertAlmostEqual(self.constrained_module_mulliken.get_population(), self.reference_population, places=4)
        self.assertAlmostEqual(self.constrained_module_mulliken.get_A_fragment_self_energy(), self.reference_self_N_fragment_energy, places=4)
        self.assertAlmostEqual(self.constrained_module_mulliken.get_B_fragment_energy(), self.reference_O_fragment_energy, places=4)
        self.assertAlmostEqual(self.constrained_module_mulliken.get_B_fragment_self_energy(), self.reference_self_O_fragment_energy, places=4)
        self.assertAlmostEqual(self.constrained_module_mulliken.get_interaction_energy(), self.reference_interaction_energy, places=4)


    def test_spinZ(self):
        """Compare the expectation value of the atomic S_z operator with values previously calculated with GQCP."""

        self.constrained_module_spin.solveMullikenDavidson(mulliken_multiplier = 0, sz_multiplier = 0.5)
        self.assertAlmostEqual(self.constrained_module_spin.get_sz(), self.atomic_sz, places=4)


if __name__ == '__main__':
    unittest.main()
