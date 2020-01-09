import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class MullikenConstrainedQCM(unittest.TestCase):

    def setUp(self):
        """ Initiates variables to be used by tests  """
        N = gqcpy.Nucleus(7, 0, 0, 0)
        O = gqcpy.Nucleus(8, 7, 0, 0)
        NO = gqcpy.Molecule([N,O], +1) # The NO+ molecule with an intramolecular distance of 7 bohr
        basis_targets = [0,1,2,3,4] 
        basis_set = "STO-3G"
        
        self.constrained_module = gqcpy.MullikenConstrainedFCI(NO, basis_set, basis_targets)

        # Reference data obtained from gqcp NO+ calculations
        self.lambda_input = -2.27812975203

        self.reference_energy = -124.450111059414
        self.reference_population = 4.37000013448507
        self.reference_entropy = 0.294476160440132
        self.reference_N_fragment_energy = -51.4005213413705
        self.reference_self_N_fragment_energy = -51.0749111177636
        self.reference_O_fragment_energy = -73.0495897180439
        self.reference_self_O_fragment_energy = -72.7239794944369
        self.reference_interaction_energy = -0.651220447214087

        # Initiate a second set to test spin constraints
        O_2b = gqcpy.Nucleus(8, 2, 0, 0)
        NO_2b = gqcpy.Molecule([N,O_2b], +1)  # the NO+ molecule with an intramolecular distance of 2 bohr

        self.constrained_module_spin = gqcpy.MullikenConstrainedFCI(NO_2b, basis_set, basis_targets)
        self.atomic_sz = 0.938336864443

    def tearDown(self):
        pass

    def test_properties(self):
        """ compare properties with reference """
        self.constrained_module.solveMullikenDavidson(self.lambda_input)
        self.assertAlmostEqual(self.constrained_module.get_energy(), self.reference_energy)
        self.assertAlmostEqual(self.constrained_module.get_population(), self.reference_population)
        #self.assertAlmostEqual(self.constrained_module.get_entropy(), self.reference_entropy) entropy is sensitive to degeneracies is disabled
        self.assertAlmostEqual(self.constrained_module.get_A_fragment_energy(), self.reference_N_fragment_energy)
        self.assertAlmostEqual(self.constrained_module.get_A_fragment_self_energy(), self.reference_self_N_fragment_energy)
        self.assertAlmostEqual(self.constrained_module.get_B_fragment_energy(), self.reference_O_fragment_energy)
        self.assertAlmostEqual(self.constrained_module.get_B_fragment_self_energy(), self.reference_self_O_fragment_energy)
        self.assertAlmostEqual(self.constrained_module.get_interaction_energy(), self.reference_interaction_energy)

    def test_spinZ(self):
        """ compare spin-z with reference """
        self.constrained_module_spin.solveMullikenDavidson(0, 0.5)
        self.assertAlmostEqual(self.constrained_module_spin.get_sz(), self.atomic_sz)

    
if __name__ == '__main__':
    unittest.main()
