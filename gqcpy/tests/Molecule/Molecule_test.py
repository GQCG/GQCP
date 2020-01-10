import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class MoleculeTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testReadXYZ(self):
        """Check if the .ReadXYZ method works as expected"""

        self.molecule = gqcpy.Molecule.ReadXYZ("data/ch4_crawdad.xyz", 0)
        representation_string = "Number of electrons: 10\nC  (-0, 0, 0)\nH  (1.18377, -1.18377, -1.18377)\nH  (1.18377, 1.18377, 1.18377)\nH  (-1.18377, 1.18377, -1.18377)\nH  (-1.18377, -1.18377, 1.18377)\n"
        self.assertEqual(self.molecule.__repr__(), representation_string)


if __name__ == '__main__':
    unittest.main()
