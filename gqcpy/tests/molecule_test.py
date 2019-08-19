import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class Molecule(unittest.TestCase):

    def setUp(self):
         """ iniates variables to be used by tests """
        pass

    def tearDown(self):
        pass

    def test_constructor(self):
        """ compare properties with reference """
        self.molecule = gqcpy.Molecule.ReadXYZ("data/ch4_crawdad.xyz", 0)
        representation_string = "Number of electrons: 10\nC  (-0, 0, 0)\nH  (1.18377, -1.18377, -1.18377)\nH  (1.18377, 1.18377, 1.18377)\nH  (-1.18377, 1.18377, -1.18377)\nH  (-1.18377, -1.18377, 1.18377)\n"
        self.assertEqual(self.molecule.__repr__(), representation_string)
    
if __name__ == '__main__':
    unittest.main()
