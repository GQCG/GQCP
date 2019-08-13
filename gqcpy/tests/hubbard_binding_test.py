import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class HubbardQCM(unittest.TestCase):

    def setUp(self):
        """ Iniates variables to be used by tests """
        self.csline_input = "-0.999984,-0.736924,0.511211,-0.082700,0.065534,-0.562082,-0.905911,0.357729,0.358593,0.869386"
        self.ref_energy1 = -3.49379514792384
        self.ref_energy2 = -3.01890254187003

        self.ref_rdm1 = np.array([[ 1.490606279396268,  0.3170311736848896, -0.7605319143378515,   0.2597814689331624],
                                   [ 0.3170311736848896,  1.181067587214437,  0.5392166542874786,   0.7488699568813794],
                                   [-0.7605319143378515,  0.5392166542874786,  0.8414178179099283,  0.2968598589512055],
                                   [ 0.2597814689331624,  0.7488699568813794,  0.2968598589512055,  0.4869083154793552]])

        self.ref_rdm2 = np.array([[ 1.542976978440379, 0.6158331134876798, -0.2525249643990615, 0.05861966911629712],
                                   [ 0.6158331134876798, 0.8454119950562293,   0.2400544876984531, 0.467304942268566],
                                   [-0.2525249643990615, 0.2400544876984531, 0.7928298261486432, -0.2741228345206897],
                                   [ 0.05861966911629712, 0.467304942268566, -0.2741228345206897, 0.8187812003547599]])
        
        self.hubbardQCM = gqcpy.Hubbard(self.csline_input, 2, 2, 2)
        self.hubbardQCM.solve()
        energies = self.hubbardQCM.get_energies()
        self.energy1 = energies[0]
        self.energy2 = energies[1]
        rdms = self.hubbardQCM.get_one_rdms()
        self.rdm1 = rdms[0]
        self.rdm2 = rdms[1]

    def tearDown(self):
        pass

    def test_energies(self):
        """ Compare energes with reference values """
        self.assertAlmostEqual(self.energy1, self.ref_energy1)
        self.assertAlmostEqual(self.energy2, self.ref_energy2)

    def test_1rdms(self):
        """ Compare RDMs with reference """
        self.assertTrue(np.allclose(self.rdm1, self.ref_rdm1))
        self.assertTrue(np.allclose(self.rdm2, self.ref_rdm2))

    
if __name__ == '__main__':
    unittest.main()
