import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class FukuiDysonAnalysisTest(unittest.TestCase):

    def setUp(self):
        """Initiate variables to be used by the tests"""

        # Create a water molecule
        O = gqcpy.Nucleus(8, 0.0,  -0.07579 , 0.0)
        H1 = gqcpy.Nucleus(1, 0.86681, 0.60144 , 0.0)
        H2 = gqcpy.Nucleus(1, -0.86681, 0.60144 , 0.0)
        water = gqcpy.Molecule([O,H1,H2])

        self.fukui_dyson_module = gqcpy.FukuiDysonAnalysis(water, "STO-3G", True)

        # Provide data generated with the C++ code
        self.fukui_vectors = np.array( [[8.86155e-17, -0.00707287, 0.998623, 0.0515051, -0.00704785, -7.01003e-17, 5.9811e-33],
                                        [7.39776e-15, -0.59471, 0.0168469, -0.494966, -0.633281, -9.07943e-15, 1.33927e-28],
                                        [-0.686138, -1.07911e-14, -9.61874e-17, 1.26337e-15, -9.42968e-15, 0.727472, 1.00769e-15],
                                        [-6.06117e-15, 0.345791, 0.0496853, -0.867346, 0.354501, 5.56009e-15, 6.80163e-29],
                                        [4.96848e-16, -1.9784e-29, -3.55624e-30, 6.16852e-29, -2.22645e-29, -3.90642e-16, 1.0],
                                        [2.22311e-14, -0.725741, 0.00013585, -0.00816175, 0.68792, 1.91338e-14, 5.81785e-31],
                                        [-0.727472, -2.09718e-14, -2.33172e-17, 6.91984e-16, 2.04771e-14, -0.686138, 3.90105e-17]])

        self.dyson_coefficients = np.array( [ 3.82513e-33, -9.30535e-31, 7.83798e-16, 7.61891e-29, 0.989971, -7.60146e-31, 1.3487e-17] )
        
        self.fukui_matrix = np.array(  [[6.54691e-07, 1.34911e-05, 1.48693e-19, -3.60703e-05, 8.29104e-33, -0.00118111, -1.98424e-17],
                                        [1.34911e-05, 0.00577327, -2.5339e-16, -0.00194855, -4.88262e-31, -0.102596, -1.49698e-15],
                                        [1.48693e-19, -2.5339e-16, 0.0077906, 3.14585e-17, 7.79157e-16, 2.14022e-15, -0.130879],
                                        [-3.60703e-05, -0.00194855, 3.14585e-17, 0.00122453, 7.68942e-29, 0.0585436, 6.49211e-16],
                                        [8.29104e-33, -4.88262e-31, 7.79157e-16, 7.68942e-29, 0.999057, -7.36496e-30, 1.34367e-16],
                                        [-0.00118111, -0.102596, 2.14022e-15, 0.0585436, -7.36496e-30, -0.00631616, 2.60986e-16],
                                        [-1.98424e-17, -1.49698e-15, -0.130879, 6.49211e-16, 1.34367e-16, 2.60986e-16, -0.00753]] )

        self.fukui_naturals = np.array( [-0.130973, -0.118294, -1.07302e-06, 0.000665593, 0.118312, 0.131233, 0.999057] )

        self.canonical_matrix = np.array( [[-0.993153, -0.235029, -5.27503e-17, -0.106599, 4.03116e-30, -0.157663, -2.30371e-15],
                                            [-0.0415334, 0.802792, 5.64077e-16, 0.50831, -1.98139e-29, 2.01176, 3.06243e-14],
                                            [3.96439e-17, 4.33903e-16, 0.681466, 1.20757e-15, -3.30016e-16, -2.17049e-14, 1.38215],
                                            [-0.0102528, 0.319206, 1.39008e-15, -0.878869, 3.40497e-29, 0.823904, 1.28786e-14],
                                            [1.10934e-34, -3.711e-31, 1.34971e-16, 3.72505e-29, 1.0, 1.24769e-30, -6.63605e-17],
                                            [0.0129601, 0.0985403, 0.39729, -0.0783755, 7.99344e-18, -1.30754, -1.56836],
                                            [0.0129601, 0.0985403, -0.39729, -0.0783755, -7.99344e-18, -1.30754, 1.56836]])


    def tearDown(self):
        pass


    def test_analysis(self):
        """Compare the various analysis parameters with their C++ reference values"""
        self.assertTrue(np.allclose(np.absolute(self.fukui_dyson_module.get_fukui_naturals()), np.absolute(self.fukui_naturals), rtol=1e-05, atol=1e-06))
        self.assertTrue(np.allclose(np.absolute(self.fukui_dyson_module.get_fukui_matrix()), np.absolute(self.fukui_matrix), rtol=1e-05, atol=1e-06))
        self.assertTrue(np.allclose(np.absolute(self.fukui_dyson_module.get_dyson_coefficients()), np.absolute(self.dyson_coefficients), rtol=1e-05, atol=1e-06))
        self.assertTrue(np.allclose(np.absolute(self.fukui_dyson_module.get_fukui_vectors()), np.absolute(self.fukui_vectors), rtol=1e-05, atol=1e-06))
        self.assertTrue(np.allclose(np.absolute(self.fukui_dyson_module.get_canonical_matrix()), np.absolute(self.canonical_matrix), rtol=1e-05, atol=1e-06))


if __name__ == '__main__':
    unittest.main()
