import gqcpy
import unittest

class TestVersion(unittest.TestCase):

    def testMajorVersion(self):
        self.assertIsNotNone(gqcpy.Version.majorVersion())

    def testMinorVersion(self):
        self.assertIsNotNone(gqcpy.Version.minorVersion())

    def testPatchVersion(self):
        self.assertIsNotNone(gqcpy.Version.patchVersion())

    def testFullVersion(self):
        self.assertIsNotNone(gqcpy.Version.fullVersion())

    def testGitSHA1(self):
        self.assertIsNotNone(gqcpy.Version.gitSHA1())

if __name__ == '__main__':
    unittest.main()