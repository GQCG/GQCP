from distutils.core import setup

import sys
if sys.version_info < (3,0):
    sys.exit('Sorry, Python < 3.0 is not supported')

setup(
    name='gqcpy',
    version='${PACKAGE_VERSION}',
    packages=['gqcpy'],
    package_dir={
        '': '${CMAKE_BINARY_DIR}'
    },
    package_data={
        '': ['gqcpy.so']
    }
)
