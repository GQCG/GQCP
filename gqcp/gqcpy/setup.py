from distutils.core import setup

import sys
if sys.version_info < (3,0):
    sys.exit('Sorry, Python < 3.0 is not supported')

setup(
    name='gqcpy',
    version='',
    packages=['gqcpy'],
    package_dir={
        '': '/Users/daria/ugent/gqcp/gqcp'
    },
    package_data={
        '': ['gqcpy.so']
    }
)
