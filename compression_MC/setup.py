from distutils.core import setup
from Cython.Build import cythonize
import numpy
#import os

#os.environ['CFLAGS'] = '-O3'
setup(
    ext_modules=cythonize(['Capsid.pyx',        # Cython code file with primes() function
                           'monte_carlo.pyx',
                           'energy.pyx',
                           'geometry.pyx',
                           'helpers.pyx',
                           'parameters.py',
                           'data_process.py'],  # Python code file with primes_python_compiled() function
                          annotate=True),        # enables generation of the html annotation file
    include_dirs=[numpy.get_include()]
)
