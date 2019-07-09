from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="PyTurbo",
    ext_modules = cythonize("PyTurbo.pyx"),
)
