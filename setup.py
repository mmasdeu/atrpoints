import os.path

from Cython.Build import cythonize
from setuptools import setup

setup(
    ext_modules=cythonize(
        os.path.join("atrpoints", "*.pyx"),
        compiler_directives={"language_level": 3},
    )
)
