from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    name='zfits',
    version='0.0.2',
    description='a pure python zfits/factfits reader',
    url='https://github.com/fact-project/zfits',
    author='Dominik Neise',
    author_email='neised@phys.ethz.ch',
    license='MIT',
    packages=['zfits'],
    install_requires=[
        'fitsio',
        'numpy',
    ],
    entry_points={},
    package_data={'zfits': ['test_data/*']},
    zip_safe=False,
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("zfits.cython_tools",
                             sources=["zfits/cython_tools.pyx", "zfits/example.cpp"],
                             include_dirs=[numpy.get_include(), "zfits"],
                             language="c++",
                             extra_compile_args=['-std=c++0x'])],
)
