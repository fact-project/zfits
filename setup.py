from setuptools import setup, Extension
from Cython.Distutils import build_ext

import numpy

setup(
    name='zfits',
    version='0.0.3',
    description='a pure python zfits/factfits reader',
    url='https://github.com/fact-project/zfits',
    author='Dominik Neise',
    author_email='neised@phys.ethz.ch',
    license='MIT',
    packages=['zfits'],
    install_requires=[
        'fitsio>=23',
        'numpy',
    ],
    dependency_links=[
        'git+https://github.com/dneise/fitsio#egg=fitsio-23.0.1'
    ],
    entry_points={},
    package_data={'zfits': ['test_data/*']},
    zip_safe=False,
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension(
            'zfits.cython_tools',
            sources=['zfits/cython_tools.pyx', 'zfits/fact++wrapper.cpp'],
            include_dirs=[numpy.get_include(), 'zfits'],
            language="c++",
            extra_compile_args=['-std=c++0x']
        )
    ],
    tests_require=['pytest>=3.0.0'],
    setup_requires=['pytest-runner'],
)
