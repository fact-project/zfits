from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


def get_numpy_includes():
    import numpy
    return numpy.get_include()

setup(
    name='zfits',
    version='0.1.0',
    description='a pure python zfits/factfits reader',
    url='https://github.com/fact-project/zfits',
    author='Dominik Neise',
    author_email='neised@phys.ethz.ch',
    license='MIT',
    packages=['zfits'],
    install_requires=[
        'numpy>=1.12.1',
        'Cython>=0.25.2',
        'fitsio>=0.9.11',  # for `.headers()` and `drs.fits.gz` file.
        'pyfact>=0.12.1',
    ],
    entry_points={},
    package_data={'zfits': ['test_data/*']},
    zip_safe=False,
    ext_modules=cythonize([
        Extension(
            name="*",
            sources=["zfits/*.pyx"],
            extra_compile_args=['-std=c++0x'],
            language='c++',
            include_dirs=[get_numpy_includes(), 'zfits'],
        )
        ]),
    tests_require=['pytest>=3.0.0'],
    setup_requires=[
        'numpy',
        'pytest-runner'],
)
