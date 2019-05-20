from setuptools import setup, find_packages
import os

# make sure users without cython can install our extensions
try:
    from Cython.Distutils.extension import Extension
    from Cython.Distutils import build_ext as _build_ext
    USE_CYTHON = True
except ImportError:
    from setuptools import Extension
    from setuptools.command.build_ext import build_ext as _build_ext
    USE_CYTHON = False


print('using cython', USE_CYTHON)


# make sure numpy is installed before we try to build
# the extenion
class build_ext(_build_ext):
    def finalize_options(self):
        super().finalize_options()
        import numpy
        self.include_dirs.append(numpy.get_include())


cmdclass = {'build_ext': build_ext}
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [
    Extension(
        name='zfits.' + name,
        sources=[os.path.join('zfits', name + ext)],
        extra_compile_args=['-std=c++0x'],
        language='c++',
        include_dirs=['zfits'],
    )
    for name in ['factfits', 'remove_spikes']
]


setup(
    name='zfits',
    version='0.2.0',
    description='a pure python zfits/factfits reader',
    url='https://github.com/fact-project/zfits',
    author='Dominik Neise',
    author_email='neised@phys.ethz.ch',
    license='MIT',
    packages=find_packages(),
    ext_modules=extensions,
    cmdclass=cmdclass,
    install_requires=[
        'numpy>=1.12.1',
        'fitsio>=0.9.11',  # for `.headers()` and `drs.fits.gz` file.
        'pyfact>=0.12.1',
    ],
    setup_requires=[
        'numpy',
    ],
    tests_require=['pytest>=3.0.0'],
    package_data={'zfits': ['test_data/*', '*.cpp']},
    zip_safe=False,
)
