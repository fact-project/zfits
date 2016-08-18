from setuptools import setup
from Cython.Build import cythonize

setup(
    name='zfits',
    version='0.0.1',
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
    ext_modules=cythonize("zfits/cython_tools.pyx"),
    entry_points={},
    package_data={'zfits': ['test_data/*']},
    zip_safe=False,
)