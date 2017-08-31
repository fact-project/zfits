from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='zfits',
    version='0.0.4',
    description='a pure python zfits/factfits reader',
    url='https://github.com/fact-project/zfits',
    author='Dominik Neise',
    author_email='neised@phys.ethz.ch',
    license='MIT',
    packages=['zfits'],
    install_requires=[
        'numpy',
    ],
    entry_points={},
    package_data={'zfits': ['test_data/*']},
    zip_safe=False,
    ext_modules=cythonize(
           "zfits/*.pyx",                 # our Cython source
    ),
    tests_require=['pytest>=3.0.0'],
    setup_requires=['pytest-runner'],
)
