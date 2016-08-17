from setuptools import setup

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
    entry_points={},
    zip_safe=False,
)