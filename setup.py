#!/usr/bin/env python

from setuptools import setup, Extension
import numpy

StisPixCteCorr_module = Extension('stis_cti.StisPixCte_FixY',
    sources=['src/StisFixYCte.c', 'src/StisPixCteCorr_funcs.c', 'src/StisPixCte_FixY.c'],
    include_dirs=[numpy.get_include(), 'src/'])

setup(
    name = 'stis_cti',
    url = 'https://grit.stsci.edu/lockwood/cti_wrapper/repository/archive',
    version = '0.4_beta3',
    description = 'Apply pixel-based CTI-correction to HST/STIS CCD data',
    author = 'Sean Lockwood, Phil Hodge, Pey Lian Lim and ' + \
             'W.J. Hack (Python), J. Anderson (Fortran), Matt Davis',
    author_email = 'help@stsci.edu',
    maintainer = 'Sean Lockwood',
    maintainer_email = 'lockwood@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Development Status :: 4 - Beta',
                   'Operating System :: OS Independent',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: BSD License',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = ['stis_cti'],
    install_requires = ['setuptools', 'numpy', 'astropy>= 1.0.1', 'stistools>= 1.0.2', 
        'refstis', 'crds>= 1.3.0', 'stistools>= 1.0.2', 'stsci.tools>= 3.2.2'],
    scripts = ['scripts/stis_cti', 'scripts/archive_dark_query'],
    package_data = {'stis_cti': ['data/*.fits']},
    ext_modules = [StisPixCteCorr_module],
    )
