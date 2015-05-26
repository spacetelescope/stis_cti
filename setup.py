#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy

StisPixCteCorr_module = Extension('stis_cti.StisPixCte_FixY',
    sources=['src/StisFixYCte.c', 'src/StisPixCteCorr_funcs.c', 'src/StisPixCte_FixY.c'],
    include_dirs=[numpy.get_include(), 'src/'])

setup(
    name = 'stis_cti',
    url = 'https://grit.stsci.edu/lockwood/cti_wrapper/repository/archive',
    version = '0.3_alpha',
    description = 'Apply pixel-based CTI-correction to HST/STIS CCD data',
    author = 'Sean Lockwood, Phil Hodge, Pey Lian Lim and ' + \
             'W.J. Hack (Python), J. Anderson (Fortran), Matt Davis',
    author_email = 'help@stsci.edu',
    maintainer = 'Sean Lockwood',
    maintainer_email = 'lockwood@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Development Status :: 2 - Pre-Alpha',
                   'Operating System :: OS Independent',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = ['stis_cti'],
    requires = ['numpy', 'astropy', 'refstis'],
    scripts =  ['scripts/stis_cti', 'scripts/archive_dark_query'],
    ext_modules = [StisPixCteCorr_module],
    license = 'TBD -- all rights reserved'
    )
