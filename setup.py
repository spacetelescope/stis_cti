#!/usr/bin/env python

from setuptools import setup, Extension
import os
import numpy

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst')) as f:
    long_description = f.read()

StisPixCteCorr_module = Extension('stis_cti.StisPixCte_FixY',
    sources=['src/StisFixYCte.c', 'src/StisPixCteCorr_funcs.c', 'src/StisPixCte_FixY.c'],
    include_dirs=[numpy.get_include(), 'src/'],
    define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],)

setup(
    name = 'stis_cti',
    url = 'https://www.stsci.edu/hst/instrumentation/stis/data-analysis-and-software-tools/pixel-based-cti',
    project_urls={
        'Homepage': 'https://www.stsci.edu/hst/instrumentation/stis/data-analysis-and-software-tools/pixel-based-cti',
        'Documentation': 'https://stis-cti.readthedocs.io',
        'Help Desk': 'https://hsthelp.stsci.edu',
        'Source Code': 'https://github.com/spacetelescope/stis_cti',
        'Issues': 'https://github.com/spacetelescope/stis_cti/issues',},
    version = '1.6.1',
    description = 'Pixel-based CTI-correction for HST/STIS CCD data',
    long_description = long_description,
    author = 'Sean Lockwood, Phil Hodge, Pey Lian Lim, W.J. Hack, J. Anderson, Matt Davis',
    author_email = 'lockwood@stsci.edu',
    maintainer = 'Sean Lockwood',
    maintainer_email = 'lockwood@stsci.edu',
    license = 'BSD-new',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python :: 3',
                   'Programming Language :: C',
                   'Development Status :: 5 - Production/Stable',
                   'Operating System :: OS Independent',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: BSD License',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = ['stis_cti'],
    install_requires = ['setuptools', 'numpy', 'astropy>= 4.0', 'stistools>= 1.2', 
        'refstis>= 0.8.1', 'crds', 'stsci.tools>= 3.2.2', 'six'],
    entry_points = {'console_scripts': [
        'stis_cti = stis_cti.stis_cti:call_stis_cti',
        'archive_dark_query = stis_cti.archive_dark_query:call_archive_dark_query', ]},
    package_data = {'stis_cti': ['data/*.fits']},
    ext_modules = [StisPixCteCorr_module],
    )
