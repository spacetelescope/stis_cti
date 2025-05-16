#!/usr/bin/env python

from setuptools import setup, Extension
import os
import numpy

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

StisPixCteCorr_module = Extension(
    'stis_cti.StisPixCte_FixY',
    sources=[
        os.path.join('src', 'StisFixYCte.c'),
        os.path.join('src', 'StisPixCteCorr_funcs.c'),
        os.path.join('src', 'StisPixCte_FixY.c'),
    ],
    include_dirs=[numpy.get_include(), 'src/'],
    define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
)

setup(
    name='stis_cti',
    use_scm_version={'write_to': 'stis_cti/_version.py'},
    setup_requires=['setuptools>=42', 'setuptools-scm'],
    description='Pixel-based CTI-correction for HST/STIS CCD data',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    author='Sean Lockwood, Phil Hodge, Pey Lian Lim, W.J. Hack, J. Anderson, Matt Davis',
    author_email='lockwood@stsci.edu',
    maintainer='Sean Lockwood',
    maintainer_email='lockwood@stsci.edu',
    license='BSD-new',
    url='https://www.stsci.edu/hst/instrumentation/stis/data-analysis-and-software-tools/pixel-based-cti',
    project_urls={
        'Homepage': 'https://www.stsci.edu/hst/instrumentation/stis/data-analysis-and-software-tools/pixel-based-cti',
        'Documentation': 'https://stis-cti.readthedocs.io',
        'Help Desk': 'https://hsthelp.stsci.edu',
        'Source Code': 'https://github.com/spacetelescope/stis_cti',
        'Issues': 'https://github.com/spacetelescope/stis_cti/issues',
    },
    keywords=['astronomy'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
        'Development Status :: 5 - Production/Stable',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages=['stis_cti'],
    install_requires=[
        'setuptools',
        'numpy',
        'astropy>=4.0',
        'stistools>=1.2',
        'refstis>=0.8.1',
        'crds',
        'stsci.tools>=3.2.2',
        'six',
    ],
    extras_require={"docs": ["sphinx",]},
    entry_points={
        'console_scripts': [
            'stis_cti = stis_cti.stis_cti:call_stis_cti',
            'archive_dark_query = stis_cti.archive_dark_query:call_archive_dark_query',
        ]
    },
    package_data={'stis_cti': ['data/*.fits']},
    ext_modules=[StisPixCteCorr_module],
)
