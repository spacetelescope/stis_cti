#! /usr/bin/env python

import archive_dark_query
from astropy.io import fits
import datetime
import re

__author__  = 'Sean Lockwood'
__version__ = '0.1'


def viable_ccd_file(file, 
                    earliest_date_allowed=None, 
                    amplifiers_allowed=None, 
                    gains_allowed=None, 
                    offsts_allowed=None):
    
    # Set defaults:
    if earliest_date_allowed is None:
        earliest_date_allowed = datetime.datetime(2009, 5, 1, 0, 0, 0)  # *** Or when? ***
    if type(earliest_date_allowed) is not datetime.datetime:
        raise TypeError('earliest_date_allowed must be a datetime.datetime, not a %s' % type(earliest_date_allowed))
    
    if amplifiers_allowed is None:
        amplifiers_allowed = ['D']
    
    if gains_allowed is None:
        gains_allowed = [1]  # *** Or include 4? ***
    
    if offsts_allowed is None:
        offsts_allowed = [3]
    
    # Read ext=0 header of observation:
    with fits.open(file) as f:
        hdr0 = f[0].header
    
    # Parse date/time of observation:
    date = hdr0['TDATEOBS'].strip()
    time = hdr0['TTIMEOBS'].strip()
    dt = datetime.datetime.strptime(date + ' ' + time, '%Y-%m-%d %H:%M:%S')
    
    return \
        hdr0['TELESCOP'].strip() == 'HST'              and \
        hdr0['INSTRUME'].strip() == 'STIS'             and \
        hdr0['DETECTOR'].strip() == 'CCD'              and \
        dt >= earliest_date_allowed                    and \
        'ACQ' not in hdr0['OBSMODE'].strip()           and \
        hdr0['CCDAMP'].strip() in amplifiers_allowed   and \
        (hdr0['CCDGAIN'] in gains_allowed  and             \
         hdr0['CCDOFFST'] in offsts_allowed)            or \
        hdr0['TARGNAME'].strip() == 'CCDFLAT'


def resolve_iraf_file(file):
    # Email sent to phil to get this routine...
    
    dir = ''
    rootname = file
    
    if '$' in file and file[0] == '$':
        dir, rootname = file[1:].split('/', 1)  # What about '\' ?
    elif '$' in file:
        dir, rootname = file.split('$', 1)
    
    if dir != '':
        dir_resolved = os.getenv(dir)
        if dir_resolved is None:
            raise IOError('Can\'t resolve environmental variable \'%s\'.' % dir)
        else:
            dir = dir_resolved
    
    new_file = os.path.normpath(os.path.join(dir, rootname))
    
    return new_file


if __name__ == '__main__':
    import pickle
    import argparse
    import os
    import glob
    #import sys
    #import pdb
    
    # For prettier implementation of help text, see:
    # http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
    
    parser = argparse.ArgumentParser( 
        description='Run STIS/CCD pixel-based CTI-correction on data specified in science_dir.', 
        epilog='Author:   ' + __author__ + '; ' + 
               'Version:  ' + __version__)
    parser.add_argument(dest='science_dir', action='store', default='./', nargs='?', 
                        metavar='SCIENCE_DIR', 
                        help='directory containing raw science data (default = \"./\")')
    parser.add_argument('-d', dest='dark_dir', action='store', default=None, 
                        help='directory of dark data (default = \"science_dir/../darks/\")')
    parser.add_argument('-r', dest='ref_dir', action='store', default=None, 
                        help='directory of CTI-corrected reference files ' + \
                             '(default = \"science_dir/../ref/\")')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False, 
                        help='Print more information')
    #parser.add_argument('-vv', dest='very_verbose', action='store_true', 
    #                    help='Very verbose')
    parser.add_argument('--allow', dest='allow', action='store_true', default=False, 
                        help=argparse.SUPPRESS)  # Allow any file into reduction
    args = parser.parse_args()
    
    verbose = args.verbose
    
    #if args.very_verbose:
    #    verbose = 2
    
    # Determine default dark_dir and ref_dir relative to science_dir:
    if args.dark_dir is None:
        args.dark_dir = os.path.join(args.science_dir, '../darks')
    if args.ref_dir is None:
        args.ref_dir = os.path.join(args.science_dir, '../ref')
    
    # Normalize redundant relative path variables:
    science_dir = os.path.normpath(args.science_dir) + os.path.sep
    dark_dir    = os.path.normpath(args.dark_dir)    + os.path.sep
    ref_dir     = os.path.normpath(args.ref_dir)     + os.path.sep
    
    # Check that directories exist:
    # (*** Convert to 'raise IOError("...")'? ***)
    assert os.path.isdir(science_dir), 'science_dir does not exist:  ' + science_dir
    assert os.path.isdir(dark_dir),    'dark_dir does not exist:  '    + dark_dir
    assert os.path.isdir(ref_dir),     'ref_dir does not exist:  '     + ref_dir
    
    if verbose:
        print 'science_dir = ', science_dir
        print 'dark_dir    = ', dark_dir
        print 'ref_dir     = ', ref_dir
        print
    
    # Find science rootnames:
    raw_files = glob.glob(os.path.join(science_dir, '*_raw.fits*'))
    if verbose:
        print 'Input _raw.fits files:'
        print '   ' + '\n   '.join(raw_files)
        print
    
    if args.allow:
        # Allow any HST/STIS non-acqs:
        if verbose:
            print 'Leniently filtering input _raw.fits files...'
        filtered_raw_files = [file for file in raw_files if \
            viable_ccd_file( file, \
                earliest_date_allowed = datetime.datetime(1990,1,1,0,0,0), \
                amplifiers_allowed = ['A','B','C','D'], \
                gains_allowed = [1,2,4,8], \
                offsts_allowed = range(9))]
    else:
        # Normal filtering on files:
        filtered_raw_files = [file for file in raw_files if viable_ccd_file(file)]
    
    # *** What about CCDFLATs? ***
    # e.g. lockwood:/Users/lockwood/frontend/stis_12512/obpz06020_raw.fits
    #     {CCDAMP=D, CCDGAIN=4, CCDOFFST=3} with valid science data
    #     Allowed for now.
    
    not_correcting = list(set(raw_files) - set(filtered_raw_files))
    if len(not_correcting) > 0:
        not_correcting.sort()
        print 'WARNING:  Not running the correction on these files:'
        print '   ' + '\n   '.join(not_correcting)
        print
    
    # Make sure some files exist:
    if len(filtered_raw_files) == 0:
        raise IOError('No viable STIS CCD files were present in the directory!\n' + \
                       '   ' + science_dir + '\n'                                 + \
                       'Please make sure science data are:\n'                     + \
                       '   -- post-May 01 2009\n'                                 + \
                       '   -- CCDAMP = D\n'                                       + \
                       '   -- CCDGAIN = 1\n'                                      + \
                       '   -- Not OBSMODE = ACQ*')
    
    if verbose:
        print 'Input _raw.fits files being corrected:'
        print '   ' + '\n   '.join(filtered_raw_files)
        print
    
    # Check for original MAST data results:
    # ...
    
    # Check for results from previous run:
    # ...
    
    # Check science files for uncorrected super-darks:
    superdark_remakes = []
    for file in filtered_raw_files:
        with fits.open(file) as f:
            hdr0 = f[0].header
        superdark = hdr0['DARKFILE']
        superdark_resolved = resolve_iraf_file(superdark)
        #if verbose:
        #    print superdark, ' = ', superdark_resolved
        try:
            with fits.open(superdark_resolved) as sd:
                hdr0 = sd[0].header
                try:
                    if hdr0['PCTECORR'].strip().upper() == 'PERFORMED':  # *** Or whatever keyword/value combination we decide. ***
                        pass # filtered_raw_files.pop(0)  # *** TEST THIS!!! ***
                    else:
                        if verbose:
                            print 'Superdark %s is not CTI-corrected.' % superdark_resolved
                        superdark_remakes.extend(file)
                except KeyError:
                    if verbose:
                        print 'Superdark %s is not CTI-corrected.' % superdark_resolved
                    superdark_remakes.extend(file)
        except IOError:
            if verbose:
                print 'Superdark %s not found!' % superdark_resolved
            superdark_remakes.extend(file)
    
    if verbose:
        print
        print 'Making superdarks for:'
        print '   ' + '\n   '.join(filtered_raw_files)
        print
    
    # Determine component darks used to make superdarks:
    anneal_data = archive_dark_query.get_anneal_boundaries()  # *** Allow user-options here? ***
    anneals = archive_dark_query.archive_dark_query( \
                      filtered_raw_files, anneal_data=anneal_data, print_url=verbose)  # print_url?
    # Set of unique component darks at expected dark_dir/:
    expected_dark_expnames = set()
    for anneal in anneals:
        for dark in anneal.darks:
            expected_dark_expnames.add(dark.exposure)
    #expected_dark_expnames = list(expected_dark_expnames)
    #expected_dark_expnames.sort()
    
    # Get list of EXPNAMEs from files in the dark_dir.
    found_dark_files = glob.glob(os.path.join(dark_dir, '*_raw.fits*'))
    # *** Filter darks through viable_ccd_file()? ***
    # What about expected expnames that are excluded due to {CCDAMP, CCDGAIN, ...}?
    
    found_dark_expnames = set()
    for file in found_dark_files:
        with fits.open(file) as f:
            hdr0 = f[0].header
        found_dark_expnames.add(hdr0['ROOTNAME'].strip().upper())
    #found_dark_expnames = list(found_dark_expnames)
    #found_dark_expnames.sort()
    # *** Add error handling! ***
    
    # Compare found_dark_expnames and expected_dark_expnames:
    if verbose:
        print 'found_dark_expnames:'
        print ', '.join(found_dark_expnames)
        print
        print 'expected_dark_expnames:'
        print ', '.join(expected_dark_expnames)
        print
    
    missing_darks = expected_dark_expnames - found_dark_expnames
    print 'ERROR:  These raw component darks are missing from %s:' % dark_dir
    print ', '.join(missing_darks)
    print
    # print useful URL to download darks and throw an exception...
    
    
