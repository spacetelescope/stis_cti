#! /usr/bin/env python

import archive_dark_query
from astropy.io import fits
from collections import defaultdict
import numpy as np
import datetime
import re
import ipdb

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
        hdr0['CCDGAIN'] in gains_allowed               and \
        hdr0['CCDOFFST'] in offsts_allowed


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
    #import ipdb
    
    # Get information about the user's system:
    try:
        import multiprocessing
        num_available_cores = multiprocessing.cpu_count()
        cores_str = "; number of available CPU cores on your system = " + str(num_available_cores)
    except (ImportError, NotImplementedError):
        num_available_cores = 1
        cores_str = ""
    
    # For prettier implementation of help text, see:
    # http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
    
    parser = argparse.ArgumentParser( 
        description='Run STIS/CCD pixel-based CTI-correction on data specified in SCIENCE_DIR. ' + \
                    'Uncorrected component darks are read from DARK_DIR, and ' + \
                    'corrected component darks are written there too. ' + \
                    'Corrected super-darks are read from and stored to REF_DIR.', 
        epilog='Author:   ' + __author__ + '; ' + \
               'Version:  ' + __version__)
    parser.add_argument(dest='science_dir', action='store', default='./', nargs='?', 
                        metavar='SCIENCE_DIR', 
                        help='directory containing RAW science data (default=\"./\")')
    parser.add_argument('-d', dest='dark_dir', action='store', default=None, 
                        help='directory of dark FLT data (default=\"[SCIENCE_DIR]/../darks/\")')
    parser.add_argument('-r', dest='ref_dir', action='store', default=None, 
                        help='directory of CTI-corrected reference files ' + \
                             '(default=\"[SCIENCE_DIR]/../ref/\")')
    parser.add_argument('-n', dest='num_processes', action='store', default=1, metavar='NUM_PROCESSES', 
                        help='maximum number of parallel processes to run (default=1)' + \
                              cores_str)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False, 
                        help='print more information')
    parser.add_argument('-vv', dest='very_verbose', action='store_true', 
                        help='Very verbose')
    # Allow any amp/gain/offst/date through processing (not well-tested):
    parser.add_argument('--allow', dest='allow', action='store_true', default=False, 
                        help=argparse.SUPPRESS)
    args = parser.parse_args()
    
    verbose = args.verbose
    
    if args.very_verbose:
        verbose = 2
    
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
    
    # Print script name/version, system date/time.
    # ...
    
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
                        pass
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
                print 'Superdark %s not found!  Will make a CTI-corrected version.' % superdark_resolved
            superdark_remakes.extend(file)
    
    if verbose:
        print
        if len(superdark_remakes) >= 1:
            print 'Making superdarks for:'
            print '   ' + '\n   '.join(filtered_raw_files)
        else:
            print 'Not remaking any superdarks.'
        print
    
    # Determine component darks used to make superdarks:
    #anneal_data = archive_dark_query.get_anneal_boundaries()  # *** Allow user-options here? ***
    with open('/Users/lockwood/stis_cte/wrapper/anneals.p', 'rb') as p:
        anneal_data = pickle.load(p)  # *** For testing only!!! ***
    print 'WARNING:  *** Using pickle file anneals.p for testing! ***'
    print
    anneals = archive_dark_query.archive_dark_query( \
                      filtered_raw_files, anneal_data=anneal_data, print_url=False)  # print_url?
    
    # Get list of EXPNAMEs from files in the dark_dir.
    found_dark_files = glob.glob(os.path.join(dark_dir, '*_flt.fits*'))  # *** Do something to allow RAW files? ***
    if verbose >= 2:
        max_file_length = str(max(map(lambda x: len(x), found_dark_files)))
    found = {}
    for file in found_dark_files:
        with fits.open(file) as f:
            hdr0 = f[0].header
        found[hdr0['ROOTNAME'].strip().upper()] = (file, hdr0)
    
    # Loop over expected darks within anneals and populate keyword data, noting missing files:
    get_keywords = ['CCDAMP', 'CCDGAIN', 'CCDOFFST', 'BIASFILE', 'DARKFILE', 'TELESCOP', 'INSTRUME', 'DETECTOR']
    missing_darks = set()
    weekdark_types = defaultdict(list)
    for anneal in anneals:
        # For each amp, count weeks based on old darkfile; reset for each annealing period:
        all_weeks = {'A':{}, 'B':{}, 'C':{}, 'D':{}}
        for dark in anneal['darks']:
            if found.has_key(dark['exposure']):
                dark['file'] = found[dark['exposure']][0]
                for get_keyword in get_keywords:
                    dark[get_keyword] = found[dark['exposure']][1][get_keyword]
                    try:
                        # Strip leading/trailing spaces off of strings:
                        dark[get_keyword] = dark[get_keyword].strip()
                    except AttributeError:
                        pass
                
                # Determine week number (for weekdark separation) for each amp:
                if len(all_weeks[dark['CCDAMP']]) == 0:
                    all_weeks[dark['CCDAMP']][dark['DARKFILE']] = 1
                elif dark['DARKFILE'] not in all_weeks[dark['CCDAMP']].keys():
                    all_weeks[dark['CCDAMP']][dark['DARKFILE']] = len(all_weeks[dark['CCDAMP']]) + 1
                dark['week_num'] = all_weeks[dark['CCDAMP']][dark['DARKFILE']]
                # Construct a unique tag for each weekdark and assign it to the component darks:
                dark['weekdark_tag'] = '{}{:03.0f}_{:1.0f}'.format(dark['CCDAMP'].lower(), anneal['index'], dark['week_num'])
                weekdark_types[dark['weekdark_tag']].append(dark)
                
                if verbose >= 2:
                    print ('{:s} {:5s} {:' + max_file_length + 's} {} {} {} {} {:5.0f} {} {} {}').format( \
                        dark['exposure'], dark['DETECTOR'], dark['file'], dark['CCDAMP'], \
                        dark['CCDGAIN'], dark['CCDOFFST'], dark['datetime'].isoformat(' '), \
                        dark['exptime'], dark['week_num'], dark['weekdark_tag'], dark['DARKFILE'])
            else:
                missing_darks.add(dark['exposure'])
    if verbose >= 2:
        print
    
    if len(missing_darks) != 0:
        missing_darks = list(missing_darks)
        missing_darks.sort()
        print 'ERROR:  These FLT component darks are missing from %s:' % dark_dir
        print ', '.join(missing_darks)
        print
        print 'Please download the missing darks via this link:'
        print '(or specify the proper dark_dir [' + dark_dir + '])'
        print
        print archive_dark_query.darks_url(missing_darks)
        print
        raise IOError('Missing component dark FLT files.')
    
    if verbose:
        print 'All required component dark FLT files for annealing periods have been located on disk.'
        print
    
    # Determine week_num time boundaries (approximate for now)
    # ...
    
    weekdarks = {}
    for weekdark_tag in weekdark_types:
        weekdarks[weekdark_tag] = {
            'amp'          : weekdark_tag[0].upper(), 
            'start'        : min([x['datetime'] for x in weekdark_types[weekdark_tag]]), 
            'end'          : max([x['datetime'] for x in weekdark_types[weekdark_tag]]), 
            'darks'        : list(weekdark_types[weekdark_tag]), 
            'anneal_num'   : int(weekdark_tag[1:4]), 
            'week_num'     : int(weekdark_tag[5]), 
            'weekdark_tag' : weekdark_tag }
    
    # Loop over each {amp, anneal_num} combination and adjust week boundaries:
    amps_anneals = list(set([(weekdarks[x]['amp'], weekdarks[x]['anneal_num']) for x in weekdarks]))
    for amp_anneal in amps_anneals:
        # Determine which data matches the (amp, anneal) specified in loop iteration:
        amp_anneal_data = [weekdarks[x] for x in weekdarks 
            if weekdarks[x]['amp'] == amp_anneal[0] and weekdarks[x]['anneal_num'] == amp_anneal[1]]
        weeks = list(amp_anneal_data)
        weeks.sort(key=lambda x: x['week_num'])
        # Extend start/end boundaries to anneal boundaries:
        anneal = filter(lambda x: weeks[0]['anneal_num'] == x['index'], anneals)[0]
        weekdarks[weeks[ 0]['weekdark_tag']]['start'] = anneal['start']  # start of annealing period
        weekdarks[weeks[-1]['weekdark_tag']]['end']   = anneal['end']    # end of annealing period
        # Make inter-week boundaries contiguous:
        if len(weeks) >= 2:
            for i, week in enumerate(weeks[0:-1]):
                delta = datetime.timedelta(seconds = round((weeks[i+1]['start'] - weeks[i]['end']).total_seconds() / 2.0))
                midpt = weeks[i]['end'] + delta
                weekdarks[weeks[i  ]['weekdark_tag']]['end']   = midpt
                weekdarks[weeks[i+1]['weekdark_tag']]['start'] = midpt
    
    # Determine which anneal_period/week_num each science file belongs to
    # ...
    
    # Filter anneals/darks based on science data quantities (time, CCDAMP, week_num, ...)
    # ...
    
    if verbose:
        # This should be after filtering by amp...
        print 'Could make weekdarks for:'
        for weekdark_tag in weekdarks:
            print '   {}:  {} - {}'.format(weekdark_tag, 
                weekdarks[weekdark_tag]['start'], 
                weekdarks[weekdark_tag]['end'] )
        print
    
    