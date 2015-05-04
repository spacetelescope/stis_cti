#! /usr/bin/env python

import os
import sys
import multiprocessing
import datetime
from astropy.io import fits
from collections import defaultdict
import archive_dark_query
import refstis
from stistools import StisPixCteCorr, basic2d
#try:
#    import ipdb as pdb
#except ImportError:
#    import pdb

__author__  = 'Sean Lockwood'
__version__ = '0.1'


# cti_wrapper()
# determine_input_science()
# viable_ccd_file()
# resolve_iraf_file()
# superdark_hash()
# perform_cti_correction()
# copy_dark_keywords()
# generate_basedark()
# generate_weekdark()
# populate_darkfiles()
# class Logger
# main


def cti_wrapper(science_dir, dark_dir, ref_dir, pctetab, num_processes, 
                all_weeks_flag=False, allow=False, verbose=False):
    '''
    Run STIS/CCD pixel-based CTI-correction on data specified in SCIENCE_DIR.
    
    Uncorrected component darks are read from DARK_DIR, and corrected component
    darks are written there too. Corrected super-darks are read from and stored to
    REF_DIR.
    
    From the command line:
        usage: cti_wrapper.py [-h] [-d DARK_DIR] [-r REF_DIR] [-n NUM_PROCESSES]
                              [-p PCTETAB] [-v | -vv]
                              [SCIENCE_DIR]
        
        positional arguments:
          SCIENCE_DIR       directory containing RAW science data (default="./")
        
        optional arguments:
          -h, --help        show this help message and exit
          -d DARK_DIR       directory of dark FLT data
                            (default="[SCIENCE_DIR]/../darks/")
          -r REF_DIR        directory of CTI-corrected reference files
                            (default="[SCIENCE_DIR]/../ref/")
          -n NUM_PROCESSES  maximum number of parallel processes to run (default=X);
                            number of available CPU cores on your system = Y
          -p PCTETAB        name of PCTETAB to use in pixel-based correction
                            (default="[REF_DIR]/test_pcte.fits")
          -v, --verbose     print more information
          -vv               very verbose
    
    Note that 'all_week_flag' and 'allow' options have not been fully tested!
    '''
    # Open a log file and copy {STDOUT, STDERR} to it:
    log = Logger(os.path.join(science_dir, 'cti_{}.log'.format(datetime.datetime.now().isoformat('_'))))
    
    try:
        from platform import node, platform
        sys_info = '{}; {}'.format(node(), platform())
    except ImportError:
        sys_info = 'Indeterminate'
    
    # Print system information:
    print 'Running CTI-correction script:  {} v{}'.format(os.path.basename(__file__), __version__)
    print 'System:                         {}'.format(sys_info)
    print 'Start time:                     {}\n'.format(datetime.datetime.now().isoformat(' '))
    
    if verbose:
        print 'science_dir = {}'.format(science_dir)
        print 'dark_dir    = {}'.format(dark_dir)
        print 'ref_dir     = {}'.format(ref_dir)
        print 'PCTETAB     = {}\n'.format(pctetab)
    log.flush()
    
    raw_files = determine_input_science(science_dir, allow, verbose)
    log.flush()
    
    # Check for original MAST data results:
    # ...
    
    # Check for results from previous run:
    # ...
    
    # Check science files for uncorrected super-darks:
    populate_darkfiles(raw_files, dark_dir, ref_dir, num_processes, all_weeks_flag, verbose)
    
    bias_corrected = bias_correct_science_files(raw_files, verbose)
    
    log.close()


def determine_input_science(science_dir, allow=False, verbose=False):
    import glob
    
    # Find science rootnames:
    raw_files = glob.glob(os.path.join(science_dir, '*_raw.fits*'))
    if verbose:
        print 'Input _raw.fits files:'
        print '   ' + '\n   '.join(raw_files) + '\n'
    
    if allow:
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
        print '   ' + '\n   '.join(not_correcting) + '\n'
    
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
        print '   ' + '\n   '.join(filtered_raw_files) + '\n'
    
    return filtered_raw_files


def viable_ccd_file(file, 
                    earliest_date_allowed=None, 
                    amplifiers_allowed=None, 
                    gains_allowed=None, 
                    offsts_allowed=None):
    
    # Set defaults:
    if earliest_date_allowed is None:
        earliest_date_allowed = datetime.datetime(2009, 5, 1, 0, 0, 0)  # Or when? ***
    if type(earliest_date_allowed) is not datetime.datetime:
        raise TypeError('earliest_date_allowed must be a datetime.datetime, not a %s' % type(earliest_date_allowed))
    
    if amplifiers_allowed is None:
        amplifiers_allowed = ['D']
    
    if gains_allowed is None:
        gains_allowed = [1]  # Or include 4? ***
    
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


def bias_correct_science_files(raw_files, verbose):
    if verbose:
        tmp = basic2d.basic2d('', print_revision=True)
    
    outnames = [f.replace('_raw.fits', '_blt.fits', 1) for f in raw_files]
    
    # Check for previous _blt.fits files first:
    for file in outnames:
        if os.path.exists(outname):
            raise IOError('File {} already exists!'.format(outname))
    
    for raw_file, outname in zip(raw_files, outnames):
        if verbose:
            print 'Running basic2d on {} --> {}.'.format(raw_file, outname)
        status = basic2d.basic2d(raw_file, output=outname, dqicorr=True, 
            blevcorr=True, biascorr=True, doppcorr=False, lorscorr=False, glincorr=False, 
            lflgcorr=False, darkcorr=False, flatcorr=False, photcorr=False, statflag=True, 
            verbose=(verbose >= 2))
        if status != 0:
            raise RuntimeError('basic2d returned non-zero status on {}:  {}'.format(raw_file, status))


def resolve_iraf_file(file):
    # Email sent to phil to get the proper version of this routine...
    
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


def superdark_hash(sim_nit=None, shft_nit=None, rn_clip=None, nsemodel=None, subthrsh=None,
                            pctetab=None,
                            superdark=None, files=None):
    from numpy import where
    
    if pctetab is not None:
        with fits.open(pctetab) as f:
            sim_nit  = f[0].header['SIM_NIT']
            shft_nit = f[0].header['SHFT_NIT']
            rn_clip  = f[0].header['RN_CLIP']
            nsemodel = f[0].header['NSEMODEL']
            subthrsh = f[0].header['SUBTHRSH']
        if files is None:
            raise IOError('Must specify files with pctetab.')
    elif superdark is not None:
        with fits.open(superdark) as f:
            sim_nit  = f[0].header['PCTESMIT']
            shft_nit = f[0].header['PCTESHFT']
            rn_clip  = f[0].header['PCTERNCL']
            nsemodel = f[0].header['PCTENSMD']
            subthrsh = f[0].header['PCTETRSH']
            if files is None:
                # Parse superdark header for files used:
                hist = f[0].header['HISTORY']
                beginning = where(['The following input files were used:' in line for line in hist])[0][0] + 1
                end = where([line == '' for line in hist[beginning:]])[0][0] + beginning
                files = hist[beginning:end]
    else:
        if sim_nit  == None or \
           shft_nit == None or \
           rn_clip  == None or \
           nsemodel == None or \
           subthrsh == None:
            raise IOError('Please specify either pctetab or the proper set of parameters!')
        if files is None:
            raise IOError('Must specify files with pctetab.')
    
    # Turn list of filenames into a sorted list of exposures:
    exposures = [os.path.basename(file).rsplit('_',1)[0].upper() for file in files]
    exposures.sort()
    exposures = ','.join(exposures)
    
    hash_str = '{};{};{};{};{};{}'.format( \
        exposures, \
        sim_nit, shft_nit, rn_clip, nsemodel, subthrsh)
    
    #return hash_str
    return hash(hash_str)


def perform_cti_correction(files, pctetab, num_cpu=1, verbose=False):
    # The call to StisPixCteCorr should be done in parallel!
    
    perform_files = []
    outnames = []
    for file in files:
        outname = file.replace('_flt.fits', '_cte.fits', 1)  # Or whatever intermediate extension is designated ***
        #outname = file.replace('_raw.fits', '_cte.fits', 1) # Needed as well? ***
        outnames.append(outname)
        if os.path.exists(outname):
            if superdark_hash(pctetab=pctetab, files=[]) == superdark_hash(superdark=outname, files=[]):
                if verbose:
                    print 'Skipping regeneration of CTI-corrected component dark:  {}'.format(outname)
        else:
            with fits.open(file, 'update') as f:
                #Update PCTETAB header keyword:
                try:
                    old_pctetab = f[0].header['PCTETAB']
                except KeyError:
                    old_pctetab = None
                
                if old_pctetab is None:
                    f[0].header.insert('DARKFILE', ('PCTETAB', pctetab))  # Insert after DARKFILE
                else:
                    f[0].header['PCTETAB'] = pctetab
                f.flush()
                
                if verbose:
                    print 'Updated hdr0 PCTETAB  of {}:  {} --> {}'.format(file, old_pctetab, pctetab)
                
                # Run the pixel-based correction on these component darks:
                perform_files.append(file)
    
    # Run the CTI-correction:
    p = multiprocessing.Pool(processes = num_cpu)
    p.map_async(StisPixCteCorr.CteCorr, perform_files)
    p.close()
    p.join()
    
    # Single-threaded version:
    #StisPixCteCorr.CteCorr(perform_files)
    
    # *** Do the corrected files need to be fed through DQICORR again to fix flags? ***
    
    if len(outnames) == 0:
        outnames = None
    elif len(outnames) == 1:
        outnames = outnames[0]
    return outnames


def copy_dark_keywords(superdark, dark_hdr0, history=None, basedark=None):
    # Copy these keywords from the last component dark to the new superdark:
    keywords = ['PCTECORR', 'PCTETAB', 'PCTEFRAC', 'PCTERNCL', 'PCTERNCL', 'PCTENSMD', 
                'PCTETRSH', 'PCTESMIT', 'PCTESHFT', 'CTE_NAME', 'CTE_VER']
    # Note:  PCTEFRAC is time-dependent; PCTECOR = COMPLETE
    keywords.reverse()
    
    with fits.open(superdark, 'update') as s:
        for keyword in keywords:
            value = dark_hdr0.get(keyword, default='unknown')
            try:
                comment = dark_hdr0.comments[keyword]
            except KeyError:
                comment = None
            
            if keyword in s[0].header:
                s[0].header[keyword] = value
            else:
                s[0].header.insert('DRK_VS_T', (keyword, value, comment), after=True)
            s.flush()
        
        if basedark is not None:
            if 'BASEDARK' in s[0].header:
                s[0].header['BASEDARK'] = basedark
            else:
                s[0].header.insert(keywords[0], ('BASEDARK', basedark, 'Used to make weekdark'), after=True)
        
        # Add HISTORY to superdark hdr0 here:
        if history is not None:
            s[0].header['HISTORY'] = ' '
            for h in history:
                s[0].header['HISTORY'] = h


def generate_basedark(files, outname, pctetab, num_cpu, verbose=False):
    # Don't make a basedark if it already exists:
    if os.path.exists(outname):
        if superdark_hash(pctetab=pctetab, files=files) == superdark_hash(superdark=outname):
            if verbose:
                print 'Skipping regeneration of basedark:  {}'.format(outname)
            return
    
    if verbose >= 2:
        print 'Working on basedark {}:'.format(outname)
        print '   ' + '\n   '.join(files) + '\n'
    
    # Correct files component darks, if necessary:
    corrected_files = perform_cti_correction(files, pctetab, num_cpu, verbose)
    
    # Make a basedark from the corrected darks:
    refstis.basedark.make_basedark(corrected_files, refdark_name=outname)
    # Print results of "[REF_DIR]/[basedark_name - '.fits']_joined_bd_calstis_log.txt"? ***
    
    if verbose:
        calstis_log = outname.replace('.fits','_joined_bd_calstis_log.txt', 1)
        print 'Calstis log file for CRJ processing [{}]:'.format(calstis_log)
        with open(calstis_log, 'r') as cs_log:
            for line in cs_log.readlines():
                print '     ' + line.strip()
        print
    #os.remove(calstis_log)  # ***
    
    # Copy the last file's ext=0 header into a variable to use in populating the basedark header:
    # (Maybe we should read all headers and make sure the keywords are the same [except PCTEFRAC]? ***)
    with fits.open(corrected_files[-1]) as file:
        dark_hdr0 = file[0].header
    
    # Update keywords in new basedark from component dark:
    history = [
        'Basedark created from CTI-corrected component darks by script',
        'cti_wrapper.py on {}.'.format(datetime.datetime.now().isoformat(' ')) ]
    copy_dark_keywords(outname, dark_hdr0, history=history)
        
    if verbose:
        print 'Basedark complete:  {}\n'.format(outname)


def generate_weekdark(files, outname, pctetab, basedark, num_cpu, verbose=False):
    # Don't make a weekdark if it already exists:
    if os.path.exists(outname):
        if superdark_hash(pctetab=pctetab, files=files) == superdark_hash(superdark=outname):
            if verbose:
                print 'Skipping regeneration of weekdark:  {}'.format(outname)
            return
    
    if verbose >= 2:
        print 'Working on weekdark {}:'.format(outname)
        print '   ' + '\n   '.join(files) + '\n'
    
    # Correct files component darks, if necessary:
    corrected_files = perform_cti_correction(files, pctetab, num_cpu, verbose)
    
    # Make a weekdark from the corrected darks:
    refstis.weekdark.make_weekdark(corrected_files, outname, basedark)
    
    if verbose:
        calstis_log = outname.replace('.fits','_joined_bd_calstis_log.txt', 1)
        print 'Calstis log file for CRJ processing [{}]:'.format(calstis_log)
        with open(calstis_log, 'r') as cs_log:
            for line in cs_log.readlines():
                print '     ' + line.strip()
        print
    #os.remove(calstis_log)  # ***
    
    # Copy the last file's ext=0 header into a variable to use in populating the basedark header:
    # (Maybe we should read all headers and make sure the keywords are the same [except PCTEFRAC]? ***)
    with fits.open(corrected_files[-1]) as file:
        dark_hdr0 = file[0].header
    
    # Update keywords in new basedark from component dark:
    history = [
        'Weekdark created from CTI-corrected component darks by script',
        'cti_wrapper.py on {}'.format(datetime.datetime.now().isoformat(' ')),
        'using basedark file {}.'.format(basedark)]
    copy_dark_keywords(outname, dark_hdr0, history=history, basedark=basedark)
    
    if verbose:
        print 'Weekdark complete:  {}\n'.format(outname)


def populate_darkfiles(raw_files, dark_dir, ref_dir, num_processes, all_weeks_flag=False, verbose=False):
    '''Check science files for uncorrected super-darks; and, if necessary, generate them and
       populate the science file headers.'''
    
    #import pickle
    import glob
    
    superdark_remakes = []
    for file in raw_files:
        with fits.open(file) as f:
            hdr0 = f[0].header
        superdark = hdr0['DARKFILE']
        superdark_resolved = resolve_iraf_file(superdark)
        #if verbose:
        #    print superdark, ' = ', superdark_resolved
        try:
            with fits.open(superdark_resolved) as sd:
                hdr0 = sd[0].header
                if hdr0.get('PCTECORR', default='unknown').strip().upper() == 'COMPLETE':
                    if verbose:
                        print 'Superdark {} is already CTI-corrected.'.format(superdark_resolved)
                    pass
                else:
                    if verbose:
                        print 'Superdark {} is not CTI-corrected.'.format(superdark_resolved)
                    superdark_remakes.extend(file)
        except IOError:
            if verbose:
                print 'Superdark {} not found!  Will make a CTI-corrected version.'.format(superdark_resolved)
            superdark_remakes.extend(file)
    
    if len(superdark_remakes) == 0:
        if verbose:
            print 'Not remaking any superdarks.\n'
        return  # Nothing to remake!
    
    if verbose:
        print 'Making superdarks for:'
        print '   ' + '\n   '.join(raw_files) + '\n'
    
    # Determine component darks used to make superdarks:
    anneal_data = archive_dark_query.get_anneal_boundaries()
    #with open('/Users/lockwood/stis_cte/wrapper/anneals.p', 'rb') as p:
    #    anneal_data = pickle.load(p)  # *** For testing only!!! ***
    #print 'WARNING:  *** Using pickle file anneals.p for testing! ***\n'
    anneals = archive_dark_query.archive_dark_query( \
                      raw_files, anneal_data=anneal_data, print_url=False)  # print_url?
    
    # Get list of EXPNAMEs from files in the dark_dir.
    found_dark_files = glob.glob(os.path.join(dark_dir, '*_flt.fits*'))  # Do something to allow RAW files? ***
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
        print 'ERROR:  These FLT component darks are missing from {}:'.format(dark_dir)
        print ', '.join(missing_darks) + '\n'
        print 'Please download the missing darks via this link:'
        print '(or specify the proper dark_dir [{}])\n'.format(dark_dir)
        print archive_dark_query.darks_url(missing_darks) + '\n'
        raise IOError('Missing component dark FLT files.')
    
    if verbose:
        print 'All required component dark FLT files for annealing periods have been located on disk.\n'
    
    # Determine week_num time boundaries (approximate for now):
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
    
    if verbose >= 2:
        # This should be after filtering by amp...
        print 'Could make weekdarks for:'
        for weekdark_tag in weekdarks:
            print '   {}:  {} - {}  [{}]'.format(
                weekdark_tag, 
                weekdarks[weekdark_tag]['start'], 
                weekdarks[weekdark_tag]['end'], 
                len(weekdarks[weekdark_tag]['darks']) )
        print
    
    weekdark = defaultdict(list)
    for file in raw_files:
        with fits.open(file, 'update') as f:
            hdr0 = f[0].header
            
            dt = datetime.datetime.strptime( \
                hdr0['TDATEOBS'].strip() + ' ' + hdr0['TTIMEOBS'].strip(), '%Y-%m-%d %H:%M:%S')
            weekdark_tag = [wd['weekdark_tag'] for wd in weekdarks.values() if \
                wd['start'] <= dt and wd['end'] > dt and \
                wd['amp'] == hdr0['CCDAMP']][0]
            weekdark[weekdark_tag].append(file)
            
            # Update hdr0 of science file:
            # DARKFILE:
            darkfile = os.path.join(ref_dir, weekdark_tag + '_drk.fits')  # Do something smart with system variables here? ***
            old_darkfile = f[0].header['DARKFILE']
            f[0].header['DARKFILE'] = darkfile
            
            #PCTETAB:
            try:
                old_pctetab = f[0].header['PCTETAB']
            except KeyError:
                old_pctetab = None
            if old_pctetab is None:
                f[0].header.insert('DARKFILE', ('PCTETAB', pctetab))  # Insert after DARKFILE
            else:
                f[0].header['PCTETAB'] = pctetab
            
            f.flush()
            if verbose:
                print 'Updated hdr0 DARKFILE of {}:  {} --> {}'.format(file, old_darkfile, darkfile)
                print 'Updated hdr0 PCTETAB  of {}:  {} --> {}'.format(file, old_pctetab, pctetab)
    
    if verbose:
        print '\nWeekdarks needed:'
        for weekdark_tag in weekdark.keys():
            print '   {} [{}]:'.format(weekdark_tag, len(weekdarks[weekdark_tag]['darks']))
            print '      ' + '\n      '.join(weekdark[weekdark_tag])
        print
    
    if all_weeks_flag:
        weekdark_tags = weekdarks.keys()
        if verbose:
            print 'Generating superdarks for all amps/weeks available.\n'
    else:
        weekdark_tags = weekdark.keys()  # missing 's' -- Need better variable names! ***
    
    # Generate specified basedarks and weekdarks, based on weekdark_tags:
    for weekdark_tag in weekdark_tags:
        # Make basedark:
        amp = weekdarks[weekdark_tag]['amp'].strip().lower()
        basedark = os.path.join(ref_dir, 'basedark_{}{}_drk.fits'.format( \
            weekdarks[weekdark_tag]['amp'].strip().lower(), weekdarks[weekdark_tag]['anneal_num']))
        anneal = [a for a in anneals if a['index'] == weekdarks[weekdark_tag]['anneal_num']][0]
        files = [f['file'] for f in anneal['darks'] if f['CCDAMP'] == weekdarks[weekdark_tag]['amp']]
        generate_basedark(files, basedark, pctetab, num_processes, verbose)  # Or, do we want to compartmentalize the amp selection?
        
        # Make weekdark:
        weekdark_name = os.path.join(ref_dir, weekdark_tag + '_drk.fits')
        files = [f['file'] for f in weekdarks[weekdark_tag]['darks']]  # Already selected for amp
        generate_weekdark(files, weekdark_name, pctetab, basedark, num_processes, verbose)


class Logger(object):
    '''Lumberjack class - Duplicates sys.stdout to a log file and it's okay
       source: http://stackoverflow.com/a/24583265
       Modified to include STDERR.
    '''
    def __init__(self, filename='cti_{}.log'.format(datetime.datetime.now().isoformat('_')), mode="a", buff=0):
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        self.file = open(filename, mode, buff)
        sys.stdout = self
        sys.stderr = self
    
    def __del__(self):
        self.close()
    
    def __enter__(self):
        pass
    
    def __exit__(self, *args):
        pass
    
    def write(self, message):
        self.stdout.write(message)  # Both STDOUT and STDERR get directed to STDOUT!
        self.file.write(message)
    
    def flush(self):
        self.stdout.flush()
        self.stderr.flush()
        self.file.flush()
        os.fsync(self.file.fileno())
    
    def close(self):
        if self.stdout != None:
            sys.stdout = self.stdout
            self.stdout = None
        
        if self.stderr != None:
            sys.stderr = self.stderr
            self.stderr = None
        
        if self.file != None:
            self.file.close()
            self.file = None


if __name__ == '__main__':
    import argparse
    
    # Get information about the user's system:
    num_available_cores = multiprocessing.cpu_count()
    
    # The suggested number of cores is num_available_cores - 2, within [1, default_max_cores]:
    default_max_cores = 15
    default_cores = min([max([1, num_available_cores - 2]), default_max_cores])
    
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
    parser.add_argument('-n', dest='num_processes', action='store', default=default_cores, \
                        metavar='NUM_PROCESSES', type=int, \
                        help='maximum number of parallel processes to run ' + \
                             '(default=' + str(default_cores) + ')' + \
                              "; number of available CPU cores on your system = " + str(num_available_cores))
    parser.add_argument('-p', dest='pctetab', action='store', metavar='PCTETAB', default=None, \
                        help='name of PCTETAB to use in pixel-based correction ' + \
                             '(default=\"[REF_DIR]/test_pcte.fits\")')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False, 
                        help='print more information')
    parser.add_argument('-vv', dest='very_verbose', action='store_true', 
                        help='very verbose')
    # Allow any amp/gain/offst/date through processing (not well-tested):
    parser.add_argument('--allow', dest='allow', action='store_true', default=False, 
                        help=argparse.SUPPRESS)
    parser.add_argument('--all_weeks', dest='all_weeks_flag', action='store_true', default=False, 
                        help=argparse.SUPPRESS)
    args = parser.parse_args()
    
    verbose = args.verbose
    if args.very_verbose:
        verbose = 2
    
    # Determine default dark_dir and ref_dir relative to science_dir:
    if args.dark_dir is None:
        args.dark_dir = os.path.join(args.science_dir, os.path.pardir + os.path.sep + 'darks')
    if args.ref_dir is None:
        args.ref_dir = os.path.join(args.science_dir, os.path.pardir + os.path.sep + 'ref')
    
    # Normalize redundant relative path variables:
    science_dir = os.path.normpath(args.science_dir) + os.path.sep
    dark_dir    = os.path.normpath(args.dark_dir)    + os.path.sep
    ref_dir     = os.path.normpath(args.ref_dir)     + os.path.sep
    
    # Check that directories exist:
    if not os.path.isdir(science_dir):
        raise IOError('science_dir does not exist:  ' + science_dir)
    if not os.path.isdir(dark_dir):
        raise IOError('dark_dir does not exist:  '    + dark_dir)
    if not os.path.isdir(ref_dir):
        raise IOError('ref_dir does not exist:  '     + ref_dir)
    
    # Determine default PCTETAB:
    if args.pctetab is None:
        pctetab = os.path.join(ref_dir, 'test_pcte.fits')
    else:
        pctetab = args.pctetab
    if not os.path.exists(pctetab):
        raise IOError('PCTETAB does not exist:  ' + pctetab)
    
    cti_wrapper(science_dir, dark_dir, ref_dir, pctetab, args.num_processes, 
                args.all_weeks_flag, args.allow, verbose)
    
    #pdb.set_trace()
