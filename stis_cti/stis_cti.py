#! /usr/bin/env python

import os
import sys
import glob
import multiprocessing
import datetime
from astropy.io import fits
from collections import defaultdict
import refstis
from stistools import basic2d, calstis
import archive_dark_query
import StisPixCteCorr
#try:
#    import ipdb as pdb
#except ImportError:
#    import pdb

__author__  = 'Sean Lockwood'
__version__ = '0.3_alpha'


# stis_cti()
# determine_input_science()
# viable_ccd_file()
# bias_correct_science_files()
# run_calstis_on_science()
# resolve_iraf_file()
# superdark_hash()
# check_pctetab_version()
# func()
# func_star()
# perform_cti_correction()
# copy_dark_keywords()
# generate_basedark()
# generate_weekdark()
# populate_darkfiles()
# check_for_old_output_files()
# map_outputs()
# class Logger
# main


class FileError(Exception):
    pass

class VersionError(Exception):
    pass

def stis_cti(science_dir, dark_dir, ref_dir, pctetab, num_processes, 
             all_weeks_flag=False, allow=False, clean=False, verbose=False):
    '''
    Run STIS/CCD pixel-based CTI-correction on data specified in SCIENCE_DIR.
    
    Uncorrected component darks are read from DARK_DIR, and corrected component
    darks are written there too. Corrected super-darks are read from and stored to
    REF_DIR.
    
    From the command line:
        usage: stis_cti.py [-h] [-d DARK_DIR] [-r REF_DIR] [-n NUM_PROCESSES]
                           [-p PCTETAB] [--clean] [-v | -vv]
                           [SCIENCE_DIR]
        
        Run STIS/CCD pixel-based CTI-correction on data specified in SCIENCE_DIR.
        Uncorrected component darks are read from DARK_DIR, and corrected component
        darks are written there too. Corrected super-darks are read from and stored to
        REF_DIR.
        
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
                            (default="[REF_DIR]/[MOST_RECENT]_pcte.fits")
          --clean           remove intermediate and final products from previous runs
                            of this script ('*.txt' files are skipped and clobbered)
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
    print 'Number of parallel processes:   {}'.format(num_processes)
    print 'Start time:                     {}\n'.format(datetime.datetime.now().isoformat(' '))
    
    if verbose:
        print 'cwd         = {}'.format(os.getcwd() + os.path.sep)
        print 'science_dir = {}'.format(science_dir)
        print 'dark_dir    = {}'.format(dark_dir)
        print 'ref_dir     = {}'.format(ref_dir)
        print 'PCTETAB     = {}\n'.format(pctetab)
    log.flush()
    
    # Put dark_dir and PCTETAB's dir in environmental variables:
    if '$' not in ref_dir:
        os.environ['ctirefs'] = os.path.abspath(ref_dir) + os.path.sep
        ref_dir = '$ctirefs' + os.path.sep
    if '$' not in pctetab:
        os.environ['pctetab'] = os.path.abspath(os.path.dirname(pctetab)) + os.path.sep
        # dir$file format needed for StisPixCteCorr()...
        pctetab = 'pctetab$' + os.path.basename(pctetab)  # This won't work with os.path.expandvars()!
    
    # Check PCTETAB's version against min/max allowed by this code:
    if not allow:
        check_pctetab_version(pctetab, verbose)
    
    raw_files = determine_input_science(science_dir, allow, verbose)
    log.flush()
    
    # Rename output files according to this:
    output_mapping = {
        'cte_flt.fits' : 'flc.fits' ,
        'cte_crj.fits' : 'crc.fits' ,
        'cte_sx2.fits' : 's2c.fits' ,
        'cte_x2d.fits' : 'x2c.fits' ,
        'cte_sx1.fits' : 's1c.fits' ,
        'cte_x1d.fits' : 'x1c.fits' ,
        'blt_tra.txt'  : 'trb.txt'  ,
        'cte_tra.txt'  : 'trc.txt'  ,
        'blt.fits'     : '<pass>'   ,
        'cte.fits'     : '<pass>'   }
    
    # Check for original MAST data results (*** Not needed? ***):
    # ...
    
    # Check for results from previous runs:
    rootnames = [os.path.basename(f).split('_',1)[0] for f in raw_files]
    check_for_old_output_files(rootnames, science_dir, output_mapping, clean, verbose)
    # (Note that the 'clean' option doesn't remove '*.txt' files.)
    
    # Check science files for uncorrected super-darks:
    populate_darkfiles(raw_files, dark_dir, ref_dir, pctetab, num_processes, all_weeks_flag, verbose)
    log.flush()
    
    # Bias-correct the science files:
    bias_corrected = bias_correct_science_files(raw_files, verbose)
    log.flush()
    
    # Perform the CTI correction in parallel on the science data:
    cti_corrected = perform_cti_correction(bias_corrected, pctetab, num_processes, verbose)
    log.flush()
    
    # Finish running CalSTIS on the CTI-corrected science data:
    flts = run_calstis_on_science(cti_corrected, verbose)
    log.flush()
    
    # Delete intermediate products and rename final products:
    map_outputs(rootnames, science_dir, output_mapping, verbose)
    log.flush()
    
    print '\nCompletion time:                {}'.format(datetime.datetime.now().isoformat(' '))
    print 'stis_cti.py complete!\n'
    
    log.close()


def determine_input_science(science_dir, allow=False, verbose=False):
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
    
    not_correcting = list(set(raw_files) - set(filtered_raw_files))
    if len(not_correcting) > 0:
        not_correcting.sort()
        print 'WARNING:  Not running the correction on these files:'
        print '   ' + '\n   '.join(not_correcting) + '\n'
    
    # Make sure some files exist:
    if len(filtered_raw_files) == 0:
        raise FileError('No viable STIS CCD files were present in the directory!\n' + \
                        '   ' + science_dir + '\n'                                  + \
                        'Please make sure science data are:\n'                      + \
                        '   -- post-May 01 2009\n'                                  + \
                        '   -- CCDAMP = D\n'                                        + \
                        '   -- CCDGAIN = {1, 4}\n'                                  + \
                        '   -- Not OBSMODE = ACQ*\n'                                + \
                        '   -- Full-frame and not binned')
    
    if verbose:
        print 'Input _raw.fits files being corrected:'
        print '   ' + '\n   '.join(filtered_raw_files) + '\n'
    
    # Check to see if IMPHTTAB is defined for files that will run PHOTCORR:
    missing_imphttab = []
    for file in filtered_raw_files:
        with fits.open(file) as f:
            if (f[0].header.get('PHOTCORR', default='').strip() == 'PERFORM' and \
                f[0].header.get('IMPHTTAB', default=None) is None):
                missing_imphttab.append(file)
    if len(missing_imphttab) > 0:
        raise IOError('These files will fail CalSTIS reduction:\n' + \
                      '   ' + '\n   '.join(missing_imphttab) + \
                      '\nPlease set the IMPHTTAB keyword in files that will perform PHOTCORR.')
    if verbose:
        print 'All files running PHOTCORR have the IMPHTTAB set.\n'
    
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
        raise TypeError('earliest_date_allowed must be a datetime.datetime, not a {}.'.format(type(earliest_date_allowed)))
    
    if amplifiers_allowed is None:
        amplifiers_allowed = ['D']
    
    if gains_allowed is None:
        gains_allowed = [1, 4]
    
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
        hdr0['CCDOFFST'] in offsts_allowed             and \
        not hdr0['SUBARRAY']                           and \
        hdr0['BINAXIS1'] == 1 and hdr0['BINAXIS2'] == 1


def bias_correct_science_files(raw_files, verbose):
    if verbose:
        print 'Bias-correcting science files...\n'
    
    outnames = [f.replace('_raw.fits', '_blt.fits', 1) for f in raw_files]
    trailers = [os.path.abspath(os.path.expandvars(f.replace('_raw.fits', '_blt_tra.txt', 1))) for f in raw_files]
    
    # Check for previous _blt.fits files first:
    for outname in outnames:
        if os.path.exists(outname):
            raise IOError('File {} already exists!'.format(outname))
    
    for raw_file, outname, trailer in zip(raw_files, outnames, trailers):
        if verbose:
            print 'Running basic2d on {} --> {}.'.format(raw_file, outname)
        if os.path.exists(trailer):
            os.remove(trailer)
        cwd = os.getcwd()
        try:
            # Need to change to science directory to find associated epc files.
            os.chdir(os.path.dirname(raw_file))
            status = basic2d.basic2d(os.path.basename(raw_file), output=os.path.basename(outname), 
                dqicorr='perform', atodcorr='omit', blevcorr='perform', biascorr='perform', 
                doppcorr='omit', lorscorr='omit', glincorr='omit', lflgcorr='omit', 
                darkcorr='omit', flatcorr='omit', shadcorr='omit', photcorr='omit', 
                statflag=True, verbose=(verbose >= 2), timestamps=False, trailer=trailer)
        finally:
            os.chdir(cwd)
        
        if verbose or status != 0:
            with open(os.path.expandvars(trailer)) as tra:
                for line in tra.readlines():
                    print '     ' + line.strip()
            print
        
        if status != 0:
            raise RuntimeError('basic2d returned non-zero status on {}:  {}'.format(raw_file, status))
        

        # Remove trailer file?
        
    return outnames


def run_calstis_on_science(files, verbose):
    if verbose:
        print 'Running CalSTIS on science files...\n'
    
    outnames = [f.replace('_cte.fits', '_cte_flt.fits', 1) for f in files]  # Replicating CalSTIS' behavior
    trailers = [os.path.abspath(os.path.expandvars(f.replace('_cte.fits', '_cte_tra.txt', 1))) for f in files]
    #trailers = [os.path.basename(f.replace('_cte.fits', '_cte_tra.txt',  1)) for f in files]  # Also remove path
    
    # Check for previous _blt.fits files first:
    for outname in outnames:
        if os.path.exists(outname):
            raise IOError('File {} already exists!'.format(outname))
    
    for file, outname, trailer in zip(files, outnames, trailers):
        if verbose:
            print 'Running calstis on {} --> {}.'.format(file, outname)
        if os.path.exists(trailer):
            os.remove(trailer)
        # Note that outname is determined by CalSTIS, and not this call:
        cwd = os.getcwd()
        try:
            # Need to change to science directory to find associated wavecals.
            os.chdir(os.path.dirname(file))
            status = calstis.calstis(os.path.basename(file), verbose=(verbose >= 2), 
                timestamps=False, trailer=trailer)
        finally:
            os.chdir(cwd)
        
        if verbose or status != 0:
            with open(os.path.expandvars(trailer)) as tra:
                for line in tra.readlines():
                    print '     ' + line.strip()
            print
        
        if status != 0:
            raise RuntimeError('CalSTIS returned non-zero status on {}:  {}'.format(file, status))
    
    return outnames


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
            raise IOError('Can\'t resolve environmental variable \'{}\'.'.format(dir))
        else:
            dir = dir_resolved
    
    new_file = os.path.normpath(os.path.join(dir, rootname))
    
    return new_file


def check_pctetab_version(pctetab, verbose=False, version_min='0.1', version_max='1.999'):
    '''
    Make sure the version keyword in the PCTETAB is within the acceptable boundaries
    for this version of the stis_cti.py code.
    
    This comparison works for version = "<int>.<int>" (e.g. "1.10" > "1.1").
    '''
    if type(version_min) != str or type(version_max) != str:
        raise TypeError('Versions must be strings.')
    
    # Handle environmental variables in the PCTETAB:
    pctetab = os.path.normpath(resolve_iraf_file(pctetab))
    
    # Read header keyword VERSION and strip anything after '_':
    with fits.open(pctetab) as p:
        pctetab_version_raw = p[0].header.get('PCTE_VER', default='').strip()
        pctetab_version = pctetab_version_raw.split('_')[0]
    
    if pctetab_version == '':
        raise VersionError('PCTE_VER keyword not found in PCTETAB {}.'.format(pctetab))
    
    pctetab_version_major = int(pctetab_version.split('.')[0])
    pctetab_version_minor = int(pctetab_version.split('.')[1])
    
    if (pctetab_version_major   < int(version_min.split('.')[0])  or
        (pctetab_version_major == int(version_min.split('.')[0]) and
         pctetab_version_minor  < int(version_min.split('.')[1]))):
        raise VersionError(('Version mismatch between PCTETAB {} and stis_cti.py code.\n' + 
            'Please download a more recent version of this reference file!\n').format(pctetab))
       
    if (pctetab_version_major   > int(version_max.split('.')[0])  or
        (pctetab_version_major == int(version_max.split('.')[0]) and
         pctetab_version_minor  > int(version_max.split('.')[1]))):
        raise VersionError(('Version mismatch between PCTETAB {} and stis_cti.py code.\n' +
            'Please upgrade this code!\n').format(pctetab))\
    
    if verbose:
        print 'Version of PCTETAB:  {}\n'.format(pctetab_version_raw)
    
    return True


def superdark_hash(sim_nit=None, shft_nit=None, rn_clip=None, nsemodel=None, subthrsh=None, pcte_ver=None,
                            pctetab=None,
                            superdark=None, files=None):
    from numpy import where
    
    if pctetab is not None:
        with fits.open(resolve_iraf_file(pctetab)) as f:
            sim_nit  = f[0].header['SIM_NIT' ]
            shft_nit = f[0].header['SHFT_NIT']
            rn_clip  = f[0].header['RN_CLIP' ]
            nsemodel = f[0].header['NSEMODEL']
            subthrsh = f[0].header['SUBTHRSH']
            pcte_ver = f[0].header['PCTE_VER'].strip()
        if files is None:
            raise IOError('Must specify file list with pctetab.')
    elif superdark is not None:
        with fits.open(os.path.expandvars(superdark)) as f:
            sim_nit  = f[0].header.get('PCTESMIT', default='undefined')
            shft_nit = f[0].header.get('PCTESHFT', default='undefined')
            rn_clip  = f[0].header.get('PCTERNCL', default='undefined')
            nsemodel = f[0].header.get('PCTENSMD', default='undefined')
            subthrsh = f[0].header.get('PCTETRSH', default='undefined')
            pcte_ver = f[0].header.get('PCTE_VER', default='0.0').strip()
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
           subthrsh == None or \
           pcte_ver == None:
            raise IOError('Please specify either pctetab or the proper set of parameters!')
        if files is None:
            raise IOError('Must specify files with pctetab.')
    
    # Ignore PCTETAB version text after underscore:
    pcte_ver = pcte_ver.split('_',1)[0]
    
    # Turn list of filenames into a sorted list of exposures:
    exposures = [os.path.basename(file).rsplit('_',1)[0].upper() for file in files]
    exposures.sort()
    exposures = ','.join(exposures)
    
    hash_str = '{};{};{};{};{};{};{}'.format( \
        exposures, \
        sim_nit, shft_nit, rn_clip, nsemodel, subthrsh, pcte_ver)
    
    #return hash_str
    return hash(hash_str)


# These functions are needed here for calling from a multiprocessing pool.
def func(in_f, out_f):
    StisPixCteCorr.CteCorr(in_f, outFits=out_f)

def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return func(*a_b)


def perform_cti_correction(files, pctetab, num_cpu=1, verbose=False):
    # The call to StisPixCteCorr should be done in parallel!
    
    perform_files = []
    outnames = []
    perform_outnames = []
    for file in files:
        outname = file.replace('_flt.fits', '_cte.fits', 1).replace('_blt.fits', '_cte.fits', 1)
        outnames.append(outname)
        if os.path.exists(outname) and (superdark_hash(pctetab=pctetab, files=[]) == superdark_hash(superdark=outname, files=[])):
            if verbose:
                print 'Skipping regeneration of CTI-corrected file:  {}'.format(outname)
        else:
            with fits.open(file, 'update') as f:
                #Update PCTETAB header keyword:
                old_pctetab = f[0].header.get('PCTETAB', default=None)
                f[0].header.set('PCTETAB', pctetab, 'Pixel-based CTI param table', after='DARKFILE')
                
                # Turn off the empirical CTI correction flag if it is set to PERFORM:
                old_ctecorr = f[0].header.get('CTECORR', default='unknown').strip()
                if old_ctecorr == 'COMPLETE':
                    raise ValueError('Empirical CTI correction correction flag, CTECORR, is set in file {}.'.format(file))
                elif old_ctecorr == 'PERFORM':
                    f[0].header['CTECORR'] = 'OMIT'
                
                f.flush()
                if verbose:
                    print 'Updated hdr0 PCTETAB  of {}:  {} --> {}'.format(file, old_pctetab, pctetab)
                    if old_ctecorr != 'unknown':
                        print 'Updated hdr0 CTECORR  of {}:  {} --> {}'.format(
                            file, old_ctecorr, f[0].header['CTECORR'])
                
                # Run the pixel-based correction on these files:
                perform_files.append(file)
                perform_outnames.append(outname)
    
    # Run the CTI-correction:
    p = multiprocessing.Pool(processes = num_cpu)
    # Iterator of input/output files:
    file_args = ((x, y) for x, y in zip(perform_files, perform_outnames))
    p.map_async(func_star, file_args)
    p.close()
    p.join()
    
    ## Single-threaded version:
    #for perform_file, outname in zip(perform_files, outnames):
    #    StisPixCteCorr.CteCorr(perform_file, outFits=outname)
    
    # *** Do the corrected files need to be fed through DQICORR again to fix flags? ***
    
    if len(outnames) == 0:
        outnames = None
    
    return outnames


def copy_dark_keywords(superdark, dark_hdr0, pctetab, history=None, basedark=None):
    # Copy these keywords from the last component dark to the new superdark:
    keywords = ['PCTECORR', 'PCTETAB', 'PCTEFRAC', 'PCTERNCL', 'PCTERNCL', 'PCTENSMD', 
                'PCTETRSH', 'PCTESMIT', 'PCTESHFT', 'CTE_NAME', 'CTE_VER', 'PCTE_VER']
    # Note:  PCTEFRAC is time-dependent; PCTECORR = COMPLETE
    keywords.reverse()
    
    with fits.open(os.path.expandvars(superdark), 'update') as s:
        for keyword in keywords:
            value = dark_hdr0.get(keyword, default='unknown')
            try:
                comment = dark_hdr0.comments[keyword]
            except KeyError:
                comment = None
            
            if keyword in s[0].header:
                s[0].header[keyword] = value
            else:
                s[0].header.set(keyword, value, comment, after='DRK_VS_T')
            s.flush()
        
        if basedark is not None:
            if 'BASEDARK' in s[0].header:
                s[0].header['BASEDARK'] = basedark
            else:
                s[0].header.set('BASEDARK', basedark, 'Used to make weekdark', after=keywords[0])
        
        # Add HISTORY to superdark hdr0 here:
        if history is not None:
            s[0].header['HISTORY'] = ' '
            for h in history:
                s[0].header['HISTORY'] = h


def generate_basedark(files, outname, pctetab, num_cpu, verbose=False):
    # Don't make a basedark if it already exists:
    if os.path.exists(os.path.expandvars(outname)):
        if superdark_hash(pctetab=pctetab, files=files) == superdark_hash(superdark=outname):
            if verbose:
                print 'Skipping regeneration of basedark:  {}'.format(outname)
            return
    
    if verbose >= 2:
        print 'Working on basedark {}:'.format(outname)
        print '   ' + '\n   '.join(files) + '\n'
    
    # Correct component darks, if necessary:
    corrected_files = perform_cti_correction(files, pctetab, num_cpu, verbose)
    
    # Make a basedark from the corrected darks:
    refstis.basedark.make_basedark(corrected_files, refdark_name=os.path.normpath(os.path.expandvars(outname)))
    # Print results of "[REF_DIR]/[basedark_name - '.fits']_joined_bd_calstis_log.txt"? ***
    
    if verbose:
        calstis_log = outname.replace('.fits','_joined_bd_calstis_log.txt', 1)
        print 'Calstis log file for CRJ processing [{}]:'.format(calstis_log)
        with open(os.path.expandvars(calstis_log), 'r') as cs_log:
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
        'stis_cti.py on {}.'.format(datetime.datetime.now().isoformat(' ')) ]
    copy_dark_keywords(os.path.expandvars(outname), dark_hdr0, pctetab, history=history)
        
    if verbose:
        print 'Basedark complete:  {}\n'.format(outname)


def generate_weekdark(files, outname, pctetab, basedark, num_cpu, verbose=False):
    # Don't make a weekdark if it already exists:
    if os.path.exists(os.path.expandvars(outname)):
        if superdark_hash(pctetab=pctetab, files=files) == superdark_hash(superdark=outname):
            if verbose:
                print 'Skipping regeneration of weekdark:  {}'.format(outname)
            return
    
    if verbose >= 2:
        print 'Working on weekdark {}:'.format(outname)
        print '   ' + '\n   '.join(files) + '\n'
    
    # Correct component darks, if necessary:
    corrected_files = perform_cti_correction(files, pctetab, num_cpu, verbose)
    
    # Make a weekdark from the corrected darks:
    refstis.weekdark.make_weekdark(corrected_files, os.path.normpath(os.path.expandvars(outname)), 
        os.path.abspath(os.path.expandvars(basedark)))
    
    if verbose:
        calstis_log = outname.replace('.fits','_joined_bd_calstis_log.txt', 1)
        print 'Calstis log file for CRJ processing [{}]:'.format(calstis_log)
        with open(os.path.expandvars(calstis_log), 'r') as cs_log:
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
        'stis_cti.py on {}'.format(datetime.datetime.now().isoformat(' ')),
        'using basedark file {}.'.format(basedark)]
    copy_dark_keywords(os.path.expandvars(outname), dark_hdr0, pctetab, history=history, basedark=basedark)
    
    if verbose:
        print 'Weekdark complete:  {}\n'.format(outname)


def populate_darkfiles(raw_files, dark_dir, ref_dir, pctetab, num_processes, all_weeks_flag=False, verbose=False):
    '''Check science files for uncorrected super-darks; and, if necessary, generate them and
       populate the science file headers.'''
    
    #import pickle
    
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
        max_file_length = str(max(map(lambda x: len(x), found_dark_files) + [1]))
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
        print 'Please download the missing darks (calibrated FLTs) via this link:'
        print '(or specify the proper dark_dir [{}])\n'.format(dark_dir)
        print archive_dark_query.darks_url(missing_darks) + '\n'
        raise FileError('Missing component dark FLT files.')
    
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
            matched_weekdark_tags = [wd['weekdark_tag'] for wd in weekdarks.values() if \
                wd['start'] <= dt and wd['end'] > dt and \
                wd['amp'] == hdr0['CCDAMP']]
            if len(matched_weekdark_tags) == 0:
                raise IOError(('No matching component darks found to make needed weekdark for amp {}\n' + \
                    'within annealing period boundaries corresponding to science set {}.')\
                    .format(hdr0['CCDAMP'], os.path.basename(file)))
            weekdark_tag = matched_weekdark_tags[0]
            weekdark[weekdark_tag].append(file)
            
            # Update hdr0 of science file:
            # DARKFILE:
            darkfile = os.path.join(ref_dir, weekdark_tag + '_drk.fits')  # Do something smart with system variables here? ***
            old_darkfile = f[0].header['DARKFILE']
            f[0].header['DARKFILE'] = darkfile
            
            #PCTETAB:
            #old_pctetab = f[0].header.get('PCTETAB', default=None)
            #f[0].header.set('PCTETAB', pctetab, after='DARKFILE')
            
            f.flush()
            if verbose:
                print 'Updated hdr0 DARKFILE of {}:  {} --> {}'.format(file, old_darkfile, darkfile)
                #print 'Updated hdr0 PCTETAB  of {}:  {} --> {}'.format(file, old_pctetab, pctetab)
    
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
        weekdark_name = os.path.abspath(os.path.expandvars(os.path.join(ref_dir, weekdark_tag + '_drk.fits')))
        files = [f['file'] for f in weekdarks[weekdark_tag]['darks']]  # Already selected for amp
        generate_weekdark(files, weekdark_name, pctetab, basedark, num_processes, verbose)


def check_for_old_output_files(rootnames, science_dir, output_mapping, clean=False, verbose=False):
    # Combine science path with rootnames:
    path_rootnames = [os.path.join(science_dir, r) for r in rootnames]
    
    # Combine rootnames with temporary and final output names:
    old_files = []
    for ext in output_mapping.keys():
        old_files.extend([f + '_' + ext for f in path_rootnames
                          if '.txt' not in ext])
        old_files.extend([f + '_' + output_mapping[ext] for f in path_rootnames
                          if '<' not in output_mapping[ext] and '.txt' not in output_mapping[ext]])
    
    # Include gzipped versions of filenames:
    tmp = [f + '.gz' for f in old_files]
    tmp.extend([f + '.GZ' for f in old_files])
    
    # Check for the existence of these files:
    error_files = []
    old_files.extend(tmp)
    for old_file in old_files:
        if os.path.exists(old_file) or os.path.islink(old_file):
            if clean:
                if verbose:
                    print 'Removing old file:  {}'.format(old_file)
                os.remove(old_file)
            else:
                error_files.append(old_file)
    if len(error_files) > 0:
        error_files.sort()
        raise IOError('Files exist from previous run of stis_cti.py:  \n{}'.format( \
                      '   ' + '\n   '.join(error_files) + \
                      '\nYou might consider running with the \'--clean\' option specified.'))
    if verbose:
        print 'The science_dir is clear of files from previous runs.\n'
    
    return True


def map_outputs(rootnames, science_dir, output_mapping, verbose=False):
    for rootname in rootnames:
        if verbose:
            print 'Renaming files with rootname {}:'.format(rootname)
        files = glob.glob(os.path.join(science_dir, rootname + '*'))
        for file in files:
            cwd = os.getcwd()
            try:
                os.chdir(os.path.dirname(file))
                file = os.path.basename(file)
                
                # Determine file extension:
                try:
                    ext = file.split('_',1)[1]
                except IndexError:
                    ext = '<undefined>'  # No underscore in name
                
                # Handle gzipped filenames:
                if ext[-3:] in ['.gz', '.GZ']:
                    gzip = ext[-3:]
                    ext = ext[-3:]
                else:
                    gzip = ''
                
                # Map action to file:
                if ext not in output_mapping.keys():
                    pass
                elif output_mapping[ext] == '<pass>':
                    pass
                elif output_mapping[ext] == '<remove>':
                    if verbose:
                        print 'Deleting:  {}'.format(file)
                    os.remove(file)
                else:
                    new_file = rootname + '_' + output_mapping[ext] + gzip
                    if verbose:
                        print 'Renaming:  {}\t-->\t{}'.format(file, new_file)
                    os.rename(file, new_file)
                    # *** What about the filename keyword? ***
            finally:
                os.chdir(cwd)


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
                             '(default=\"[REF_DIR]/[MOST_RECENT]_pcte.fits\")')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False, 
                        help='remove intermediate and final products from previous runs of ' + \
                             'this script (\'*.txt\' files are skipped and clobbered)')
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
        raise FileError('science_dir does not exist:  ' + science_dir)
    if not os.path.isdir(dark_dir):
        raise FileError('dark_dir does not exist:  '    + dark_dir)
    if not os.path.isdir(ref_dir):
        raise FileError('ref_dir does not exist:  '     + ref_dir)
    
    # Determine default PCTETAB (default = last alphabetical _pcte.fits from ref_dir):
    if args.pctetab is None:
        pctetabs = glob.glob(os.path.join(ref_dir, '*_pcte.fits*'))
        pctetabs.sort()
        if len(pctetabs) > 0:
            pctetab = pctetabs[-1]
        else:
            raise FileError('No PCTETAB file found in {}.'.format(ref_dir))
    else:
        pctetab = args.pctetab
    if not os.path.exists(pctetab):
        raise FileError('PCTETAB does not exist:  ' + pctetab)
    
    stis_cti(science_dir, dark_dir, ref_dir, pctetab, args.num_processes, 
             args.all_weeks_flag, args.allow, args.clean, verbose)
    
