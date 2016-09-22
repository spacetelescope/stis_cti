import os
import sys
import glob
import multiprocessing
import datetime
from astropy.io import fits
from collections import defaultdict
from copy import deepcopy
import refstis
from stistools import basic2d, calstis
import archive_dark_query
import StisPixCteCorr
from crds.bestrefs import BestrefsScript

__author__  = 'Sean Lockwood'
__version__ = '1.1'

crds_server_url = 'https://hst-crds.stsci.edu'

class FileError(Exception):
    pass

class VersionError(Exception):
    pass

def stis_cti(science_dir, dark_dir, ref_dir, num_processes, pctetab=None, 
             all_weeks_flag=False, allow=False, clean=False, clean_all=False, 
             crds_update=False, ignore_missing=False, verbose=1):
    '''
    Runs the HST/STIS/CCD pixel-based CTI-correction on science data and 
    component darks, generating and applying a CTI-corrected super-dark in 
    the process.
    
    Documentation is available at http://pythonhosted.org/stis_cti/
    
    :param science_dir:
        Directory containing uncalibrated science data to be corrected.
    :param dark_dir: str
        Directory containing calibrated component darks to be corrected and
        location where CTI-corrected component darks are placed.
    :param ref_dir: str
        Directory where CTI-corrected super-darks are placed.  These will be
        used again unless they are deleted or clean_all=True.
    :param num_processes:
        Max number of parallel processes to use when running the CTI-correction 
        algorithm.
    :param pctetab:
        The path + name of the PCTETAB reference file to use in the CTI-correction.  
        If not specified, one is selected from (1) the ref_dir, or (2) from the 
        package data directory.  The last file (alphabetically) is chosen.
    :param all_weeks_flag:
        UNTESTED.  Generates weekdarks for all weeks within each annealing period 
        to be processed.
    :param allow:
        UNTESTED.  Use more lenient filtering when determining which files should 
        be allowed to be corrected.
    :param clean:
        Remove intermediate and final products in the science_dir from previous 
        runs of this script.
    :param clean_all:
        'clean' + remove CTI-corrected super-darks and component darks before 
        reprocessing.
    :param crds_update:
        Runs crds.bestrefs script to update science file headers and download 
        pipeline reference files.
    :param ignore_missing:
        Ignore missing dark files.  This is useful for when annealing month contains
        amp=A RAW files, but no associated FLT files.
    :param verbose:
        Verbosity of text printed to the screen and saved in the log file.
    
    :type science_dir: str
    :type dark_dir: str
    :type ref_dir: str
    :type num_processes: int
    :type pctetab: str, optional 
    :type all_weeks_flag: bool, optional
    :type allow: bool, optional
    :type clean: bool, optional
    :type clean_all: bool, optional
    :type crds_update: bool, optional
    :type ignore_missing: bool, optional, default=False
    :type verbose: int {0,1,2}, optional, default=1
    
    .. note::
      Unless crds_update is True, the $oref shell variable must be set to the directory of 
      STIS standard pipeline reference files.
      
    .. note::
      Note that the 'all_week_flag' and 'allow' options have not been tested!
    '''
    # Open a log file and copy {STDOUT, STDERR} to it:
    log = Logger(os.path.join(science_dir, 'cti_{}.log'.format(datetime.datetime.now().isoformat('_'))))
    
    try:
        import pkg_resources
    except ImportError:
        pass
    try:
        from platform import node, platform
        sys_info = '{}; {}'.format(node(), platform())
    except ImportError:
        sys_info = 'Indeterminate'
    
    start_time = datetime.datetime.now()
    
    # Print system information:
    print 'Running CTI-correction script:  {} v{}'.format(os.path.basename(__file__), __version__)
    print 'System:                         {}'.format(sys_info)
    print 'Number of parallel processes:   {}'.format(num_processes)
    if clean:
        print '--clean'
    if clean_all:
        print '--clean_all'
    if crds_update:
        print '--crds_update'
    if ignore_missing:
        print '--ignore_missing'
    if verbose:
        print 'verbose mode:                   {}'.format(verbose)

    print 'Start time:                     {}\n'.format(start_time.isoformat(' '))
    
    # Check that directories exist:
    if not os.path.isdir(science_dir):
        raise FileError('science_dir does not exist:  ' + science_dir)
    if not os.path.isdir(dark_dir):
        raise FileError('dark_dir does not exist:  '    + dark_dir)
    if not os.path.isdir(ref_dir):
        raise FileError('ref_dir does not exist:  '     + ref_dir)
    
    if clean_all:
        clean = True
    
    # Check PCTETAB:
    if pctetab is None:
        pctetabs = glob.glob(os.path.join(ref_dir, '*_pcte.fits*'))
        pctetabs.sort()
        if len(pctetabs) > 0:
            pctetab = pctetabs[-1]
        else:
            # Couldn't find PCTETAB in ref/ directory, so use the PCTETAB included with the stis_cti package:
            try:
                default_pcte_str = pkg_resources.resource_filename(stis_cti.__name__, 'data/*_pcte.fits')
            except NameError:
                default_pcte_str = os.path.join(os.path.dirname(stis_cti.__file__), 'data', '*_pcte.fits')
            package_pctes = glob.glob(default_pcte_str)
            if len(package_pctes) == 0:
                raise FileError('Couldn\'t find a PCTETAB in {}\nor package default in {}'.format(
                    ref_dir, os.path.dirname(default_pcte_str)))
            package_pctes.sort()
            pctetab = package_pctes[-1]
            print 'WARNING:  Using package-default PCTETAB:\n{}\n'.format(pctetab)
    if not os.path.exists(pctetab):
        raise FileError('PCTETAB does not exist:  ' + pctetab)
    
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
    
    # Setup environmental variables for CRDS code to update headers:
    if crds_update:
        setup_crds(ref_dir, verbose)
    
    # Test that $oref is properly defined:
    oref = os.environ.get('oref', failobj='Undefined')
    if oref is 'Undefined' or not os.access(oref, os.R_OK):
        raise OSError('Cannot read $oref directory!\n    {}\n'.format(oref) + \
            '    Please set $oref environmental variable appropriately or run with --crds_update.')
    if verbose:
        print '$oref      = {}\n'.format(oref)
    
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
        'blt.fits'     : '<remove>' ,
        'cte.fits'     : '<pass>'   }
    
    # Check for results from previous runs:
    rootnames = [os.path.basename(f).split('_',1)[0] for f in raw_files]
    check_for_old_output_files(rootnames, science_dir, output_mapping, clean, verbose)
    # (Note that the 'clean' option doesn't remove '*.txt' files.)
    
    # Run crds.BestrefsScript on science files:
    if crds_update:
        # Run correction on matching _wav files too, if they exist:
        raw_path = os.path.dirname(raw_files[0])
        wav_files = []
        for file in raw_files:
            with fits.open(file) as f:
                wav_files.append(f[0].header.get('WAVECAL', default='N/A').strip())
        wav_files = [file for file in wav_files if 'N/A' not in file]
        wav_files = [os.path.join(raw_path, file) for file in wav_files]
        wav_files = [file for file in wav_files if os.path.exists(file)]
        
        crds_files = deepcopy(raw_files)
        crds_files.extend(wav_files)
        
        if verbose:
            print 'Running crds.BestrefsScript on science and wav files...'
        errors = BestrefsScript('BestrefsScript --update-bestrefs -s 1 -f ' + ' '.join(crds_files))()
        if int(errors) > 0:
            raise Exception('CRDS BestrefsScript:  Call returned errors!')
    
    # Check science files for uncorrected super-darks:
    populate_darkfiles(raw_files, dark_dir, ref_dir, pctetab, num_processes, all_weeks_flag, 
        clean_all, crds_update, ignore_missing, verbose)
    log.flush()
    
    # Bias-correct the science files:
    bias_corrected = bias_correct_science_files(raw_files, verbose)
    log.flush()
    
    # Perform the CTI correction in parallel on the science data:
    cti_corrected = perform_cti_correction(bias_corrected, pctetab, num_processes, False, verbose)
    log.flush()
    
    # Finish running CalSTIS on the CTI-corrected science data:
    flts = run_calstis_on_science(cti_corrected, verbose)
    log.flush()
    
    # Delete intermediate products and rename final products:
    map_outputs(rootnames, science_dir, output_mapping, verbose)
    log.flush()
    
    end_time = datetime.datetime.now()
    print '\nCompletion time:                {}'.format(end_time.isoformat(' '))
    print 'Run time:                       {}'.format(end_time - start_time)
    print 'stis_cti.py complete!\n'
    log.close()


def setup_crds(ref_dir, verbose=False):
    '''
    Setup $CRDS_PATH and $oref environmental variables to house CRDS reference files.
    
    Sets the $CRDS_PATH environmental variable according to:
        (1) If undefined, Central Store location (if possible)
        (2) If not on Central Store (local), read/writable directory defined by system's $CRDS_PATH
        (3) If local and system's $CRDS_PATH directory not read/writable, nested within 'ref'
    
    Unless $CRDS_PATH is on Central Store, $oref is redefined to be nested within new $CRDS_PATH.
    
    $oref value must end in the path separator (e.g. '/'), so this is appended if necessary.
    '''
    
    if os.environ.get('CRDS_SERVER_URL') is None:
        os.environ['CRDS_SERVER_URL'] = crds_server_url
    
    if verbose:
        print 'Setting up CRDS environmental variables...'
        print '$CRDS_SERVER_URL = {}'.format(os.environ['CRDS_SERVER_URL'])
    
    ref_dir = resolve_iraf_file(ref_dir)
    
    if os.environ.get('CRDS_PATH') is None and \
       os.path.exists('/grp/crds/cache') and os.path.exists('/grp/crds/cache/references/hst'):
            # On-site CRDS cache on Central Store:
            # (This location is special because we can assume it does not need updating.)
            os.environ['CRDS_PATH'] = '/grp/crds/cache'
            os.environ['oref'] = os.path.join(os.environ['CRDS_PATH'], 'references', 'hst') + \
                                 os.path.sep
    elif os.path.exists('/grp/crds/cache') and \
         os.path.samefile(os.environ.get('CRDS_PATH'), '/grp/crds/cache'):
        # CRDS_PATH is already defined as Central Store:
        os.environ['oref'] = os.path.join(os.environ['CRDS_PATH'], 'references', 'hst') + \
                             os.path.sep
    else:
        if os.environ.get('CRDS_PATH') is None:
            # $CRDS_PATH is undefined and Central Store is not available:
            os.environ['CRDS_PATH'] = os.path.abspath(ref_dir)
            os.environ['oref'] = os.path.join(os.environ['CRDS_PATH'], 'references', 'hst', 'stis') + \
                                 os.path.sep
        elif not os.access(os.environ.get('CRDS_PATH'), os.R_OK | os.W_OK):
            # $CRDS_PATH was already defined by the user, but is NOT read/writable:
            print 'WARNING:  Local $CRDS_PATH is not read/writable.'
            print "         Using 'ref' directory for local CRDS cache."
            os.environ['CRDS_PATH'] = os.path.abspath(ref_dir)
            os.environ['oref'] = os.path.join(os.environ['CRDS_PATH'], 'references', 'hst', 'stis') + \
                                 os.path.sep
        elif os.access(os.environ.get('CRDS_PATH'), os.R_OK | os.W_OK):
            # $CRDS_PATH was already defined by the user and is read/writable:
            
            # Check that $oref is defined within $CRDS_PATH properly for either a flat or
            # instrument-specific cache:
            if os.path.abspath(os.environ.get('oref')) not in \
                [os.path.abspath(os.path.join(os.environ.get('CRDS_PATH'), \
                                              'references', 'hst')), 
                 os.path.abspath(os.path.join(os.environ.get('CRDS_PATH'), \
                                              'references', 'hst', 'stis'))]:
                raise OSError('$oref is not properly defined within $CRDS_PATH.\n' + 
                              '    Consider fixing $oref or running without --crds_update option.')
        else:
            raise Exception('Unexpected condition.')
        
        if not os.access(os.path.abspath(os.environ.get('oref')), os.F_OK):
            try:
                os.makedirs(os.path.abspath(os.environ.get('oref')))
                if verbose:
                    print 'Created $oref directory.'
            except OSError as err:
                print 'Could not create directory structure:'
                print '    {}'.format(os.path.abspath(os.environ.get('oref')))
                raise err
        
        # Check that local $CRDS_PATH and $oref are read/writable:
        if not os.access(os.environ.get('CRDS_PATH'), os.R_OK | os.W_OK) or \
           not os.access(os.environ.get('oref'), os.R_OK | os.W_OK):
            raise OSError(('Local $CRDS_PATH and/or $oref are not read/writable!\n' + 
                           '    $CRDS_PATH = {}\n' + 
                           '    $oref      = {}\n' +
                           '    Either fix these or run without --crds_update option.').format( \
                           os.environ.get('CRDS_PATH'), os.environ.get('oref')))
    
    # Check that $oref ends in the path separator:
    if os.environ['oref'][-1] != os.path.sep:
        os.environ['oref'] += os.path.sep
    
    if not os.access(os.environ['oref'], os.R_OK):
        raise OSError('Cannot read $oref!\n    {}'.format(os.environ['oref']))
    
    if verbose:
        print '$CRDS_PATH = {}'.format(os.environ.get('CRDS_PATH'))
    
    return True


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
    '''
    :param earliest_date_allowed: Beginning datetime allowed.
    :param amplifiers_allowed: CCDAMP values allowed.
    :param gains_allowed:  CCDGAIN values allowed.
    :param offsts_allowed: CCDOFFST values allowed
    
    :type earliest_date_allowed: datetime.datetime object (default=2009-May-01 00:00:00 UTC)
    :type amplifiers_allowed: list of strings (default=['D'])
    :type gains_allowed: list of ints (default=[1,4])
    :type offsts_allowed: list of ints (default=[3])
    
    :returns:  bool -- Run stis_cti.stis_cti() on the file?
    '''
    
    # Set defaults:
    if earliest_date_allowed is None:
        earliest_date_allowed = datetime.datetime(2009, 5, 1, 0, 0, 0)
    if type(earliest_date_allowed) is not datetime.datetime:
        raise TypeError('earliest_date_allowed must be a datetime.datetime, not a {}.'.format(
            type(earliest_date_allowed)))
    
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


def bias_correct_science_files(raw_files, verbose=False):
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
    
    return outnames


def run_calstis_on_science(files, verbose):
    if verbose:
        print 'Running CalSTIS on science files...\n'
    
    # Note that outname is determined by CalSTIS, and not this function.
    outnames = [f.replace('_cte.fits', '_cte_flt.fits', 1) for f in files]  # Replicating CalSTIS' behavior
    trailers = [os.path.abspath(os.path.expandvars(f.replace('_cte.fits', '_cte_tra.txt', 1))) for f in files]
    
    # Check for previous output files:
    for outname in outnames:
        if os.path.exists(outname):
            raise IOError('File {} already exists!'.format(outname))
    
    for file, outname, trailer in zip(files, outnames, trailers):
        if verbose:
            print 'Running calstis on {} --> {}.'.format(file, outname)
        if os.path.exists(trailer):
            os.remove(trailer)
        
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
    '''Resolves pathvar$filename into usable format.'''
    dir = ''
    rootname = file
    
    if '$' in file and file[0] == '$':
        dir, rootname = file[1:].split(os.path.sep, 1)
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
    
    return hash(hash_str)


# These functions are needed here for calling from a multiprocessing pool.
def func(in_f, out_f):
    StisPixCteCorr.CteCorr(in_f, outFits=out_f)

def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return func(*a_b)


def perform_cti_correction(files, pctetab, num_cpu=1, clean_all=False, verbose=False):
    '''Run StisPixCteCorr when needed.'''
    
    perform_files = []
    outnames = []
    perform_outnames = []
    for file in files:
        outname = file.replace('_flt.fits', '_cte.fits', 1).replace('_blt.fits', '_cte.fits', 1)
        outnames.append(outname)
        
        if clean_all and os.path.exists(outname):
            if verbose:
                print 'Deleting file:  {}'.format(outname)
            os.remove(outname)
        
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
    
    # Single-threaded version:
    #for perform_file, outname in zip(perform_files, outnames):
    #    StisPixCteCorr.CteCorr(perform_file, outFits=outname)
    
    if len(outnames) == 0:
        outnames = None
    
    return outnames


def copy_dark_keywords(superdark, dark_hdr0, pctetab, history=None, basedark=None):
    '''Copy header keywords from a used component dark to the new super-dark.'''
    
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


def generate_basedark(files, outname, pctetab, num_cpu, clean_all=False, verbose=False):
    '''Generate a basedark for an annealing period, if it doesn't already exist.'''
    
    if os.path.exists(os.path.expandvars(outname)):
        if clean_all:
            if verbose:
                print 'Deleting basedark:  {}'.format(os.path.expandvars(outname))
            os.remove(os.path.expandvars(outname))
        elif superdark_hash(pctetab=pctetab, files=files) == superdark_hash(superdark=outname):
            # Don't make a basedark if it already exists:
            if verbose:
                print 'Skipping regeneration of basedark:  {}'.format(outname)
            return
    
    if verbose >= 2:
        print 'Working on basedark {}:'.format(outname)
        print '   ' + '\n   '.join(files) + '\n'
    
    # Correct component darks, if necessary:
    corrected_files = perform_cti_correction(files, pctetab, num_cpu, clean_all, verbose)
    
    # Make a basedark from the corrected darks:
    refstis.basedark.make_basedark(corrected_files, refdark_name=os.path.normpath(os.path.expandvars(outname)))
    
    if verbose:
        calstis_log = outname.replace('.fits','_joined_bd_calstis_log.txt', 1)
        print 'Calstis log file for CRJ processing [{}]:'.format(calstis_log)
        with open(os.path.expandvars(calstis_log), 'r') as cs_log:
            for line in cs_log.readlines():
                print '     ' + line.strip()
        print
    
    # Copy the last file's ext=0 header into a variable to use in populating the basedark header:
    with fits.open(corrected_files[-1]) as file:
        dark_hdr0 = file[0].header
    
    # Update keywords in new basedark from component dark:
    history = [
        'Basedark created from CTI-corrected component darks by script',
        'stis_cti.py on {}.'.format(datetime.datetime.now().isoformat(' ')) ]
    copy_dark_keywords(os.path.expandvars(outname), dark_hdr0, pctetab, history=history)
        
    if verbose:
        print 'Basedark complete:  {}\n'.format(outname)


def generate_weekdark(files, outname, pctetab, basedark, num_cpu, clean_all=False, verbose=False):
    '''Generate a weekdark for part of an annealing period, if it doesn't already exist.'''
    
    if os.path.exists(os.path.expandvars(outname)):
        if clean_all:
            if verbose:
                print 'Deleting weekdark:  {}'.format(os.path.expandvars(outname))
            os.remove(os.path.expandvars(outname))
        elif superdark_hash(pctetab=pctetab, files=files) == superdark_hash(superdark=outname):
            # Don't make a weekdark if it already exists:
            if verbose:
                print 'Skipping regeneration of weekdark:  {}'.format(outname)
            return
    
    if verbose >= 2:
        print 'Working on weekdark {}:'.format(outname)
        print '   ' + '\n   '.join(files) + '\n'
    
    # Correct component darks, if necessary:
    corrected_files = perform_cti_correction(files, pctetab, num_cpu, False, verbose)
      # Don't delete component darks here, as this has already been done in generate_basedark().
    
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
    
    # Copy the last file's ext=0 header into a variable to use in populating the basedark header:
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


def populate_darkfiles(raw_files, dark_dir, ref_dir, pctetab, num_processes, all_weeks_flag=False, 
    clean_all=False, crds_update=False, ignore_missing=False, verbose=False):
    '''Check science files for uncorrected super-darks; and, if necessary, generate them and
       populate the science file headers.'''
    
    superdark_remakes = []
    for file in raw_files:
        with fits.open(file) as f:
            hdr0 = f[0].header
        superdark = hdr0['DARKFILE']
        superdark_resolved = resolve_iraf_file(superdark)
        try:
            with fits.open(superdark_resolved) as sd:
                hdr0 = sd[0].header
                if hdr0.get('PCTECORR', default='unknown').strip().upper() == 'COMPLETE':
                    if not clean_all:
                        if verbose:
                            print 'Superdark {} is already CTI-corrected.'.format(superdark_resolved)
                        pass
                    else:
                        if verbose:
                            print 'Remaking CTI-corrected superdark {}'.format(superdark_resolved)
                        if os.path.abspath(os.path.dirname(superdark_resolved)) == os.path.abspath(ref_dir):
                            os.remove(superdark_resolved)
                            if verbose:
                                print 'Deleted old superdark:  {}'.format(superdark_resolved)
                        superdark_remakes.extend(file)
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
    anneals = archive_dark_query.archive_dark_query( \
                      raw_files, anneal_data=anneal_data, print_url=False)
    
    # Get list of EXPNAMEs from files in the dark_dir.
    found_dark_files = glob.glob(os.path.join(dark_dir, '*_flt.fits*'))
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
        # Populate file locations:
        for dark in anneal['darks']:
            if found.has_key(dark['exposure']):
                dark['file'] = found[dark['exposure']][0]
            else:
                missing_darks.add(dark['exposure'])
        
        # Update file headers:
        if crds_update:
            dark_files_to_run = [d['file'] for d in anneal['darks'] if d.has_key('file')]
            if len(dark_files_to_run) > 0:
                if verbose:
                    print 'Running crds.BestrefsScript on dark files from anneal {}...'.format(anneal)
                
                errors = BestrefsScript('BestrefsScript --update-bestrefs -s 1 -f ' + ' '.join(dark_files_to_run))()
                if int(errors) > 0:
                    raise Exception('CRDS BestrefsScript (darks):  Call returned errors!')
                
                # Update hdr0 data after BestrefsScript is run:
                for dark in anneal['darks']:
                    if dark.has_key('file'):
                        file = dark['file']
                        with fits.open(file) as f:
                            hdr0 = f[0].header
                        found[hdr0['ROOTNAME'].strip().upper()] = (file, hdr0)
        
        # For each amp, count weeks based on old darkfile; reset for each annealing period:
        # Get keywords about each dark:
        all_weeks = {'A':{}, 'B':{}, 'C':{}, 'D':{}}
        for dark in anneal['darks']:
            if found.has_key(dark['exposure']):
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
        if not ignore_missing:
            print 'Please download the missing darks (calibrated FLTs) via this link:'
            print '(or specify the proper dark_dir [{}])\n'.format(dark_dir)
            print 'If missing files are expected (e.g. amp=A darks), then run with'
            print '--ignore_missing flag.\n'
            print archive_dark_query.darks_url(missing_darks) + '\n'
            sys.exit(1)
        else:
            print '"--ignore_missing" flag set:  Ignoring these missing FLT files...'
    elif verbose:
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
            
            # Update ext=0 hdr of science file:
            # DARKFILE:
            darkfile = os.path.join(ref_dir, weekdark_tag + '_drk.fits')
            old_darkfile = f[0].header['DARKFILE']
            f[0].header['DARKFILE'] = darkfile
            
            f.flush()
            if verbose:
                print 'Updated hdr0 DARKFILE of {}:  {} --> {}'.format(file, old_darkfile, darkfile)
    
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
        weekdark_tags = weekdark.keys()  # Note the missing 's'
    
    # Generate specified basedarks and weekdarks, based on weekdark_tags:
    already_deleted = set()
    for weekdark_tag in weekdark_tags:
        # Make basedark:
        amp = weekdarks[weekdark_tag]['amp'].strip().lower()
        basedark = os.path.join(ref_dir, 'basedark_{}{}_drk.fits'.format( \
            weekdarks[weekdark_tag]['amp'].strip().lower(), weekdarks[weekdark_tag]['anneal_num']))
        anneal = [a for a in anneals if a['index'] == weekdarks[weekdark_tag]['anneal_num']][0]
        files = [f['file'] for f in anneal['darks'] if f.get('CCDAMP') == weekdarks[weekdark_tag]['amp']]
        if clean_all and basedark not in already_deleted:
            clean_this = True
            already_deleted.add(basedark)
        else:
            clean_this = False
        generate_basedark(files, basedark, pctetab, num_processes, clean_this, verbose)
        
        # Make weekdark:
        weekdark_name = os.path.abspath(os.path.expandvars(os.path.join(ref_dir, weekdark_tag + '_drk.fits')))
        files = [f['file'] for f in weekdarks[weekdark_tag]['darks']]  # Already selected for amp
        if clean_all and weekdark_name not in already_deleted:
            clean_this = True
            already_deleted.add(weekdark_name)
        else:
            clean_this = False
        generate_weekdark(files, weekdark_name, pctetab, basedark, num_processes, clean_this, verbose)


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
    '''
    Rename outputs according to output_mapping dictionary.
    
    Special commands for output products:
        '<pass>'   -- Don't do anything with product; denotes it as an intermediate product.
        '<remove>' -- Delete the intermediate product.
    '''
    
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
                    ext = ext[:-3]
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
            finally:
                os.chdir(cwd)


class Logger(object):
    '''Copies sys.stdout to a log file
       source: http://stackoverflow.com/a/24583265
       Modified to include STDERR.
    '''
    def __init__(self, filename='cti_{}.log'.format(datetime.datetime.now().isoformat('_')), mode="a", buff=0, disable=False):
        self.disable = disable
        
        if not self.disable:
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
        if not self.disable:
            self.stdout.write(message)  # Both STDOUT and STDERR get directed to STDOUT!
            self.file.write(message)
    
    def flush(self):
        if not self.disable:
            self.stdout.flush()
            self.stderr.flush()
            self.file.flush()
            os.fsync(self.file.fileno())
    
    def close(self):
        if not self.disable:
            if self.stdout != None:
                sys.stdout = self.stdout
                self.stdout = None
            
            if self.stderr != None:
                sys.stderr = self.stderr
                self.stderr = None
            
            if self.file != None:
                self.file.close()
                self.file = None
