#!/usr/bin/env python

import sys
import pickle
from astropy.io import fits
import glob
from collections import Counter
import datetime
from numpy import size, shape
from six.moves.urllib import request as urlrequest, parse as urlparse, error as urlerror

__author__  = 'Sean Lockwood'
__version__ = '0.3.0'


# Data container for our dark exposure search results:
class darkType(dict):
    def __str__(self):
        try:
            return self['exposure']
        except KeyError:
            return 'undefined'
    def __repr__(self):
        return '<' + self.__str__() + '>'


class annealType(dict):
    def __str__(self):
        try:
            return '<' + '{:03d}'.format(self['index']) + '; ' + \
                   self['start'].isoformat(' ') + ' - ' + self['end'].isoformat(' ') + '>'
        except KeyError:
            return '<undefined>'
    def __repr__(self):
        return self.__str__()


def get_proposal_ids(abstract='stis, +ccd', title='dark, +monitor'):
    '''
    Perform an HST Abstract search for STIS/CCD Dark Monitor programs.
    '''
    url = 'https://archive.stsci.edu/hst/abstract.html'
    form_data = urlparse.urlencode({
        'abstract': abstract,      # String to be searched for within the abstract
        'atitle':   title,         # String to be searched for within the title
        'checkbox': 'no',          # Display abstract?
        'submit':   'submit'}).encode()

    # Submit POST request:
    try:
        response = urlrequest.urlopen(url, form_data)  # Needs better error handling!
        lines = response.readlines()
        lines = [line.decode('utf-8') for line in lines]
        text = [line.rsplit('\n',1)[0] for line in lines]
        response.close()
    except (urlerror.URLError, urlerror.HTTPError) as e:
        print('Please check your internet connection!')
        raise e

    # Parse returned text for proposal IDs:
    text = [x for x in text if 'proposal_search.php' in x]
    proposals = [int(y.split('>')[1].split('<')[0]) for y in text]

    return proposals


def read_dark_exposures():
    '''
    Query the HST archive for STIS/CCD dark exposures.
    '''
    # URL for HTTP GET request to the HST archive:
    url = 'https://archive.stsci.edu/hst/search.php?' + urlparse.urlencode([
        ('sci_instrume'         , 'STIS'                                         ),
        ('sci_instrument_config', 'STIS/CCD'                                     ),
        ('sci_targname'         , 'DARK'                                         ),
        ('sci_aec'              , 'C'                                            ),
        ('resolve'              , 'don%27tresolve'                               ),
        ('max_records'          , '50000'                                        ),
        ('selectedColumnsCsv'   , 'sci_data_set_name,sci_pep_id,sci_start_time,sci_actual_duration'),
        ('ordercolumn1'         , 'sci_start_time'                               ),
        ('showquery'            , 'off'                                          ),
        ('skipformat'           , 'on'                                           ),
        ('nonull'               , 'on'                                           ),
        ('outputformat'         , 'CSV'                                          ),
        ('action'               , 'Search'                                       )])

    # Also interested in:
    #    stis_ref_data..ssr_ccdgain
    #    stis_ref_data..ssr_ccdamp
    #    stis_ref_data..ssr_ccdoffst

    # Submit HTTP GET request:
    try:
        response = urlrequest.urlopen(url)
        data = response.readlines()
        response.close()
    except (urlerror.URLError, urlerror.HTTPError) as e:
        print('Please check your internet connection!')
        raise e

    darks = []
    for line in data:
        line = line.decode('utf-8')
        if line != '\n' and 'Dataset' not in line and 'string' not in line:
            exposure, proposid, time, exptime = line.split(',')
            exptime = float(exptime.split('\n')[0])
            proposid = int(proposid)
            time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
            darks.append(darkType(exposure=exposure, proposid=proposid, datetime=time, exptime=exptime))

    return darks


def get_anneal_boundaries(delta_days=5, min_exptime=None, verbose=False):
    '''
    Determines the annealing period boundaries and all associated darks within.
    '''

    if min_exptime is None:
        min_exptime = 960.0  # seconds
    else:
        min_exptime = float(min_exptime)

    if verbose:
        print('Minimum exposure time = {} s\n'.format(min_exptime))

    print('Querying MAST archive for dark and anneal program IDs...')
    dark_programs = get_proposal_ids()
    dark_programs.extend([7092, 7601, 7802, 7926, 7948, 7949])  # Older programs that didn't match this search pattern
    dark_programs = list(set(dark_programs))
    dark_programs.sort()

    anneal_programs = get_proposal_ids(abstract='', title='STIS CCD Hot Pixel Annealing')
    anneal_programs.sort()

    if verbose:
        print('\nDark proposal IDs:       ')
        print(', '.join([str(x) for x in dark_programs]))
        print('\nAnnealing proposal IDs:  ')
        print(', '.join([str(x) for x in anneal_programs]))
        print()

    print('Querying MAST archive for darks...') 
    darks = read_dark_exposures()

    print('Parsing archive results...\n')

    anneal_exposures = [x for x in darks if x['proposid'] in anneal_programs]

    # Find delta time between neighboring anneal exposures:
    anneal_start_boundary = [anneal_exposures[0]['datetime'] - datetime.timedelta(days=60)]
    anneal_end_boundary = [anneal_exposures[0]['datetime']]  # Or, use a later exposure < delta_days?

    previous_anneal = anneal_exposures[0]
    for anneal_exposure in anneal_exposures[1:]:
        delta = anneal_exposure['datetime'] - previous_anneal['datetime']
        if (delta.total_seconds() / 3600. / 24.) >= delta_days:
            anneal_start_boundary.append(previous_anneal['datetime'])
            anneal_end_boundary.append(anneal_exposure['datetime'])  # Or, use start time to avoid gaps?
            previous_anneal = anneal_exposure

    # Add an additional annealing period after the last anneal (end date is approximate):
    anneal_start_boundary.append(anneal_exposure['datetime'])
    anneal_end_boundary.append(anneal_exposure['datetime'] + datetime.timedelta(days=60))

    # Determine which darks fall within each annealing period:
    anneals = []

    for i, (a, b) in enumerate(zip(anneal_start_boundary, anneal_end_boundary)):
        darks_within_anneal = [x for x in darks if \
                      x['datetime'] >= a and 
                      x['datetime'] < b  and
                      x['proposid'] in dark_programs and
                      x['exptime'] >= min_exptime]
        
        anneals.append(annealType(index=i, start=a, end=b, darks=list(darks_within_anneal)))

    if (verbose >= 2):
        # Print information about all annealing periods:
        print("Dark data found for all annealing periods:")
        for i, a in enumerate(anneals):
            print(('{i:3d} {start:} - {end:} ({length:4d} days), {total_exposures:3d} total: ' + \
                   '{exptimes} {proposals}').format( 
                i = i, 
                start = a['start'].isoformat(' '), 
                end = a['end'].isoformat(' '), 
                length = (a['end'] - a['start']).days, 
                total_exposures = shape(a['darks'])[0], 
                exptimes = Counter([int(x['exptime']) for x in a['darks']]), 
                proposals = list(set([x['proposid'] for x in a['darks']]))))
        print()

    return anneals


def darks_url(exposures):
    # Build a URL to retrieve matched datasets from MAST:
    url = 'https://archive.stsci.edu/hst/search.php?' + urlparse.urlencode([
        ('sci_instrume'         , 'STIS'              ),
        ('sci_instrument_config', 'STIS/CCD'          ),
        ('sci_targname'         , 'DARK'              ),
        ('sci_aec'              , 'C'                 ),
        ('resolve'              , 'don\'tresolve'     ),
        ('sci_data_set_name'    , ','.join(exposures) ),
        ('max_records'          , '50000'             ),
        ('max_rpp'              , '5000'              ),
        ('ordercolumn1'         , 'sci_start_time'    ),
        ('action'               , 'Search'            )])

    return url


def archive_dark_query(files, anneal_data=None, min_exptime=None, verbose=False, print_results=True):
    '''
    Queries the MAST archive to determine which component darks are needed to 
    generate a CTI-corrected super-dark.

    :param files: list of str, str
        Science files to re-reduce (or search string).
    :param anneal_data:
        Previously saved anneal boundary data, as restored from a Pickle file.
    :param min_exptime: float or None
        Minimum exposure time of component darks (in seconds) [default=960]
    :param verbose: bool
        Print more information about matching annealing periods.
    :param print_url: bool
        Print data retrieval URL [default=Truee]

    :returns: list of matching dark rootnames
    '''
    # Convert to a list of filenames if only a string is provided:
    if isinstance(files, str):
        files = [files]

    # Apply glob to every element of files to evaluate search patterns:
    files_parsed = []
    for f in files:
        files_parsed.extend(glob.glob(f))
    files = files_parsed

    if verbose:
        print('Input files:')
        for f in files:
            print(f)
        print()

    if len(files) == 0:
        raise IOError('No files specified/found.')

    # Determine time boundaries for each annealing period:
    # (This is costly, so pass the value into this function when calling multiple times.)
    if anneal_data == None:
        anneal_data = get_anneal_boundaries(verbose=verbose, min_exptime=min_exptime)
    else:
       if verbose:
           print('Skipping re-population of anneal_data.\n')

    # Populate matched exposures uniquely:
    matches = list()

    for file in files:
        with fits.open(file) as data:
            hdr0 = data[0].header
            if hdr0['INSTRUME'].strip() != 'STIS':
                raise ValueError('Not a STIS file!')
            date = hdr0['TDATEOBS']
            time = hdr0['TTIMEOBS']
            dt = datetime.datetime.strptime(date + ' ' + time, '%Y-%m-%d %H:%M:%S')
            match = [x for x in anneal_data if dt >= x['start'] and dt < x['end']][0]

            #matches.add(match)  # WON'T WORK WITH MUTABLE SUB-TYPES! ***
            matches.append(match)

            if verbose:
                print(file)
                print('Annealing period:  ', match)
                if verbose >= 2:
                    print(', '.join([str(i) for i in match['darks']]))
                print()

    # Remove duplicate entries based on the index field:
    matches = list(dict((item['index'], item) for item in matches).values())
    matches.sort(key=lambda x: x['index'])

    if verbose:
        print('Unique annealing periods (' + str(shape(matches)[0]) + '):')
        print('\n'.join([str(i) for i in matches]))
        print()

    # Get the unique list of exposure names:
    all_exposures = set()
    for m in matches:
        for file in m['darks']:
            all_exposures.add(file['exposure'])
    all_exposures = list(all_exposures)
    all_exposures.sort()

    if verbose:
        print('Unique exposures (' + str(shape(all_exposures)[0]) + '):')
        print(', '.join(all_exposures))
        print()

    if print_results:
        print(', '.join(all_exposures))
        print('\nDownload dark FLT files via Astroquery or MAST form:')
        print('https://mast.stsci.edu/search/ui/#/hst')
        print()

    return matches


def call_archive_dark_query():
    import argparse

    parser = argparse.ArgumentParser( 
        description='Determines which STIS/CCD darks are necessary to regenerate a ' + 
            'science dataset\'s super-dark.', 
        epilog='Author:   ' + __author__ + '; ' + 
               'Version:  ' + __version__)
    parser.add_argument(dest='file', action='store', nargs='+', 
                        help='Science files to re-reduce (or search string).')
    parser.add_argument('-t', dest='min_exptime', action='store', default=None, 
                        help='Minimum exposure time of component darks (in seconds) [default=960]')
    parser.add_argument('-p', dest='pickle_file', action='store', default=None, 
                        help='Optional pickle file containing anneal data (skip querying MAST)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', 
                        help='Print more information about matching annealing periods.')
    parser.add_argument('-vv', dest='very_verbose', action='store_true', 
                        help='Very verbose')
    args = parser.parse_args()

    verbose = args.verbose

    if args.very_verbose:
        verbose = 2

    if args.pickle_file != None:
        if args.min_exptime is not None:
            print('Error!  Pickle file does not support non-default min_exptime.')
            sys.exit()
        anneal_data = pickle.load(open(args.pickle_file, 'rb'))
        if verbose:
            print('Reading anneal_data from pickle file:  {}\n'.format(args.pickle_file))
    else:
        anneal_data = None

    if len(args.file) <= 0:
        print('Error:  No files found!')
        sys.exit()

    tmp = archive_dark_query(args.file, anneal_data=anneal_data, 
        min_exptime=args.min_exptime, verbose=verbose, print_results=True)


if __name__ == '__main__':
    call_archive_dark_query()
