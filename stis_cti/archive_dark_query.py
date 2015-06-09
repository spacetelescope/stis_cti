from astropy.io import fits
import glob
from collections import Counter
import datetime
import urllib, urllib2
from numpy import size, shape

__author__  = 'Sean Lockwood'
__version__ = '0.1.2'


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
    
    url = 'http://archive.stsci.edu/hst/abstract.html'
    form_data = urllib.urlencode({
        'abstract': abstract,      # String to be searched for within the abstract
        'atitle':   title,         # String to be searched for within the title
        'checkbox': 'no',          # Display abstract?
        'submit':   'submit'})
    
    # Submit POST request:
    try:
        response = urllib2.urlopen(url, form_data)  # Needs better error handling!
        text = response.read().split('\n')
        response.close()
    except (urllib2.URLError, urllib2.HTTPError) as e:
        print 'Please check your internet connection!'
        raise e
    
    # Parse returned text for proposal IDs:
    text = filter(lambda x: 'proposal_search.php' in x, text)
    proposals = map(lambda y: int(y.split('>')[1].split('<')[0]), text)
    
    return proposals


def read_dark_exposures():
    '''
    Query the HST archive for STIS/CCD dark exposures.
    '''
    
    # URL for HTTP GET request to the HST archive:
    url = 'http://archive.stsci.edu/hst/search.php?' + urllib.urlencode([
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
        response = urllib2.urlopen(url)
        data = response.readlines()
        response.close()
    except (urllib2.URLError, urllib2.HTTPError) as e:
        print 'Please check your internet connection!'
        raise e
    
    darks = []
    for line in data:
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
        print 'Minimum exposure time = ', min_exptime, ' s'
        print
    
    print 'Querying MAST archive for dark and anneal program IDs...'
    dark_programs = get_proposal_ids()
    dark_programs.extend([7092, 7601, 7802, 7926, 7948, 7949])  # Older programs that didn't match this search pattern
    dark_programs = list(set(dark_programs))
    dark_programs.sort()
    
    anneal_programs = get_proposal_ids(abstract='', title='STIS CCD Hot Pixel Annealing')
    anneal_programs.sort()
    
    if verbose:
        print
        print 'Dark proposal IDs:       '
        print ', '.join(map(lambda x: str(x), dark_programs))
        print
        print 'Annealing proposal IDs:  '
        print ', '.join(map(lambda x: str(x), anneal_programs))
        print
    
    print 'Querying MAST archive for darks...' 
    darks = read_dark_exposures()
    
    print 'Parsing archive results...'
    print
    
    anneal_exposures = filter(lambda x: x['proposid'] in anneal_programs, darks)
    
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
        darks_within_anneal = filter( \
            lambda x: x['datetime'] >= a and 
                      x['datetime'] < b  and
                      x['proposid'] in dark_programs and
                      x['exptime'] >= min_exptime,
            darks)
        
        anneals.append(annealType(index=i, start=a, end=b, darks=list(darks_within_anneal)))
    
    if (verbose >= 2):
        # Print information about all annealing periods:
        print "Dark data found for all annealing periods:"
        for i, a in enumerate(anneals):
            print ('{i:3d} {start:} - {end:} ({length:4d} days), {total_exposures:3d} total: ' + \
                   '{exptimes} {proposals}').format( 
                i = i, 
                start = a['start'].isoformat(' '), 
                end = a['end'].isoformat(' '), 
                length = (a['end'] - a['start']).days, 
                total_exposures = shape(a['darks'])[0], 
                exptimes = Counter(map(lambda x: int(x['exptime']), a['darks'])), 
                proposals = list(set(map(lambda x: x['proposid'], a['darks']))))
        print
    
    return anneals


def darks_url(exposures):
    # Build a URL to retrieve matched datasets from MAST:
    url = 'http://archive.stsci.edu/hst/search.php?' + urllib.urlencode([
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


def archive_dark_query(files, anneal_data=None, min_exptime=None, verbose=False, print_url=True):
    '''
    Queries the MAST archive to determine which component darks are needed to 
    generate a CTI-corrected super-dark.
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
        print 'Input files:'
        for f in files:
            print f
        print
    
    if len(files) == 0:
        raise IOError('No files specified/found.')
    
    # Determine time boundaries for each annealing period:
    # (This is costly, so pass the value into this function when calling multiple times.)
    if anneal_data == None:
        anneal_data = get_anneal_boundaries(verbose=verbose, min_exptime=min_exptime)
    else:
       if verbose:
           print 'Skipping re-population of anneal_data.'
           print
    
    # Populate matched exposures uniquely:
    matches = list()
    
    for file in files:
        with fits.open(file) as data:
            hdr0 = data[0].header
            date = hdr0['TDATEOBS']
            time = hdr0['TTIMEOBS']
            dt = datetime.datetime.strptime(date + ' ' + time, '%Y-%m-%d %H:%M:%S')
            match = filter(lambda x: dt >= x['start'] and dt < x['end'], anneal_data)[0]
            
            #matches.add(match)  # WON'T WORK WITH MUTABLE SUB-TYPES! ***
            matches.append(match)
            
            if verbose:
                print file
                print 'Annealing period:  ', match
                if verbose >= 2:
                    print ', '.join([str(i) for i in match['darks']])
                print
    
    # Remove duplicate entries based on the index field:
    matches = dict((item['index'], item) for item in matches).values()
    matches.sort(key=lambda x: x['index'])
    
    if verbose:
        print 'Unique annealing periods (' + str(shape(matches)[0]) + '):'
        print '\n'.join([str(i) for i in matches])
        print
    
    # Get the unique list of exposure names:
    all_exposures = set()
    for m in matches:
        for file in m['darks']:
            all_exposures.add(file['exposure'])
    all_exposures = list(all_exposures)
    all_exposures.sort()
    
    if verbose:
        print 'Unique exposures (' + str(shape(all_exposures)[0]) + '):'
        print ', '.join(all_exposures)
        print
    
    if print_url:
        url = darks_url(all_exposures)
        print 'Download darks via this link:'
        print
        print url
        print
    
    return matches
