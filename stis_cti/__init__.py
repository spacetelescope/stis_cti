'''
Usage:
import stis_cti
stis_cti.stis_cti(science_dir, dark_dir, ref_dir, pctetab, num_processes, 
                  all_weeks_flag=False, allow=False, clean=False, verbose=False)

...
'''

from stis_cti import stis_cti, __version__, __author__
from stis_cti import viable_ccd_file
import StisPixCteCorr
#from StisPixCteCorr import AddYCte, CteCorr, YCte
#from StisPixCteCorr import __version__ as __StisPixCteCorr_version__
#from StisPixCteCorr import __vdate__ as __StisPixCteCorr_vdate__
from archive_dark_query import archive_dark_query

# These lines allow TEAL to print out the names of TEAL-enabled tasks
# upon importing this package.
#import os
#from stsci.tools import teal
#teal.print_tasknames(__name__, os.path.dirname(__file__))
