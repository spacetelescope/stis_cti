'''
Usage:
import stis_cti
stis_cti.stis_cti(science_dir, dark_dir, ref_dir, num_processes, pctetab=None, 
                  all_weeks_flag=False, allow=False, clean=False, clean_all=False, verbose=1)
'''

from stis_cti import stis_cti, __version__, __author__
from stis_cti import viable_ccd_file
from custom_superdark_info import custom_superdark_info
import StisPixCteCorr
from archive_dark_query import archive_dark_query
