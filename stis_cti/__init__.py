'''
Usage:
import stis_cti
stis_cti.stis_cti(science_dir, dark_dir, ref_dir, num_processes, pctetab=None, 
                  all_weeks_flag=False, allow=False, clean=False, clean_all=False, verbose=1)
'''

from .stis_cti import (
    stis_cti, __author__, viable_ccd_file, determine_input_science, resolve_iraf_file,
    check_pctetab_version, superdark_hash, check_for_old_output_files,
    VersionError, FileError
)
from .custom_superdark_info import custom_superdark_info
from . import StisPixCteCorr
from .archive_dark_query import archive_dark_query
try:
    from ._version import __version__
except ImportError:
    __version__ = 'unknown'
