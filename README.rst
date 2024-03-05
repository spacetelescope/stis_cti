HST/STIS Pixel-Based CTI-Correction Scripts
===========================================

Utilities needed to correct for Charge Transfer Inefficiency (CTI) in the Hubble
Space Telescope (HST) STIS CCD.

Scripts installed in shell:

  * stis_cti  -- Runs CTI correction on raw files, handling super-dark creation
  * archive_dark_query -- Determines component darks needed to remake a super-dark

For more information, type:

  ``stis_cti --help``

Python usage::

  import stis_cti
  stis_cti.stis_cti(...)

Other utilities:

  * stis_cti.StisPixCteCorr(...) -- Code to run pixel-based CTI-correction on intermediate products
  * stis_cti.viable_ccd_file(...) -- Test to see which FITS files on which to run the correction
  * stis_cti.archive_dark_query(...) -- Utility to query MAST for needed component darks
  * stis_cti.custom_superdark_info() -- Print information on manually creating/implementing a corrected super-dark

See documentation at https://stis-cti.readthedocs.io

Project page:  http://www.stsci.edu/hst/instrumentation/stis/data-analysis-and-software-tools/pixel-based-cti

Code repository:  https://github.com/spacetelescope/stis_cti
