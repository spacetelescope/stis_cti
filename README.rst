HST/STIS Pixel-Based CTI-Correction Scripts

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

See documentation at http://pythonhosted.org/stis_cti/

Project page:  http://www.stsci.edu/hst/stis/software/analyzing/scripts/pixel_based_CTI/
