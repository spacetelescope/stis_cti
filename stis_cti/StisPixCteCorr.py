"""
Functions to apply pixel-based CTE correction to STIS images.

The algorithm implemented in this code was
described in detail by [Anderson]_ as available
online at:

http://adsabs.harvard.edu/abs/2010PASP..122.1035A

.. note::
    * This code only works for STIS/CCD but can be modified to work on other
      detectors.
    * It was developed for use with full-frame GAIN={1,4} FLT images as input.
    * It has not been fully tested with any other formats.
    * Noise is slightly enhanced in the output (see [Anderson]_).
    * This code assumes a linear time dependence for a given set of
      coefficients.
    * This algorithm does not account for traps with very long release timescale
      but it is not an issue.
    * This code also does not account for second-exposure effect.
    * Multi-threading support was not implemented in this version.

:Optional preprocessing for nonstandard FLT input:

If you are not using a fully calibrated FLT image as input,
you might also need to do one or more of the following before running the task:

    * Convert image to unit of electrons.

    * For combined image (e.g., superdark), set `noise_model` to 0.

    * Primary FITS header (EXT 0) must have these keywords populated:
        * ROOTNAME
        * INSTRUME (must be STIS)
        * DETECTOR (must be CCD)
        * NEXTEND (number of extensions = 3 * number of image sets)
        * CCDAMP (ABCD, AD, BC, A, B, C, or D)
        * TEXPSTRT
        * ATODGAIN

    * SCI FITS header must have these keywords populated:
        * EXTNAME (must be SCI)
        * EXTVER (1 or greater)

    * ERR FITS header must have these keywords populated:
        * EXTNAME (must be ERR)
        * EXTVER (1 or greater)

    * DQ FITS header must have these keywords populated:
        * EXTNAME (must be DQ)
        * EXTVER (1 or greater)

:Example:

To correct a set of STIS FLT images, with one new CTE-corrected image
for each input.

>>> from stistools import StisPixCteCorr
>>> StisPixCteCorr.CteCorr('o*q_flt.fits')

:References:

.. [Anderson] Anderson J. & Bedin, L.R., 2010, PASP, 122, 1035

"""

#### TAKEN OUT FROM FORMAL DOC BUT KEPT FOR HISTORY
#
#:Authors: Pey Lian Lim and W.J. Hack (Python), J. Anderson (Fortran), Matt Davis
#
#:Organization: Space Telescope Science Institute
#
#:History:
#    * 2010/09/01 PLL created this module.
#    * 2010/10/13 WH added/modified documentations.
#    * 2010/10/15 PLL fixed PCTEFILE lookup logic.
#    * 2010/10/26 WH added support for multiple file processing
#    * 2010/11/09 PLL modified `YCte`, `_PixCteParams` and `_DecomposeRN` to reflect noise improvement by JA. Also updated documentations.
#    * 2011/04/26 MRD Began updates for new CTE algorithm.
#    * 2011/07/20 MRD Updated to handle new PCTETAB containing time dependent
#      CTE characterization.
#    * 2011/12/08 MRD Updated with latest CTE algorithm 3.0.
#    * 2012/01/31 MRD Updated with latest CTE algorithm 3.1.
#    * 2012/05/11 PLL updated to version 3.2 to be consistent with CALACS.
#    * 2012/05/22 PLL removed {} formatting for Python 2.5/2.6 compatibility.
####################

# External modules
import os, shutil, time
import numpy
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits

try:
    from stsci.tools import teal
except:
    teal = None

from stsci.tools import parseinput

# Local modules
import StisPixCte_FixY as pcfy # C extension


__taskname__ = "StisPixCteCorr"
__version__ = "0.1"
__vdate__ = "2013-Mar-20"

# constants related to the CTE algorithm in use
STIS_CTE_NAME = 'PixelCTE 2012'
STIS_CTE_VER = '3.2'

# number of extensions per image set (SCI, ERR, DQ)
EXTN_PER_IMSET = 3


# general error for things related to this module
class StisPixCteError(Exception):
    pass


#--------------------------
def CteCorr(input, outFits='', read_noise=None, noise_model=None,
            oversub_thresh=None, sim_nit=None, shift_nit=None):
    """
    Run all the CTE corrections on all the input files.

    This function simply calls `YCte` on each input image
    parsed from the `input` parameter, and passes all remaining
    parameter values through unchanged.

    Parameters
    ----------
    input : str or list of str
        name of FLT image(s) to be corrected. The name(s) can be specified
        either as:

          * a single filename ('o1234567q_flt.fits')
          * a Python list of filenames
          * a partial filename with wildcards ('\*flt.fits')
          * filename of an ASN table ('o12345670_asn.fits')
          * an at-file (``@input``)

    outFits : str
        *USE DEFAULT IF `input` HAS MULTIPLE FILES.*
        CTE corrected image in the same
        directory as input. If not given, will use
        ROOTNAME_cte.fits instead. Existing file will
        be overwritten.

    read_noise : float, optional
        Read noise level.
        If None, takes value from PCTETAB header RN_CLIP keyword.

    noise_model : {0, 1, 2, None}, optional
        Noise mitigation algorithm.
        If None, takes value from PCTETAB header NSEMODEL keyword.

        0: No smoothing
        1: Normal smoothing
        2: Strong smoothing

    oversub_thresh : float, optional
        Pixels corrected below this value will be re-corrected.
        If None, takes value from PCTETAB header SUBTHRSH keyword.

    sim_nit : int, optional
        Number of times readout simulation is performed per column.
        If None, takes value from PCTETAB header SIM_NIT keyword.

    shift_nit : int, optional
        Number of times column is shifted during simulated readout.
        If None, takes value from PCTETAB header SHFT_NIT keyword.

    """
    # Parse input to get list of filenames to process
    infiles, output = parseinput.parseinput(input)

    # Process each file
    for file in infiles:
        YCte(file, outFits=outFits, read_noise=read_noise,
             noise_model=noise_model, oversub_thresh=oversub_thresh,
             sim_nit=sim_nit, shift_nit=shift_nit)


#--------------------------
def YCte(inFits, outFits='', read_noise=None, noise_model=None,
         oversub_thresh=None, sim_nit=None, shift_nit=None):
    """
    Apply correction for parallel CTE loss.

    Input image that is already de-striped is desired
    but not compulsory. Using image with striping
    will enhance the stripes in output. Calibrations
    that have been applied to FLT should not
    significantly affect the result.

    Notes
    -----
    * EXT 0 header will be updated. ERR arrays will be
      added in quadrature with 10% of the correction.
      DQ not changed.

    * Does not work on RAW but can be modified
      to do so.

    Parameters
    ----------
    inFits : str
        FLT image to be corrected.

    outFits : str
        CTE corrected image in the same
        directory as input. If not given, will use
        ROOTNAME_cte.fits instead. Existing file will
        be overwritten.

    read_noise, noise_model, oversub_thresh, sim_nit, shift_nit : see `CteCorr`

    Examples
    --------
    Correct a single FLT image and write output to 'o12345678_cte.fits':

    >>> import StisPixCteCorr
    >>> StisPixCteCorr.YCte('o12345678_flt.fits')

    """
    if noise_model not in (0, 1, 2, None):
        raise ValueError("noise_model must be one of (0, 1, 2, None).")

    # Start timer
    timeBeg = time.time()

    # For output files naming.
    # Store in same path as input.
    outPath = os.path.dirname( os.path.abspath(inFits) ) + os.sep
    rootname = pyfits.getval(inFits, 'ROOTNAME')
    print 'Performing pixel-based CTE correction on', rootname
    rootname = outPath + rootname

    # Construct output filename
    if not outFits:
        outFits = rootname + '_cte.fits'

    # Copy input to output
    shutil.copyfile(inFits, outFits)

    # Open output for correction
    pf_out = pyfits.open(outFits, mode='update')

    # For checking that the detector is supported
    detector = pf_out['PRIMARY'].header['DETECTOR']

    atodgain = pf_out['PRIMARY'].header.get('ATODGAIN', default=1.)
    
    ccdamp = pf_out['PRIMARY'].header.get('CCDAMP', default='D')
    if ccdamp != 'D':
        print 'Non-standard amplifer %s being corrected!' % ccdamp
        
    if ccdamp == 'A' or ccdamp == 'B':
        readout_dir = 1
    elif ccdamp == 'C' or ccdamp == 'D':
        readout_dir = -1
    else:
        raise ValueError("Amplifier %s is unrecognized!" % ccdamp)

    # For epoch-specific operations
    expstart = pf_out['PRIMARY'].header['TEXPSTRT']

    nextend = pf_out['PRIMARY'].header.get('NEXTEND', default=EXTN_PER_IMSET)
    nimsets = nextend // EXTN_PER_IMSET
    if nextend != nimsets * EXTN_PER_IMSET:
        raise ValueError('Number of extensions is %d; must be a multiple'
                         ' of %d.' % (nextend, EXTN_PER_IMSET))

    # This is just for STIS CCD.
    if detector != 'CCD':
        raise StisPixCteError('Invalid detector: StisPixCteCorr only'
                              ' supports STIS CCD.')

    # Read CTE params from file
    pctefile = pf_out['PRIMARY'].header['PCTETAB']
    pardict = _PixCteParams(pctefile, expstart)

    cte_frac = pardict['cte_frac']
    q_dtde = pardict['q_dtde']
    dtde_l = pardict['dtde_l']
    psi_node = pardict['psi_node']
    chg_leak = pardict['chg_leak']
    levels = pardict['levels']
    col_scale = pardict['col_scale']
    pcte_ver = pardict['pcte_ver']  # version of PCTETAB reference file

    if sim_nit is None:
        sim_nit = pardict['sim_nit']

    if shift_nit is None:
        shft_nit = pardict['shft_nit']
    else:
        shft_nit = shift_nit

    if read_noise is None:
        read_noise = pardict['read_noise']

    if noise_model is None:
        noise_model = pardict['noise_model']

    if oversub_thresh is None:
        oversub_thresh = pardict['oversub_thresh']

    # N in charge tail
    chg_leak_kt, chg_open_kt = pcfy.InterpolatePsi(chg_leak, psi_node)
    del chg_leak, psi_node

    # dtde_q: Marginal PHI at a given chg level.
    dtde_q = pcfy.InterpolatePhi(dtde_l, q_dtde, shft_nit)
    del dtde_l, q_dtde

    # finish interpolation along the Q dimension and reduce arrays to contain
    # only info at the levels specified in the levels array
    chg_leak_lt, chg_open_lt, dpde_l = \
      pcfy.FillLevelArrays(chg_leak_kt, chg_open_kt, dtde_q, levels)
    del chg_leak_kt, chg_open_kt, dtde_q

    ########################################
    # perform correction for each image set
    ########################################
    for imset in range(nimsets):
        # get data
        extn = 1 + imset * EXTN_PER_IMSET
        # flip the data so the readout amp is closest to row 0, column 0
        # (this assumes STIS amp D was used)
        # Modified (SL):  Does readout_dir==+1 work for 'A' and 'B'?
        scidata = pf_out[extn].data[::readout_dir,::readout_dir].copy().astype(numpy.float)
        errdata = pf_out[extn+1].data[::readout_dir,::readout_dir].copy().astype(numpy.float)

        # convert to electrons
        scidata *= atodgain
        errdata *= atodgain

        # separate signal and noise
        if noise_model in (1,2):
            sigdata, nsedata = pcfy.DecomposeRN(scidata, read_noise,
                                                noise_model)
        elif noise_model == 0:
            sigdata = scidata.copy()
            nsedata = numpy.zeros_like(sigdata)

        # set up cte frac array; flip amp D info from col_scale
        col_scale_temp = col_scale[ccdamp][:scidata.shape[1]][::readout_dir]
        cte_frac_arr = cte_frac * numpy.ones_like(scidata) * col_scale_temp
        cte_frac_arr = cte_frac_arr * \
            numpy.arange(1,scidata.shape[0]+1).reshape((scidata.shape[0],1))
        cte_frac_arr /= 1024.

        # call CTE correction routine. data must be in units of electrons.
        print 'Performing CTE correction for science extension %d.' % extn

        t1 = time.time()
        cordata = pcfy.FixYCte(sigdata, sim_nit, shft_nit, oversub_thresh,
                               cte_frac_arr, levels, dpde_l, chg_leak_lt,
                               chg_open_lt)
        t2 = time.time()

        print 'FixYCte took {:.3f} seconds for science extension {:d}.'.format(t2-t1, extn)

        # add noise back in
        findata = cordata + nsedata

        # add error as 10% of change; note that units are still electrons
        delta = 0.1 * numpy.abs(scidata - findata)
        errdata = numpy.sqrt(errdata**2 + delta**2)

        # convert corrected data back to DN
        findata /= atodgain
        errdata /= atodgain

        # copy corrected data back to image.
        pf_out[extn].data[:,:] = findata.astype(numpy.float32)[::readout_dir,::readout_dir]

        pf_out[extn+1].data[:,:] = errdata.astype(numpy.float32)[::readout_dir,::readout_dir]
        # end for imset in range(nimsets)

    # Update header
    pf_out['PRIMARY'].header['PCTECORR'] = 'COMPLETE'
    pf_out['PRIMARY'].header['PCTEFRAC'] = (cte_frac, 'CTE time scaling value')
    pf_out['PRIMARY'].header['PCTERNCL'] = (read_noise, 'PCTE readnoise amplitude')
    pf_out['PRIMARY'].header['PCTENSMD'] = (noise_model, 'PCTE readnoise mitigation algorithm')
    pf_out['PRIMARY'].header['PCTETRSH'] = (oversub_thresh, 'PCTE over-subtraction threshold')
    pf_out['PRIMARY'].header['PCTESMIT'] = (sim_nit, 'PCTE readout simulation iterations')
    pf_out['PRIMARY'].header['PCTESHFT'] = (shft_nit, 'PCTE readout number of shifts')
    pf_out['PRIMARY'].header['CTE_NAME'] = (STIS_CTE_NAME, 'name of CTE algorithm')
    pf_out['PRIMARY'].header['CTE_VER'] = (STIS_CTE_VER, 'version of CTE algorithm')
    pf_out['PRIMARY'].header['PCTE_VER'] = (pcte_ver, 'version of PCTETAB')
    pf_out['PRIMARY'].header.add_history('PCTECORR complete ...')
    
    pf_out['PRIMARY'].header['ATODGAIN'] = atodgain
    
    # Close output file
    pf_out.close()
    print outFits, 'written'

    # Stop timer
    timeEnd = time.time()
    print 'Run time:  {:.3f} secs'.format(timeEnd - timeBeg)


#--------------------------
def _PixCteParams(fitsTable, expstart):
    """
    Read params from PCTEFILE.

    .. note:: Environment variable pointing to
              reference file directory must exist.

    Parameters
    ----------
    fitsTable : str
        PCTEFILE from header.

    expstart : float
        MJD of exposure start time, TEXPSTRT in image header

    Returns
    -------
    A dictionary containing the following:

    cte_frac : float
        Time dependent CTE scaling.

    sim_nit : int
        Number of readout simulations to do for each column of data

    shft_nit : int
        Number of shifts to break each readout simulation into

    read_noise : float
        Maximum amplitude of read noise removed by DecomposeRN.

    noise_model : {0, 1, 2}
        Read noise smoothing algorithm selection.

    oversub_thresh : float
        CTE corrected pixels taken below this value are re-corrected.

    dtde_q : ndarray
        Charge levels at which dtde_l is parameterized

    dtde_l : ndarray
        PHI(Q).

    psi_node : ndarray
        N values for PSI(Q,N).

    chg_leak : ndarray
        PSI(Q,N).

    levels : ndarray
        Charge levels at which to do CTE evaluation

    col_scale : dict of ndarray
        Dictionary containing the column-by-column CTE scaling for each amp.

    """

    # Resolve path to PCTEFILE
    refFile = _ResolveRefFile(fitsTable)
    if not os.path.isfile(refFile):
        raise IOError, 'PCTEFILE not found: %s' % refFile

    # Open FITS table
    pf_ref = pyfits.open(refFile)

    # Read RN_CLIP value from header
    read_noise = pf_ref['PRIMARY'].header['RN_CLIP']

    # read NSEMODEL value
    noise_model = pf_ref['PRIMARY'].header['NSEMODEL']

    # read SIM_NIT value from header
    sim_nit = pf_ref['PRIMARY'].header['SIM_NIT']

    # read SHFT_NIT value from header
    shft_nit = pf_ref['PRIMARY'].header['SHFT_NIT']

    # read SUBTHRSH value from header
    oversub_thresh = pf_ref['PRIMARY'].header['SUBTHRSH']

    # read number of CHG_LEAK# extensions from the header
    nchg_leak = pf_ref['PRIMARY'].header['NCHGLEAK']

    # read dtde data from DTDE extension
    dtde_l = pf_ref['DTDE'].data['DTDE']
    q_dtde = pf_ref['DTDE'].data['Q']

    # read levels data from LEVELS extension
    levels = pf_ref['LEVELS'].data['LEVEL']

    # read scale data from CTE_SCALE extension
    scalemjd = pf_ref['CTE_SCALE'].data['MJD']
    scaleval = pf_ref['CTE_SCALE'].data['SCALE']

    cte_frac = _CalcCteFrac(expstart, scalemjd, scaleval)

    # there are nchg_leak CHG_LEAK# extensions. we need to find out which one
    # is the right one for our data.
    chg_leak_names = ['CHG_LEAK%i'%(i) for i in range(1,nchg_leak+1)]

    for n in chg_leak_names:
        mjd1 = pf_ref[n].header['MJD1']
        mjd2 = pf_ref[n].header['MJD2']

        if (expstart >= mjd1) and (expstart < mjd2):
            # read chg_leak data from CHG_LEAK extension
            psi_node = pf_ref[n].data['NODE']
            chg_leak = numpy.array(pf_ref[n].data.tolist(), dtype=numpy.float32)[:,1:]
            break

    # column-by-column CTE scaling
    col_scale = {}

    col_scale['A'] = pf_ref['COL_SCALE'].data['AMPA']
    col_scale['B'] = pf_ref['COL_SCALE'].data['AMPB']
    col_scale['C'] = pf_ref['COL_SCALE'].data['AMPC']
    col_scale['D'] = pf_ref['COL_SCALE'].data['AMPD']

    pcte_ver = pf_ref['PRIMARY'].header['PCTE_VER']

    # Close FITS table
    pf_ref.close()

    d = {}
    d['cte_frac'] = cte_frac
    d['sim_nit'] = sim_nit
    d['shft_nit'] = shft_nit
    d['oversub_thresh'] = oversub_thresh
    d['read_noise'] = read_noise
    d['noise_model'] = noise_model
    d['q_dtde'] = q_dtde
    d['dtde_l'] = dtde_l
    d['psi_node'] = psi_node
    d['chg_leak'] = chg_leak
    d['levels'] = levels
    d['col_scale'] = col_scale
    d['pcte_ver'] = pcte_ver

    return d


#--------------------------
def _ResolveRefFile(refText, sep='$'):
    """
    Resolve the full path to reference file.
    This could be replaced with existing STSDAS
    library function, if necessary.

    Assume standard syntax: dir$file.fits

    Parameters
    ----------
    refText : str
        The text to process.

    sep : char
        Separator between directory and file name.

    Returns
    -------
    f : str
        Full path to reference file.
    """

    s = refText.split(sep)
    n = len(s)
    if n > 1:
        p = os.getenv(s[0])
        if p:
            p += os.sep
        else:
            p = ''
        # End if
        f = p + s[1]
    else:
        f = os.path.abspath(refText)
    # End if
    return f


#--------------------------
def _CalcCteFrac(expstart, scalemjd, scaleval):
    """
    Calculate CTE_FRAC used to scale CTE according to time dependence.

    Parameters
    ----------
    expstart : float
        TEXPSTRT from header.

    scalemjd : ndarray
        MJD points for corresponding CTE scale values in scaleval

    scaleval : ndarray
        CTE scale values corresponding to MJDs in scalemjd

    Returns
    -------
    cte_frac : float
        Time scaling factor.
    """

    # Calculate CTE_FRAC
    cte_frac = pcfy.CalcCteFrac(expstart, scalemjd, scaleval)

    return cte_frac


#--------------------------
def _InterpolatePsi(chg_leak, psi_node):
    """
    Interpolates the `PSI(Q,N)` curve at all N from
    1 to 100.

    `PSI(Q,N)` models the CTE tail profile across N
    pixels from the original pixel for a given
    charge, Q. Up to 100 pixels are tracked. For
    post-SM4 STIS/CCD, CTE loss is within 100 pixels.
    Might be worse for WFPC2 since it is older and
    has faster readout time.

    .. note:: As this model is refined, future release
              might only have PSI(N) independent of Q.

    Parameters
    ----------
    chg_leak : ndarray
        PSI table data from PCTEFILE.

    psi_node : ndarray
        PSI node data from PCTEFILE.

    Returns
    -------
    chg_leak : ndarray
        Interpolated PSI.

    chg_open : ndarray
        Interpolated tail profile data.

    """

    chg_leak, chg_open = pcfy.InterpolatePsi(chg_leak, psi_node.astype(numpy.int32))

    return chg_leak, chg_open


#--------------------------
def _InterpolatePhi(dtde_l, q_dtde, shft_nit):
    """
    Interpolates the `PHI(Q)` at all Q from 1 to
    99999 (log scale).

    `PHI(Q)` models the amount of charge in CTE
    tail, i.e., probability of an electron being
    grabbed by a charge trap.

    Parameters
    ----------
    dtde_l : ndarray
        PHI data from PCTEFILE.

    q_dtde : ndarray
        Q levels at which dtde_l is defined, read from PCTEFILE

    shft_int : int
        Number of shifts performed reading out CCD

    Returns
    -------
    dtde_q : ndarray
        dtde_l interpolated at all PHI levels

    """

    dtde_q = pcfy.InterpolatePhi(dtde_l, q_dtde, shft_nit)

    return dtde_q


def _FillLevelArrays(chg_leak, chg_open, dtde_q, levels):
    """
    Interpolates CTE parameters to the charge levels specified in levels.

    Parameters
    ----------
    chg_leak : ndarray
        Interpolated chg_leak tail profile data returned by _InterpolatePsi.

    chg_open : ndarray
        Interpolated chg_open tail profile data returned by _InterpolatePsi.

    dtde_q : ndarray
        PHI data interpolated at all PHI levels as returned by
        _InterpolatePhi.

    levels : ndarray
        Charge levels at which output arrays will be interpolated.
        Read from PCTEFILE.

    Returns
    -------
    chg_leak_lt : ndarray
        chg_leak tail profile data interpolated at the specified charge levels.

    chg_open_lt : ndarray
        chg_open tail profile data interpolated at the specified charge levels.

    dpde_l : ndarray
        dtde_q interpolated and summed for the specified charge levels.

    tail_len : ndarray
        Array of maximum tail lengths for the specified charge levels.

    """

    chg_leak_lt, chg_open_lt, dpde_l = \
      pcfy.FillLevelArrays(chg_leak, chg_open, dtde_q, levels)

    return chg_leak_lt, chg_open_lt, dpde_l


#--------------------------
def _DecomposeRN(data_e, read_noise=4.25, noise_model=1):
    """
    Separate noise and signal.

    REAL DATA = SIGNAL + NOISE

    Parameters
    ----------
    data_e : ndarray
        SCI data in electrons.

    read_noise : float, optional
        Read noise level.

    noise_model : {0, 1, 2}, optional
        Noise mitigation algorithm.

        0: No smoothing
        1: Normal smoothing
        2: Strong smoothing

    Returns
    -------
    sigArr : ndarray
        Noiseless signal component in electrons.

    nseArr : ndarray
        Noise component in electrons.

    """
    if noise_model not in (0, 1, 2):
        raise ValueError("noise_model must be one of (0, 1, 2).")

    sigArr, nseArr = pcfy.DecomposeRN(data_e, read_noise, noise_model)

    return sigArr, nseArr


def _FixYCte(detector, cte_data, sim_nit, shft_nit, oversub_thresh,
             cte_frac, levels, dpde_l, chg_leak_lt, chg_open_lt):
    """
    Perform CTE correction on input data. It is best to perform some kind
    of readnoise smoothing on the data, otherwise the CTE algorithm will
    amplify the read noise. (In the read out process readnoise is added to
    the data after CTE blurring.)

    Parameters
    ----------
    detector : str
        DETECTOR from header.
        Currently only 'CCD' is supported.

    cte_data : ndarray
        Data in need of CTE correction. For proper results cte_data[0,x] should
        be next to the readout register and cte_data[-1,x] should be furthest.
        Data are processed a column at a time, e.g. cte_data[:,x] is corrected,
        then cte_data[:,x+1] and so on.

    sim_nit : int
        Number of readout simulation iterations to perform.

    shft_nit : int
        Number of readout shifts to do.

    oversub_thresh : float
        CTE corrected pixels taken below this value are re-corrected.

    cte_frac : ndarray
        CTE scaling image combining the time dependent and column-by-column
        CTE factors. Should be same shape as `cte_data`.

    levels : ndarray
        Levels at which CTE is evaluated as read from PCTEFILE.

    dpde_l : ndarray
        Parameterized amount of charge in CTE trails as a function of
        specific charge levels, as returned by _FillLevelArrays.

    chg_leak_lt : ndarray
        Tail profile data at charge levels specified by levels, as returned
        by _FillLevelArrays.

    chg_open_lt : ndarray
        Tail profile data at charge levels specified by levels, as returned
        by _FillLevelArrays.

    Returns
    -------
    corrected : ndarray
        Data after CTE correction algorithm applied. Same size and shape as
        input cte_data.

    """

    if detector == 'CCD':
        corrected = pcfy.FixYCte(cte_data, sim_nit, shft_nit, oversub_thresh,
                                 cte_frac, levels, dpde_l,
                                 chg_leak_lt, chg_open_lt)
    else:
        raise StisPixCteError('Invalid detector: StisPixCteCorr only supports STIS CCD.')

    return corrected


def AddYCte(infile, outfile, shift_nit=None, units=None):
    """
    Add CTE blurring to input image using an inversion of the CTE correction
    code.

    .. note::

        No changes are made to the error or data quality arrays.

        Data should not have bias or prescan regions.

        Image must have PCTETAB, DETECTOR, TEXPSTRT, ATODGAIN,
        and NEXTEND header keywords.

    Parameters
    ----------
    infile : str
        Filename of image to be blurred. Should have the PCTETAB header
        keyword pointing to the PCTETAB reference file.

    outfile : str
        Filename of blurred output image.

    shift_nit : int, optional
        Number of times column is shifted during simulated readout.
        If None, takes value from PCTETAB header SHFT_NIT keyword.

    units : {None,'electrons','counts'}, optional
        If 'electrons', the input image is assumed to have units of electrons
        and no gain operations are performed.
        If 'counts', the data are assumed to be in DN and they are converted
        to electrons before CTE blurring is performed. The ATODGN* keywords
        from the primary header are used for the conversions.
        If None, the BUNIT keyword from the science extension headers is used
        to set the unit behavior.
        Defaults to None.
        For STIS data, this argument is ignored; the units are assumed
        to be counts, so the data will be converted to electrons.

    Raises
    ------
    ValueError
        If the number of extensions is not a multiple of three.

    stistools.StisPixCteCorr.StisPixCteError
        If the input image comes from an imcompatible detector.

    """

    # copy infile to outfile
    shutil.copyfile(infile, outfile)

    # open file for blurring
    pf_out = pyfits.open(outfile, mode='update')

    # For checking that the detector is supported
    detector = pf_out['PRIMARY'].header['DETECTOR']

    atodgain = pf_out['PRIMARY'].header.get('ATODGAIN', default=1.)

    # For epoch-specific operations
    expstart = pf_out['PRIMARY'].header['TEXPSTRT']
    
    nextend = pf_out['PRIMARY'].header.get('NEXTEND', default=EXTN_PER_IMSET)
    nimsets = nextend // EXTN_PER_IMSET
    if nextend != nimsets * EXTN_PER_IMSET:
        raise ValueError('Number of extensions is %d; must be a multiple'
                         ' of %d.' % (nextend, EXTN_PER_IMSET))

    # This is just for CCD.
    if detector != 'CCD':
        os.remove(outfile)
        raise StisPixCteError('Invalid detector: StisPixCteCorr only'
                              ' supports STIS CCD.')

    # Read CTE params from file
    pctefile = pf_out['PRIMARY'].header['PCTETAB']
    pardict = _PixCteParams(pctefile, expstart)

    cte_frac = pardict['cte_frac']
    q_dtde = pardict['q_dtde']
    dtde_l = pardict['dtde_l']
    psi_node = pardict['psi_node']
    chg_leak = pardict['chg_leak']
    levels = pardict['levels']
    col_scale = pardict['col_scale']

    if shift_nit is None:
        shft_nit = pardict['shft_nit']
    else:
        shft_nit = shift_nit

    # The following are only used to update the header of the output, blurred image
    sim_nit = pardict['sim_nit']
    rn_clip = pardict['read_noise']

    # N in charge tail
    chg_leak_kt, chg_open_kt = pcfy.InterpolatePsi(chg_leak, psi_node)
    del chg_leak, psi_node

    # dtde_q: Marginal PHI at a given chg level.
    # q_pix_array: Maps Q (cumulative charge) to P (dependent var).
    # pix_q_array: Maps P to Q.
    dtde_q = pcfy.InterpolatePhi(dtde_l, q_dtde, shft_nit)
    del dtde_l, q_dtde

    # finish interpolation along the Q dimension and reduce arrays to contain
    # only info at the levels specified in the levels array
    chg_leak_lt, chg_open_lt, dpde_l = \
      pcfy.FillLevelArrays(chg_leak_kt, chg_open_kt, dtde_q, levels)
    del chg_leak_kt, chg_open_kt, dtde_q
    
    ccdamp = pf_out['PRIMARY'].header.get('CCDAMP', default='D')
    if ccdamp == 'A' or ccdamp == 'B':
        readout_dir = 1
    elif ccdamp == 'C' or ccdamp == 'D':
        readout_dir = -1
    else:
        raise ValueError("Amplifier %s is unrecognized!" % ccdamp)
    
    ########################################
    # perform CTE blurring for each image set
    ########################################
    for imset in range(nimsets):
        # get data
        extn = 1 + imset * EXTN_PER_IMSET
        # [::-1,::-1] flips the data so that pixel [0, 0] is closest to
        # the readout, since that's what the algorithm expects
        # (this assumes STIS amp D was used)
        # Modified (SL):  Does readout_dir==+1 work for 'A' and 'B'?
        scidata = pf_out[extn].data[::readout_dir,::readout_dir].copy().astype(numpy.float)
        # convert to electrons
        scidata *= atodgain

        # set up cte frac array; flip amp D info from col_scale
        col_scale_temp = col_scale[ccdamp][:scidata.shape[1]][::readout_dir]
        cte_frac_arr = cte_frac * numpy.ones_like(scidata) * col_scale_temp
        cte_frac_arr = cte_frac_arr * \
            numpy.arange(1,scidata.shape[0]+1).reshape((scidata.shape[0],1))
        cte_frac_arr /= 1024.

        # call CTE blurring routine. data must be in units of electrons.
        print 'Performing CTE blurring for science extension %d.' % extn

        t1 = time.time()
        cordata = _AddYCte(detector, scidata, cte_frac_arr, shft_nit,
                           levels, dpde_l, chg_leak_lt, chg_open_lt)
        t2 = time.time()

        print 'AddYCte took {:.3f} seconds for science extension {:d}.'.format(t2-t1, extn)

        # convert blurred data back to DN
        cordata /= atodgain

        # copy blurred data back to image, flipping to its original orientation
        pf_out[extn].data[:,:] = cordata.astype(numpy.float32)[::readout_dir,::readout_dir]
        # end for imset in range(nimsets)

    # Update header
    pf_out['PRIMARY'].header['PCTEFRAC'] = cte_frac
    pf_out['PRIMARY'].header['PCTERNCL'] = rn_clip
    pf_out['PRIMARY'].header['PCTESMIT'] = sim_nit
    pf_out['PRIMARY'].header['PCTESHFT'] = shft_nit
    pf_out['PRIMARY'].header.add_history('CTE blurring performed by PixCteCorr.AddYCte')

    # close image
    pf_out.close()


def _AddYCte(detector, input_data, cte_frac, shft_nit, levels, dpde_l,
             chg_leak_lt, chg_open_lt):
    """
    Apply STIS CTE blurring to input data.

    Parameters
    ----------
    detector : str
        DETECTOR from header.
        Currently only 'CCD' is supported.

    input_data : ndarray
        Data in need of CTE correction. For proper results cte_data[0,x] should
        be next to the readout register and cte_data[-1,x] should be furthest.
        Data are processed a column at a time, e.g. cte_data[:,x] is corrected,
        then cte_data[:,x+1] and so on.

    cte_frac : ndarray
        CTE scaling image combining the time dependent and column-by-column
        CTE factors. Should be same shape as `cte_data`.

    shft_nit : int
        Number of readout shifts to do.

    levels : ndarray
        Levels at which CTE is evaluated as read from PCTEFILE.

    dpde_l : ndarray
        Parameterized amount of charge in CTE trails as a function of
        specific charge levels, as returned by _FillLevelArrays.

    chg_leak_lt : ndarray
        Tail profile data at charge levels specified by levels, as returned
        by _FillLevelArrays.

    chg_open_lt : ndarray
        Tail profile data at charge levels specified by levels, as returned
        by _FillLevelArrays.

    Returns
    -------
    blurred : ndarray
        Data CTE correction algorithm applied.
        Same size and shape as input_data.

    """

    if detector == 'CCD':
        blurred = pcfy.AddYCte(input_data, shft_nit, cte_frac,
                                levels, dpde_l, chg_leak_lt, chg_open_lt)
    else:
        raise StisPixCteError('Invalid detector: StisPixCteCorr only supports STIS CCD.')

    return blurred


#--------------------------
# TEAL Interface functions
#--------------------------
def run(configObj):
    """
    TEAL interface for the `CteCorr` function.

    """
    CteCorr(configObj['inFits'],
            outFits=configObj['outFits'],
            read_noise=configObj['read_noise'],
            noise_model=configObj['noise_model'],
            oversub_thresh=configObj['oversub_thresh'],
            sim_nit=configObj['sim_nit'],
            shift_nit=configObj['shift_nit'])


def getHelpAsString():
    """
    Returns documentation on the `CteCorr` function. Required by TEAL.

    """
    return CteCorr.__doc__
