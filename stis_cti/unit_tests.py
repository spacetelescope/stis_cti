#!/usr/bin/env python

from nose.tools import assert_equals, assert_false, assert_raises
import os
import datetime
from astropy.io import fits

#from stis_cti import resolve_iraf_file, viable_ccd_file, superdark_hash
from stis_cti import *


class TestPaths:
    '''
    Tests functionality of stis_cti.resolve_iraf_file
    '''
    def setup(self):
        os.environ['toref'] = '/grp/hst/cdbs/oref/'
    
    def teardown(self):
        os.environ.pop('toref')
    
    def test_no_dollar(self):
        assert_equals(resolve_iraf_file('filename.fits'), 'filename.fits')
        assert_equals(resolve_iraf_file('/dir/filename.fits'), '/dir/filename.fits')
        assert_equals(resolve_iraf_file('dir/filename.fits'), 'dir/filename.fits')
    
    def test_dollar_at_beginning(self):
        assert_equals(resolve_iraf_file('$toref/filename.fits'), '/grp/hst/cdbs/oref/filename.fits')
        assert_equals(resolve_iraf_file('$toref/path/filename.fits'), '/grp/hst/cdbs/oref/path/filename.fits')
    
    def test_dollar_in_middle(self):
        assert_equals(resolve_iraf_file('toref$filename.fits'), '/grp/hst/cdbs/oref/filename.fits')
        assert_equals(resolve_iraf_file('toref$path/filename.fits'), '/grp/hst/cdbs/oref/path/filename.fits')
    
    def test_environ_var_undefined(self):
        with assert_raises(IOError) as cm:
            resolve_iraf_file('undoref$filename.fits')
        ex = cm.exception
        assert_equals(ex.message, "Can't resolve environmental variable 'undoref'.")


def setup_file(test_file):
    # Create and write out a test FITS file:
    if not os.access(os.path.curdir, os.W_OK):
        raise Exception('Can\'t write test file to CWD!')
    
    hdu = fits.PrimaryHDU()
    hdu.header['TELESCOP'] = ('HST',           'telescope used to acquire data')
    hdu.header['INSTRUME'] = ('STIS',          'identifier for instrument used to acquire data')
    hdu.header['TDATEOBS'] = ('2012-07-14',    'UT date of start of first exposure in file')
    hdu.header['TTIMEOBS'] = ('04:35:52',      'UT start time of first exposure in file')
    hdu.header['OBSTYPE']  = ('SPECTROSCOPIC', 'observation type - imaging or spectroscopic')
    hdu.header['OBSMODE']  = ('ACCUM',         'operating mode')
    hdu.header['SUBARRAY'] = (False,           'data from a subarray (T) or full frame (F)')
    hdu.header['DETECTOR'] = ('CCD',           'detector in use: NUV-MAMA, FUV-MAMA, or CCD')
    hdu.header['CCDAMP']   = ('D',             'CCD amplifier read out (A,B,C,D)')
    hdu.header['CCDGAIN']  = (1,               'commanded gain of CCD')
    hdu.header['CCDOFFST'] = (3,               'commanded CCD bias offset')
    hdu.header['BINAXIS1'] = (1,               'axis1 data bin size in unbinned detector pixels')
    hdu.header['BINAXIS2'] = (1,               'axis2 data bin size in unbinned detector pixels')
    hdu.writeto(test_file, output_verify='exception', clobber=True)


class TestFileFiltering:
    '''
    Tests functionality of stis_cti.viable_ccd_file
    '''
    test_file = 'testfits_2334134234_raw.fits'
    
    def setup(self):
        setup_file(self.test_file)
    
    def teardown(self):
        os.remove(self.test_file)
    
    def test_viable_file(self):
        assert viable_ccd_file(self.test_file)
    
    def test_viable_file_lenient(self):
        assert viable_ccd_file(self.test_file, \
            earliest_date_allowed = datetime.datetime(1990,1,1,0,0,0), \
            amplifiers_allowed = ['A','B','C','D'], \
            gains_allowed = [1,2,4,8], \
            offsts_allowed = range(9))
    
    def test_acq_file_reject(self):
        fits.setval(self.test_file, 'OBSMODE', value='ACQ')
        assert_false(viable_ccd_file(self.test_file))
    
    def test_gain3_file_reject(self):
        fits.setval(self.test_file, 'CCDGAIN', value=3)
        assert_false(viable_ccd_file(self.test_file))
    
    def test_gain4_file(self):
        fits.setval(self.test_file, 'CCDGAIN', value=4)
        assert viable_ccd_file(self.test_file)

    def test_pre_sm4_reject(self):
        fits.setval(self.test_file, 'TDATEOBS', value='1990-01-01')
        assert_false(viable_ccd_file(self.test_file))

    def test_binned_data_reject1(self):
        fits.setval(self.test_file, 'BINAXIS1', value=2)
        assert_false(viable_ccd_file(self.test_file))
    
    def test_binned_data_reject2(self):
        fits.setval(self.test_file, 'BINAXIS2', value=2)
        assert_false(viable_ccd_file(self.test_file))
    
    def test_subarray_reject(self):
        fits.setval(self.test_file, 'SUBARRAY', value=True)
        assert_false(viable_ccd_file(self.test_file))


class Test_determine_input_science:
    '''
    Tests functionality of stis_cti.determine_input_science
    '''
    test_dir   = 'dir_2334134234'
    test_file  = os.path.join(test_dir, 'testfits_001_raw.fits')
    test_file2 = os.path.join(test_dir, 'testfits_002_raw.fits')
    
    def setup(self):
        if not os.access(os.path.curdir, os.W_OK):
            raise Exception('Can\'t write test dir and file to CWD!')
        
        # Setup directory structure:
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
        if os.path.exists(self.test_file2):
            os.remove(self.test_file2)
        if os.path.exists(self.test_dir):
            os.rmdir(self.test_dir)
        os.mkdir(self.test_dir)
        
        # Write FITS file:
        setup_file(self.test_file)
    
    def teardown(self):
        os.remove(self.test_file)
        if os.path.exists(self.test_file2):
            os.remove(self.test_file2)
        os.rmdir(self.test_dir)
    
    def test_file_written(self):
        assert os.path.exists(self.test_file)
    
    def test_filter_pass(self):
        assert_equals(determine_input_science(self.test_dir, False, False), [self.test_file])
    
    def test_filter_date_reject(self):
        fits.setval(self.test_file, 'TDATEOBS', value='1992-01-01')
        assert_raises(FileError, determine_input_science, self.test_dir, False, False)
    
    def test_filter_date_allow_pass(self):
        fits.setval(self.test_file, 'TDATEOBS', value='1992-01-01')
        assert determine_input_science(self.test_dir, True, False)
    
    def test_filter_date_partial_reject(self):
        setup_file(self.test_file2)
        fits.setval(self.test_file2, 'TDATEOBS', value='1992-01-01')
        assert_equals(determine_input_science(self.test_dir, False, False), [self.test_file])


class Test_check_pctetab_version:
    '''
    Tests functionality of stis_cti.check_pctetab_version
    '''
    pctetab = 'a00_stis_000test37155_pcte.fits'
    
    def setup(self):
        if not os.access(os.path.curdir, os.W_OK):
            raise Exception('Can\'t write test PCTETAB file to CWD!')
        
        hdu = fits.PrimaryHDU()
        hdu.header['FILENAME'] = self.pctetab
        hdu.header['FILETYPE'] = 'PIXCTE'                                                            
        hdu.header['TELESCOP'] = 'HST'                                                            
        hdu.header['USEAFTER'] = 'Oct 01 1996 00:00:00'                                                
        hdu.header['PEDIGREE'] = 'INFLIGHT 01/10/1996 25/06/2012'                                      
        hdu.header['DESCRIP']  = 'Parameters needed for pixel-based CTE correction ------------------' 
        hdu.header['NCHGLEAK'] = (1, 'number of chg_leak extensions')
        hdu.header['INSTRUME'] = 'STIS'                                                            
        hdu.header['DETECTOR'] = 'CCD'                                                            
        hdu.header['SIM_NIT '] = (7, 'number of readout simulations done per column')
        hdu.header['SHFT_NIT'] = (4, 'the number of shifts each column readout simula')
        hdu.header['RN_CLIP']  = (5.6, 'Read noise level in electrons.')
        hdu.header['NSEMODEL'] = (1, 'Read noise smoothing algorithm.')
        hdu.header['SUBTHRSH'] = (-30.0, 'Over-subtraction correction threshold.')
        hdu.header['PCTE_VER'] = ('0.1_alpha', 'Version of PCTETAB')
        hdu.writeto(self.pctetab, output_verify='exception', clobber=True)
    
    def teardown(self):
        os.remove(self.pctetab)
     
    def test_pctetab_ver_pass(self):
        assert check_pctetab_version(self.pctetab, False, '0.1', '1.999')
    
    def test_pctetab_ver_reject(self):
        assert_raises(VersionError, check_pctetab_version, self.pctetab, False, '1.0', '1.999')
    
    def test_pctetab_ver_file_reject(self):
        fits.setval(self.pctetab, 'PCTE_VER', value='3.0_beta')
        assert_raises(VersionError, check_pctetab_version, self.pctetab, False, '1.0', '1.999')



if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['nose', '--verbosity=2'])
