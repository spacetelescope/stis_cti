#!/usr/bin/env python

from nose.tools import assert_equals, assert_false, assert_raises, assert_not_equals
import os
import datetime
from astropy.io import fits

from stis_cti import *
from archive_dark_query import archive_dark_query


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


def write_pctetab(pctetab):
    if not os.access(os.path.curdir, os.W_OK):
        raise Exception('Can\'t write test PCTETAB file to CWD!')
    
    hdu = fits.PrimaryHDU()
    hdu.header['FILENAME'] = pctetab
    hdu.header['FILETYPE'] = 'PIXCTE'                                                            
    hdu.header['TELESCOP'] = 'HST'                                                            
    hdu.header['USEAFTER'] = 'Oct 01 1996 00:00:00'                                                
    hdu.header['PEDIGREE'] = 'INFLIGHT 01/10/1996 25/06/2012'                                      
    hdu.header['DESCRIP']  = 'Parameters needed for pixel-based CTE correction ------------------' 
    hdu.header['NCHGLEAK'] = (1, 'number of chg_leak extensions')
    hdu.header['INSTRUME'] = 'STIS'                                                            
    hdu.header['DETECTOR'] = 'CCD'                                                            
    hdu.header['SIM_NIT']  = (7, 'number of readout simulations done per column')
    hdu.header['SHFT_NIT'] = (4, 'the number of shifts each column readout simula')
    hdu.header['RN_CLIP']  = (5.6, 'Read noise level in electrons.')
    hdu.header['NSEMODEL'] = (1, 'Read noise smoothing algorithm.')
    hdu.header['SUBTHRSH'] = (-30.0, 'Over-subtraction correction threshold.')
    hdu.header['PCTE_VER'] = ('0.1_alpha', 'Version of PCTETAB')
    hdu.writeto(pctetab, output_verify='exception', clobber=True)


class Test_check_pctetab_version:
    '''
    Tests functionality of stis_cti.check_pctetab_version
    '''
    pctetab = 'a00_stis_000test37155_pcte.fits'
    
    def setup(self):
        write_pctetab(self.pctetab)
    
    def teardown(self):
        os.remove(self.pctetab)
     
    def test_pctetab_ver_pass(self):
        assert check_pctetab_version(self.pctetab, False, '0.1', '1.999')
    
    def test_pctetab_ver_reject(self):
        assert_raises(VersionError, check_pctetab_version, self.pctetab, False, '1.0', '1.999')
    
    def test_pctetab_ver_file_reject(self):
        fits.setval(self.pctetab, 'PCTE_VER', value='3.0_beta')
        assert_raises(VersionError, check_pctetab_version, self.pctetab, False, '0.1', '1.999')


def write_superdark(superdark):
    if not os.access(os.path.curdir, os.W_OK):
        raise Exception('Can\'t write test superdark file to CWD!')
    
    hdu = fits.PrimaryHDU()
    
    hdu.header['FILENAME'] = superdark
    hdu.header['FILETYPE'] = 'DARK IMAGE'
    hdu.header['TELESCOP'] = 'HST'
    hdu.header['REF_TEMP'] = 18
    hdu.header['DRK_VS_T'] = 0.07
    hdu.header['PCTECORR'] = 'COMPLETE'
    hdu.header['PCTETAB']  = 'pctetab$a01_stis_pcte.fits'
    hdu.header['USEAFTER'] = 'Jul 24 2012 00:59:27'
    hdu.header['PEDIGREE'] = 'INFLIGHT 01/10/1996 25/06/2012'
    hdu.header['DESCRIP']  = 'Weekly gain=1 dark for STIS CCD data taken after Jul 24 2012-------'
    hdu.header['INSTRUME'] = 'STIS'
    hdu.header['DETECTOR'] = ('CCD', 'detector in use: CCD')
    
    hdu.header['PCTEFRAC'] = (1.192887185841836,              'CTE time scaling value')
    hdu.header['PCTESMIT'] = (7,                              'PCTE readout simulation iterations')
    hdu.header['PCTESHFT'] = (4,                              'PCTE readout number of shifts')
    hdu.header['PCTERNCL'] = (5.6,                            'PCTE readnoise amplitude')
    hdu.header['PCTENSMD'] = (1,                              'PCTE readnoise mitigation algorithm')
    hdu.header['PCTETRSH'] = (-30.0,                          'PCTE over-subtraction threshold')
    hdu.header['PCTE_VER'] = ('0.1_alpha',                    'Version of PCTETAB')
    hdu.header['CTE_NAME'] = ('PixelCTE 2012',                'name of CTE algorithm')
    hdu.header['CTE_VER']  = ('3.2',                          'version of CTE algorithm')
    hdu.header['BASEDARK'] = ('$ctitest/bdark_d133_drk.fits', 'Used to make weekdark')
    
    hdu.header['HISTORY'] = 'blah, blah blah'
    hdu.header['HISTORY'] = 'The following input files were used:'
    hdu.header['HISTORY'] = 'abc_cte.fits'
    hdu.header['HISTORY'] = 'def_cte.fits'
    hdu.header['HISTORY'] = 'hij_cte.fits'
    hdu.header['HISTORY'] = ''
    hdu.header['HISTORY'] = 'blah2, blah2, blah2'
    
    hdu.writeto(superdark, output_verify='exception', clobber=True)


class Test_superdark_hash:
    '''
    Tests functionality of stis_cti.superdark_hash
    '''
    pctetab = 'a00_stis_000test37166_pcte.fits'
    superdark = 'd001_testfile57254_drk.fits'
    
    def setup(self):
        write_superdark(self.superdark)
        write_pctetab(self.pctetab)
    
    def teardown(self):
        os.remove(self.superdark)
        os.remove(self.pctetab)
    
    def test_manual_inputs_pass(self):
        assert_equals(
            superdark_hash(sim_nit=7, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.1_alpha', files=[]), 
            superdark_hash(pctetab=self.pctetab, files=[]))
    
    def test_manual_inputs2_pass(self):
       # superdark_hash should ignore text after '_' in pcte_ver
        assert_equals(
            superdark_hash(sim_nit=7, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.1_beta', files=[]), 
            superdark_hash(pctetab=self.pctetab, files=[]))
    
    def test_manual_inputs3_pass(self):
        # try with a file list specified
        assert_equals(
            superdark_hash(sim_nit=7, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.1_alpha', files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    
    def test_manual_inputs4_pass(self):
        assert_equals(
            superdark_hash(sim_nit=7, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.1_alpha', files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']), 
            superdark_hash(superdark=self.superdark))
    
    def test_manual_inputs_reject(self):
        # changed file list
        assert_not_equals(
            superdark_hash(sim_nit=7, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.1_alpha', files=[]), 
            superdark_hash(pctetab=self.pctetab, files=['otherfile_cte.fits']))
    
    def test_manual_inputs2_reject(self):
        # changed sim_nit
        assert_not_equals(
            superdark_hash(sim_nit=70, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.1_alpha', files=[]), 
            superdark_hash(pctetab=self.pctetab, files=[]))
    
    def test_manual_inputs3_reject(self):
        # changed pcte_ver
        assert_not_equals(
            superdark_hash(sim_nit=7, shft_nit=4, rn_clip=5.6, nsemodel=1, subthrsh=-30.0, 
                           pcte_ver='0.2_alpha', files=[]), 
            superdark_hash(pctetab=self.pctetab, files=[]))
    
    def test_superdark_pass(self):
        assert_equals(
            superdark_hash(superdark=self.superdark, files=[]), 
            superdark_hash(pctetab=self.pctetab, files=[]))
    
    def test_superdark2_pass(self):
        assert_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    
    def test_superdark3_pass(self):
        fits.setval(self.superdark, 'PCTE_VER', value='0.1_gamma')
        assert_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    
    def test_superdark_reject(self):
        assert_not_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=[]))
    
    def test_superdark2_reject(self):
        fits.setval(self.superdark, 'PCTENSMD', value=5)
        assert_not_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    
    def test_superdark3_reject(self):
        fits.setval(self.pctetab, 'nsemodel', value=6)
        assert_not_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    
    def test_superdark4_reject(self):
        fits.setval(self.superdark, 'PCTE_VER', value='0.2_alpha')
        assert_not_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    
    def test_superdark4_reject(self):
        fits.setval(self.pctetab, 'PCTE_VER', value='0.2_alpha')
        assert_not_equals(
            superdark_hash(superdark=self.superdark), 
            superdark_hash(pctetab=self.pctetab, files=['abc_cte.fits', 'def_cte.fits', 'hij_cte.fits']))
    

class Test_archive_dark_query:
    '''
    Tests functionality of stis_cti.archive_dark_query
    '''
    test_file = 'testfits_2334134667_raw.fits'
    
    @classmethod
    def setup_class(cls):
        setup_file(cls.test_file)
        cls.anneal = archive_dark_query([cls.test_file], None, None, False, False)
    
    @classmethod
    def teardown_class(cls):
        os.remove(cls.test_file)
    
    def test_number_of_anneals(self):
        assert_equals(len(self.anneal), 1)
    
    def test_anneal_index(self):
        assert_equals(self.anneal[0]['index'], 133)
    
    def test_anneal_darks(self):
        darks = [d['exposure'] for d in self.anneal[0]['darks']]
        assert_equals(set(darks), \
            set(['OBVM3XH9Q', 'OBVM3YHHQ', 'OBVM3ZN5Q', 'OBVM40NCQ', 'OBVM41T2Q', 'OBVM42TFQ', 
                 'OBVM43YEQ', 'OBVM44YOQ', 'OBVM45FQQ', 'OBVM46G7Q', 'OBVM47LNQ', 'OBVM48MCQ', 
                 'OBVM49ANQ', 'OBVM4AAVQ', 'OBVM4BH5Q', 'OBVM4CHJQ', 'OBVM4DLRQ', 'OBVM4EM2Q', 
                 'OBVM4FRJQ', 'OBVM4GRQQ', 'OBVM4HYVQ', 'OBVM4IZ5Q', 'OBVM4JG4Q', 'OBVM4KGCQ', 
                 'OBVM4LO1S', 'OBVM4MOIS', 'OBVM4NABQ', 'OBVM4OAMQ', 'OBVM4PI2Q', 'OBVM4QIAQ', 
                 'OBVM4RPDQ', 'OBVM4SPJQ', 'OBVM4TW9Q', 'OBVM4UWDQ', 'OBVM4VCOQ', 'OBVM4WCXQ', 
                 'OBVM4XJ6Q', 'OBVM4YJCQ', 'OBVM4ZOMQ', 'OBVM50POQ', 'OBVM51A4Q', 'OBVM52ABQ', 
                 'OBVM53FIQ', 'OBVM54FNQ', 'OBVM55LTQ', 'OBVM56M5Q', 'OBVM57T3S', 'OBVM58T8S', 
                 'OBVM59W6S', 'OBVM5AWDS', 'OBVM5BBQQ', 'OBVM5CBUQ', 'OBVM5DG7Q', 'OBVM5EGBQ', 
                 'OBVM5FALQ', 'OBVM5GAUQ', 'OBVM5HFTQ', 'OBVM5IFPQ', 'OBVM5JKLQ', 'OBVM5KLPQ']))
    
    def test_undefined_fits_file(self):
        undefined_filename = 'testfits_undefined_28739723_raw.fits'
        assert_raises(IOError, archive_dark_query, [undefined_filename], None, None, False, False)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['nose', '--verbosity=2'])
