def custom_superdark_info():
    '''
    Prints informative text on manually applying the CTI-correction and creating custom
    super-darks.
    '''
    import stis_cti
    import pkg_resources
    import glob
    import os
    
    # Find the stis_cti package's most recent PCTETAB:
    pctetabs = glob.glob(pkg_resources.resource_filename(stis_cti.__name__, 'data/*_pcte.fits'))
    pctetabs.sort()
    if len(pctetabs) > 0:
        pctetab = pctetabs[-1]
    pctetab_dir = os.path.dirname(pctetab) + os.path.sep
    pctetab_file = os.path.basename(pctetab)
    
    print 'How to create and apply custom super-darks'
    print '------------------------------------------\n'
    print 'See the documentation at:  http://pythonhosted.org/stis_cti/readme.html#custom-super-darks\n'
    
    print 'The current stis_cti package PCTETAB reference file is:\n    {}\n'.format(pctetab)
    
    print 'To set the $pctetab environmental variable in the shell:'
    print '    In TCSH:  setenv pctetab {}'.format(pctetab_dir)
    print '    In BASH:  export pctetab="{}"'.format(pctetab_dir)
    print "    Within Python:"
    print "        import os"
    print "        os.environ['pctetab'] = '{}'\n".format(pctetab_dir)
    
    print 'In Python, set the PCTETAB header keyword in the appropriate files to be corrected:'
    print '    from astropy.io import fits\n'
    print "    fits.setval('darks/filename_flt.fits', 'DARKFILE', value='pctetab${}')".format(pctetab_file)
    print '    -- or --'
    print "    fits.setval('science/filename_raw.fits', 'DARKFILE', value='pctetab${}')\n".format(pctetab_file)
    
    print 'Now you can run stis_cti.StisPixCteCorr.CteCorr() on the file to produce a _cte.fits intermediate file:'
    print '    import stis_cti'
    print "    stis_cti.StisPixCteCorr.CteCorr('darks/filename_flt.fits')\n"
    
    print '(Of course, you probably want to loop over all files of interest.)\n'
    
    print 'Component darks corrected this way can be fed into refstis.basedark.make_basedark() and '
    print 'refstis.weekdark.make_weekdark() in order to make custom super-darks.\n'
    print 'First, manually separate the darks within the anneal month of interest and the week'
    print 'surrounding the science files.  Then, run:'
    
    print "    import refstis"
    print "    import glob\n"
    print "    month_files = glob.glob('annealing_month/*_cte.fits')"
    print "        # Assuming only the annealing month's darks are selected"
    print "    refstis.basedark.make_basedark(month_files, refdark_name='basedark_drk.fits')\n"
    print "    week_files = glob.glob('my_week/*_cte.fits')  # Where we have moved the week's files to my_week/"
    print "    refstis.weekdark.make_weekdark(week_files, refdark_name='weekdark_drk.fits',"
    print "        thebasedark='basedark_drk.fits')\n"
    
    print 'Be sure to set the header keyword PCTECORR to COMPLETE in super-darks to use them with stis_cti:'
    print "    fits.setval('weekdark_drk.fits', 'PCTECORR', value='COMPLETE')\n"
    
    print 'Point each science file to the new super-dark:'
    print "    fits.setval('science/filename_raw.fits', 'DARKFILE', value='stisref$weekdark_drk.fits')\n"
    
    print 'Also, don\'t run stis_cti with the --crds_update option, as this will override your custom'
    print 'super-dark.'
