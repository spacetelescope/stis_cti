=================
Quick Start Guide
=================

Thank you very much for trying out this **beta** release of the
``stis_cti`` correction package.  This section is intended to provide a dialogue and inputs which will give good results for beginning users. 
Please let us know how it goes, and ask any questions you have.  Send
feedback, questions, etc. to biretta@stsci.edu.
Please note: we are not responsible for any damages resulting from
use of this experimental software.

Note: as you are working through this guide, you can copy and paste
text from this PDF file into your terminal window, if that helps you
work faster.  Or you may want to copy blocks of text from this PDF into a textfile,
to create your own customized dialog and commands.

If you don't have it already, you will need to install the standard
Ureka/STSDAS distribution available at
http://ssb.stsci.edu/ureka/. The normal version (not SSBX nor SSBDEV)
is preferred for our scripts.

After Ureka is installed, open a new Ureka session (either click on the
Ureka icon on the desktop, or open a new terminal window and type either ``ur_setup common ssbrel``
or ``ssbrel``, depending on your installation details).

Next we need to install the ``stis_cti`` package. Create a directory
where the scripts will be installed, for example:

::

  mkdir ~/cti
  cd ~/cti

::

Put the tar files we emailed to you in this directory, then unpack and install the scripts:

::

  pip install crds
  tar -xzf refstis-0.0.1.tar.gz
  cd refstis-0.0.1
  ./setup.py install
  cd ..
  tar -xzf stis_cti-0.4-beta2.tar.gz
  cd stis_cti-0.4-beta2
  ./setup.py install
  cd ..

::

After that open a brand new Ureka session (see above), and continue there.  Create a working directory where the data will reside, and also create the necessary sub-directories.  For example, our disk here is called "/Users" and the account is "biretta" and we might use the proposal number "13542" so our working directory was "/Users/biretta/13542".

::

  ur_setup common ssbrel
  mkdir /Users/biretta/13542
  cd    /Users/biretta/13542
  mkdir darks
  mkdir ref
  mkdir science
  mkdir oref

::

Place the raw science data in the "science" directory.  This will
include the files you got from the archive: \*_raw.fits, \*_epc.fits, \*_spt.fits, \*_asn.fits, \*_trl.fits

Let's move any normal pipeline calibrated products out of the way to a sub-directory "pipeline":

::
  
  cd science
  mkdir pipeline
  mv *crj.fits pipeline
  mv *flt.fits pipeline
  mv *x1d.fits pipeline
  mv *x2d.fits pipeline
  mv *sx1.fits pipeline
  mv *sx2.fits pipeline

::

Next we define the location of the calibration reference
files. We will assume you want the script to automatically download the reference files from the HST
archive:  

::

  setenv oref "/Users/biretta/13542/ref/references/hst/stis/"

::

Next we do an initial run of the script from the science directory:

::

  stis_cti --crds_update

::

It is expected the script will run for a minute, and then give an
error message, along with a web url to get the dark frames that are
missing, and then stop:

::

   ERROR:  These FLT component darks are missing from ../darks/:
   .
   .
   .
   Please download the missing darks (calibrated FLTs) via this link:
  
   http://archive.stsci.edu/hst/search.php?sci_instrume=STIS&sci_instrument_config=
   STIS%2FCCD&sci_targname=DARK&sci_aec=C&resolve=don%27tresolve&sci_data_set_name=
   OC4W6XH3Q%2COC4W6YHBQ%2COC4W6ZP2Q%2COC4W70PCQ%2COC4W71TEQ%2COC4W72TOQ%2COC4W73X8Q%
   2COC4W74XJQ%2COC4W75D0Q%2COC4W76DCQ%2COC4W77HHQ%2COC4W78I0Q%2COC4W79A5Q%2COC4W7AADQ%
   2COC4W7BFGQ%2COC4W7CF9Q%2COC4W7DJNQ%2COC4W7EJRQ%2COC4W7FOAQ%2COC4W7GO4Q%2COC4W7HSNQ%
   2COC4W7ISUQ%2COC4W7JXEQ%2COC4W7KXAQ%2COC4W7LGRQ%2COC4W7MGWQ%2COC4W7NA1Q%2COC4W7OA8Q%
   2COC4W7PM6Q%2COC4W7QMDQ%2COC4W7RTJQ%2COC4W7STNQ%2COC4W7TX4Q%2COC4W7UXDQ%2COC4W7VIKQ%
   2COC4W7WIRQ%2COC4W7XNJQ%2COC4W7YNRQ%2COC4W7ZSZQ%2COC4W80TMQ%2COC4W81A4Q%2COC4W82AGQ%
   2COC4W83NMQ%2COC4W84O1Q%2COC4W85SRQ%2COC4W86SZQ%2COC4W87XWQ%2COC4W88YHQ%2COC4W89D6Q%
   2COC4W8ADJQ%2COC4W8BHWQ%2COC4W8CI2Q%2COC4W8DNUQ%2COC4W8EOAQ%2COC4W8FBPQ%2COC4W8GBTQ&
   max_records=50000&max_rpp=5000&ordercolumn1=sci_start_time&action=Search

   .
   .
   .
   stis_cti.stis_cti.FileError: Missing component dark FLT files.

::

Next we need to get the missing dark frames.  Copy the entire URL which the script generated, and paste it into a web browser (e.g. on a
Mac highlight the URL with the cursor, and hit "command-C", move the cursor
to the web browser URL window, and hit "command-V").  Then hit return.
This will generate a HST archive request for the missing dark files.  Then do this on the archive web page:

Click ``[Mark all]``

Click ``[Submit marked data for retrieval from STSDAS]``

This will bring up a new page.  Fill out your HST archive credentials,
and the computer and directory where you want the files to be sent:

Check ``[sftp the data]``

::

  plhstins2.stsci.edu (put your computer name here instead)

  /Users/biretta/13542/darks  (use the darks directory that we created earlier)

  biretta (put your computer account name)

  (put your computer password)

::

File options: check ``Calibrated`` (it will probably checked already by default)

Click ``[Send retrieval request to ST-DADS]``

Then wait for the dark frames to be delivered by the HST archive....
After you receive an email from archive.stsci.edu that the request has
completed successfully, run ``stis_cti`` again.  This time it should run
to completion.  On a typical Mac laptop, it might take an hour to
run.  Make sure you are in the science directory still, and then:

::

  stis_cti --crds_update

::

When you get the message that looks like:

.. parsed-literal:: 

   Completion time:                2015-06-12 19:51:29.862291
   Run time:                       0:23:05.068940
   stis_cti.py complete!


...it is done running.  You should find the output files in the science directory with names like \*_cte.fits, \*_flc.fits, \*_crc.fits, etc.

Good luck!  Let us know if you encounter problems, or need any help.

-- John Biretta  (biretta@stsci.edu)




