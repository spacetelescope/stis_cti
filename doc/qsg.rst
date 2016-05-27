=================
Quick Start Guide
=================

This section is intended to provide a dialogue and inputs which will give 
good results for beginning users.  Please let us know how it goes, and ask 
any questions you have.  Send feedback and questions to help@stsci.edu.

Please note that AURA/STScI is not responsible for any damages resulting 
from use of this software.

------------------------------------------------------------------------------------------

As you are working through this guide, you can copy and paste text 
from this file into your terminal window, if that helps you work 
faster.  Or you may want to copy blocks of text into a text file, 
to create your own customized dialog and commands.

If you don't have it already, you will need to install the AstroConda distribution with 
legacy IRAF support available at 
http://astroconda.readthedocs.io.

After AstroConda is installed, open a new session (open a new Bash terminal and type 
``source activate astroconda``).

To install the ``stis_cti`` package and dependencies, run:

::
  
  pip install stis_cti

Or, to upgrade from a previous version, run:

::
  
  pip install --upgrade --no-deps stis_cti
  pip install stis_cti
  pip install --upgrade --no-deps refstis
  pip install refstis

After that, open a brand new AstroConda session (see above), and continue there.  
Create a working directory where the data will reside, and also create the 
necessary sub-directories.  For example, our disk here is called ``/Users``, 
and the account is ``biretta``, and we might use the proposal number ``13542``, 
so our working directory was ``/Users/biretta/13542``.

::
  
  source activate astroconda
  mkdir /Users/biretta/13542
  cd    /Users/biretta/13542
  mkdir darks
  mkdir ref
  mkdir science

Place the raw science data in the ``science`` directory.  This will
include the files you got from the archive:
``*_raw.fits``, ``*_epc.fits``, ``*_spt.fits``, ``*_asn.fits``, ``*_wav.fits``

Note that the presense of non-CTI-corrected products from the archive 
(e.g. ``*_flt.fits``) will not interfere with the processing.

We will assume you want the script to automatically download the 
reference files from the HST archive.  We do an initial run of the 
script from the science directory:

::
  
  stis_cti --crds_update

It is expected the script will run for a minute, and then give an
error message, along with a URL to get the dark frames that are
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

Next we need to get the missing dark frames.  Copy the entire URL which 
the script generated, and paste it into a web browser (e.g. on a Mac 
highlight the URL with the cursor, and hit ``command-c``, move the cursor 
to the web browser URL bar, and hit ``command-v``).  Then hit ``return``.  
This will generate an HST archive request for the missing dark files.  
Then do this on the archive web page:

Click ``[Mark all]``

Click ``[Submit marked data for retrieval from STSDAS]``

This will bring up a new page.  Fill out your HST archive credentials,
and the computer and directory where you want the files to be sent:

Check ``[sftp the data]``

::
  
  myserver.mydomain.edu (put your computer name here instead)

  /Users/biretta/13542/darks  (use the darks directory that we created earlier)

  biretta (put your computer account name)

  (put your computer password)

Alternatively, you may stage and retrieve files from the archive's ftp server.

File options: check ``Calibrated`` (it will probably checked already by default)

Click ``[Send retrieval request to ST-DADS]``

Then wait for the dark frames to be delivered by the HST archive....
After you receive an email from archive.stsci.edu that the request has
completed successfully, put the darks in the ``darks/`` directory and run 
``stis_cti`` again.  This time it should run to completion.  On a typical Mac 
laptop, it might take an hour to run.  Make sure you are in the science 
directory still, and then:

::
  
  stis_cti --crds_update

When you get the message that looks like:

.. parsed-literal:: 
   
   Completion time:                2015-06-12 19:51:29.862291
   Run time:                       0:23:05.068940
   stis_cti.py complete!


...it is done running.  You should find the output files in the ``science/`` 
directory with names like ``*_cte.fits``, ``*_flc.fits``, ``*_crc.fits``, etc.

.. Warning::
   
   For recent STIS observations (new data taken in last 30 to 60 days) optimal dark 
   reference files will not yet be available.  This will affect the selection of data 
   being used to generate the CTI-corrected super-darks.  To get the most accurate 
   calibration, please re-reduce your data after the pipeline's new super-biases and 
   super-darks have been delivered by deleting the relevant old CTI-corrected super-darks 
   in the ``ref/`` directory and running ``stis_cti`` with the ``--clean`` and 
   ``--crds_update`` options specified.  You may need to download additional component 
   darks from MAST.
   
   To receive updates when STIS reference files are delivered to CRDS, go to 
   https://maillist.stsci.edu and subscribe to the ``stis_reffiles_upd`` mailing list.
   
   You can also check the status of super-dark and super-bias files by going to 
   https://hst-crds.stsci.edu and clicking on STIS-->darkfile and STIS-->biasfile.  Sort 
   by USEAFTER to see if the week corresponding to your science data has been delivered 
   yet.

Good luck!  Let us know if you encounter problems, or need any assistance at 
help@stsci.edu.
