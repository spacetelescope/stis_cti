=================
Quick Start Guide
=================

This section is intended to provide a dialogue and inputs which will give 
good results for beginning users.  Please let us know how it goes, and ask 
any questions you have.  Send feedback and questions to https://hsthelp.stsci.edu.

Please note that AURA/STScI is not responsible for any damages resulting 
from use of this software.

------------------------------------------------------------------------------------------

As you are working through this guide, you can copy and paste text 
from this file into your terminal window, if that helps you work 
faster.  Or you may want to copy blocks of text into a text file, 
to create your own customized dialog and commands.

If you don't have it already, you will need to install the ``stenv`` conda distribution:
https://stenv.readthedocs.io

After ``stenv`` is installed, open a new session (open a new Bash terminal and type 
``conda activate stenv``).

To install the ``stis_cti`` package and dependencies, run:

::
  
  pip install stis_cti

Or, to upgrade from a previous version, run:

::
  
  pip install --upgrade --no-deps stis_cti
  pip install stis_cti
  pip install --upgrade --no-deps refstis
  pip install refstis

Create a working directory where the data will reside, and also create the 
necessary sub-directories.  For example, our disk here is called ``/Users``, 
and the account is ``biretta``, and we might use the proposal number ``13542``, 
so our working directory was ``/Users/biretta/13542``.

::
  
  conda activate stenv
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
   Please download missing darks (calibrated FLTs) from MAST
   https://mast.stsci.edu/search/ui/#/hst
   (or specify the proper dark_dir [darks/])
   
   If missing files are expected (e.g. amp=A darks in MAST), then run with
   --ignore_missing flag.
   
   OC4W6XH3Q, OC4W6YHBQ, OC4W6ZP2Q, OC4W70PCQ, OC4W71TEQ, OC4W72TOQ, OC4W73X8Q, OC4W74XJQ,
   OC4W75D0Q, OC4W76DCQ, OC4W77HHQ, OC4W78I0Q, OC4W79A5Q, OC4W7AADQ, OC4W7BFGQ, OC4W7CF9Q,
   OC4W7DJNQ, OC4W7EJRQ, OC4W7FOAQ, OC4W7GO4Q, OC4W7HSNQ, OC4W7ISUQ, OC4W7JXEQ, OC4W7KXAQ,
   OC4W7LGRQ, OC4W7MGWQ, OC4W7NA1Q, OC4W7OA8Q, OC4W7PM6Q, OC4W7QMDQ, OC4W7RTJQ, OC4W7STNQ,
   OC4W7TX4Q, OC4W7UXDQ, OC4W7VIKQ, OC4W7WIRQ, OC4W7XNJQ, OC4W7YNRQ, OC4W7ZSZQ, OC4W80TMQ,
   OC4W81A4Q, OC4W82AGQ, OC4W83NMQ, OC4W84O1Q, OC4W85SRQ, OC4W86SZQ, OC4W87XWQ, OC4W88YHQ,
   OC4W89D6Q, OC4W8ADJQ, OC4W8BHWQ, OC4W8CI2Q, OC4W8DNUQ, OC4W8EOAQ, OC4W8FBPQ, OC4W8GBTQ
   
Next we need to get the missing dark frames.  Copy the filenames which 
the script generated, and paste it into the "Dataset Name=" field on
the MAST search page.
Select "Observations" = "All"

Click the box in the upper-left to select all rows

Click "Download Selected" --> "Choose which files to download"

Under "File category" select "Calibrated" --> "FLT"

Click the "Start Download" button.  A typical anneal month will be ~600 MB.

Once the dark FLT files are downloaded, locate them in your download directory,
unzip the archive, and move the files into the ``darks/`` directory and run 
``stis_cti`` again.  This time it should run to completion.  On a typical Mac 
laptop, it might take 15 minutes to run.  Make sure you are in the science 
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
https://hsthelp.stsci.edu.

Continue to `the full documentation <readme.html>`_...
