.. _sec-getting_started:

===============
Getting Started
===============

ONETOOL is freely available. However, we ask you to add an appropriate statement (including the NIH grant number) 
under "acknowledgments" in any publication of results obtained by using this program.

Suggested wording is:

`"(Some of)The results of this paper were obtained by using the software package ONETOOL, which was 
supported by the National Research Foundation of Korea Grant funded by Korean Government (NRF-2014S1A2A2028559)."`


License:

- Copyright 2017 ONETOOL team

----------------------
Getting ONETOOL binary
----------------------

The pre-built ONETOOL binaries can be downloaded from the `homepage <http://healthstat.snu.ac.kr/software/onetool/pages/download.php>`_ for free.

ONETOOL has been compiled and tested on the platforms indicated below.

================   =====================   ==========================
Platform           Operating system type   Operating system version
================   =====================   ==========================
64-bit (x86_64)    Linux                   Kernel 2.6+ RHEL5, Ubuntu
64-bit (x86_64)    Windows                 Vista, 7, 10
================   =====================   ==========================

Two different versions are provided, static and dynamic, for each platform.

- **Static** version does not include R library, so will run under most system environments, but pedigree plot option (``--plot``) is not available.
- **Dynamic** version requires R system library in the right place when it is executed.  This version supports pedigree plot function (``--plot``).

Depending on the version, you will need:

- Intel Fortran Compiler Redistributable is required to run Windows version of ONETOOL.
- R package kinship2 is required to use ``--plot`` function in ONETOOL for both dynamic versions.
- For Windows, ONETOOL expect R is installed on the initial installation directory (``C:\Program Files\R\R-X.X.X``).

Unzip the distribution after you download it.

For Linix version, **un-tar** ONETOOL package.

.. code-block:: text

    tar -xzf ONETOOL.v2.0_Linux_x86_64_dynamic_withR.tgz

For Windows version, double-click to **unzip** ONETOOL distribution.

Once you've set up ONETOOL, we recommend that you download the `example files <http://healthstat.snu.ac.kr/software/onetool/pages/download.php>`_. 
for some test runs.  Most of examples in this documentation are using the same example files.

----------------------------
Building ONETOOL from source
----------------------------

The ONETOOL source code.  To clone the `ONETOOL repository <https://github.com/SWonlab/onetool>`_ using `Git <https://git-scm.com/>`_, run

  .. code-block:: text

      $ git clone https://github.com/SWonlab/onetool.git
      $ cd onetool

  You can also download the source code directly from `Github <https://github.com/SWonlab/onetool/archive/master.zip>`_.

ONETOOL uses several external libraries.  These should load automatically.  If not, these must be explicitly installed.

 - `boost <http://www.boost.org/>`_

 - `blas and lapack <http://www.netlib.org/lapack/>`_

 - `R 3.3.1 <http://www.r-project.org/>`_ with packages ``kinship``.

---------------
Running ONETOOL
---------------

To run ONETTOL binary, go to the directory where ONETOOL is installed.

.. code-block:: text

    $ ./onetool --fam test_miss0.fam

Several ONETOOL tests are available on `homepage <http://healthstat.snu.ac.kr/software/onetool/pages/tutorial_new.php>`_.
