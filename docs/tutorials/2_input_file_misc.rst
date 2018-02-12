.. _sec-ifile-misc:

=====================
Additinal Input files
=====================

.. _ipar:

PAR file
========

This is S.A.G.E. parameter file which is an optional file to run a
segregation analysis in ONETOOL.  It is a text file containing a list of
S.A.G.E. SEGREG instructions written according to a specific syntax (i.e.,
instruction blocks). 

For the detailed description of the PAR file, see the `S.A.G.E.`_ user
manual.

.. note:: When a ``par`` file is not specified, ONETOOL runs the default analysis in SEGREG which is the comingling analysis.

.. code-block:: text

    $ onetool --fam test_miss00.fam --pheno test_miss00.pheno --pname t2d --segreg --par segreg.par


.. _ityp:

TYP file
========

This is S.A.G.E. trait genotype probability file which is required to run a
model-based linkage analysis in ONETOOL.  It is a text file containing the
individual specific type probabilities conditional on the model and all
pedigree information available, and individual specific penetrance
information.  This file is prodiced by a segregation analysis.

For the detailed description of the PAR file, see the `S.A.G.E.`_ user
manual.

.. code-block:: text

    $ onetool --fam test_miss00.fam --vcf test_miss00.vcf -lodlink --typ test_segreg.typ


.. _imap:

MERLIN MAP file
===============

This is MERLIN genetic map file which is required to run a model-free
linkage analysis in ONETOOL.

For the detailed description of this MAP file, go to `MERLIN`_.

.. note:: The PLINK MAP file contains 4 columns while the MERLIN MAP file contains 3 columns (sex-average) or 5 columns (sex-specific).  The PLINK MAP file with the correct genetic map information in 3rd column can be used for the model-free linkage analysis in ONETOOL.

.. code-block:: text

    $ onetool --fam test_miss00.fam --vcf test_miss00.vcf --merlin --map test_miss00.map


.. _iset:

SET file
========

This is a required file to run gene-based association analyses for rare
variants.  It contains the SNP clustering information in genes in the
following format:

    #. Gene ID
    #. SNP ID

.. code-block:: text

    $ onetool --fam test_miss00.fam --vcf test_miss00.vcf --genetest --set test_gene.txt


.. _MERLIN: http://csg.sph.umich.edu/abecasis/merlin/tour/input_files.html
.. _S.A.G.E.: http://darwin.cwru.edu/sage/pages/download.php

