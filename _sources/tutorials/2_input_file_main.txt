.. _sec-ifile-main:

====================
Main Input Data File
====================

.. _ifam:

FAM file
========

This is the FAM file in `PLINK`_.

Each lines of the FAM file describes an individual.  It contains a
white-space (space or tab) delimited records for each individual including
fields for identifier, sex, two parents and an optional trait field.

    #. Family ID
    #. Individual ID
    #. Father ID
    #. Mother ID
    #. Sex (1=male; 2=female; other=unknown)
    #. Phenotype (optional)

.. note:: Unlike in PLINK, the last column with dummy data can be omitted in the FAM file when the main phenotype is in the phenotype file.
.. note:: It is the required input file and can be used alone for certain types of analyses in ONETOOL.

    .. code-block:: text

        $ onetool --fam test_miss00.fam


.. _ipheno:

PHENO file
==========

This is the alternate phenotype file in `PLINK`_ to specify an alternate phenotype for analysis, 
i.e. other than the one in the PED (or FAM) file.

Each lines of the PHENO file contains the phenotype data for an individual. 
The first two columns contain the identifiers and the rest of columns are
phenotype data.

    #. Family ID
    #. Individual ID
    #. Phenotype1
    #. Phenotype2
    #. Phenotype3
    #. Phenotype4
    #. ...

The detailed description of the PHENO file can be found in `PLINK`_.

.. code-block:: text

    $ onetool --fam test_miss00.fam --vcf test_miss00.vcf --pheno test_miss00.pheno --pname sbp


.. _ivcf:

VCF file
========

This is the Variant Call Format(VCF) file used in `1000 Genomes Project`_.
Files in both plain text format (.vcf) or gzipped format (.bcf) are supported.
The meta information lines (starting with ##) are ignored.

The first 9 columns in header and data lines are:

    #. CHROM
    #. POS
    #. ID
    #. REF
    #. ALT
    #. QUAL
    #. FILTER
    #. INFO
    #. FORMAT

.. note:: The sample IDs in the header line (starting with #CHROM) have to match with the individual IDs in FAM file uniquely.

.. code-block:: text

    $ onetool --fam test_miss00.fam --vcf test_miss00.vcf


.. _ibed_bim:

BED/BIM file
============

BED file is the PLINK binary PED file and it is used together with BIM file
(extended MAP file: two extra cols = allele names).

The detailed description of the BED/BIM files can be found in `PLINK`_.

.. code-block:: text

    $ onetool --fam test_miss00.fam --bed test_miss00.bed --bim test_miss0.bim


.. _iimpute2:

IMPUTE2 file
============

IMPUTE2 file is the output files from IMPUTE2 program and both .impute2 and .impute2_info files are required.

The genotype file stores (.impute2) data on a one-line-per-SNP format. 
The first 5 entries of each line should be:

    #. SNP ID
    #. RS ID of the SNP
    #. base-pair position of the SNP
    #. the allele coded A
    #. the allele coded B.

The next three numbers on the line should be the probabilities of the three genotypes AA, AB and BB 
at the SNP for the first individual in the sample and the next three numbers for the second individual and so on.

The SNP-wise information file (.impute2_info) contains the following columns (header shown in parentheses): 

    #. SNP identifier from -g file (snp_id) 
    #. rsID (rs_id) 
    #. base pair position (position) 
    #. expected frequency of allele coded '1' in the -o file (exp_freq_a1) 
    #. measure of the observed statistical information associated with the allele frequency estimate (info) 
    #. average certainty of best-guess genotypes (certainty) 
    #. internal "type" assigned to SNP (type)

The detailed description of the IMPUTE2 files can be found in `IMPUTE2`_.

.. code-block:: text

    $ onetool --fam test_miss00.fam --dosage test_miss00.impute2 --mqls


.. _iped_map:

PED/MAP file
============

PED file is the default a white-space (space or tab) delimited text file format 
used in `PLINK`_.  The same first six columns as in FAM file are mandatory.

When a PED file is used, it has to be used with PLINK MAP file which
contains the following 4 columns.

    #. chromosome (1-22, X, Y or 0 if unplaced)
    #. rs# or snp identifier
    #. Genetic distance (morgans)
    #. Base-pair position (bp units)

The detailed description of the PED/MAP files can be found in `PLINK`_.

.. code-block:: text

    $ onetool --ped test_miss00.ped --map test_miss00.map


.. _PLINK: http://zzz.bwh.harvard.edu/plink/
.. _1000 Genomes Project: http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
.. _MERLIN: http://csg.sph.umich.edu/abecasis/merlin/tour/input_files.html
.. _S.A.G.E.: http://darwin.cwru.edu/sage/pages/download.php
.. _IMPUTE2: http://mathgen.stats.ox.ac.uk/impute/impute_v2.html

