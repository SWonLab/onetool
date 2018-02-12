.. _sec-input:


ONETOOL supports many different types of data in many different file formats.
For a run, however, it only requires the input file(s) that are relevant to
the specified analyses in that run.

================
Input Data Types
================

Though the main purpose of ONETOOL is for family-based big data analyses, it
can analyze unrelated individual data as well.  

Types of data supported:

-  **Sample**

    #. Family samples (related individuals)
    #. Independent samples (unrelated individuals)

-  **Phenotype**

    #. Binary
    #. Continuous

-  **Variant**

    #. genotype - common and rare SNPs
    #. genotype probability/dosage - imputed variants

.. note:: The terms 'sample', 'subject' and 'individual' are used interchangeably.
.. note:: The terms 'variant', 'SNP' and 'genotype' are used interchangeably.
.. note:: The terms 'trait' and 'phenotype' are used interchangeably.


================
Input File Types
================

The types of input file (with the expected extension in parenthesis)
that can be used for an ONETOOL run are listed for different data types.

================== ============================================  ======================
Data Type          File Type                                     Extension
================== ============================================  ======================
sample             :ref:`PLINK FAM file <ifam>`                  .fam
phenotype          :ref:`PLINK Phenotype file <ipheno>`          .pheno
variant            :ref:`Variant Call Format (VCF) file <ivcf>`  .vcf
variant            :ref:`PLINK BED/BIM file <ibed_bim>`          .bed/.bim
variant (dosage)   :ref:`IMPUTE2 file <iimpute2>`                .impute2/.impute2_info
sample + variant   :ref:`PLINK PED/MAP file <iped_map>`          .ped/.map
================== ============================================  ======================


Additinal files for a specific analysis:

  - :ref:`S.A.G.E. Parameter file (.par) <ipar>`
  - :ref:`S.A.G.E. Trait genotype probability file (.typ) <ityp>`
  - :ref:`MERLIN MAP file (.map) <imap>`
  - :ref:`Gene SET file (.set) <iset>`
  - :ref:`Script file (.script) <iscript>`


Two main input file sets are 'VCF set' and 'PLINK set'.

The **VCF set** consist of a PLINK format family file (.fam) and a Variant Call
Format file(.vcf).

The **PLINK set** consists of three files (i.e., .fam, .bed, and .bim) that 
are used to run PLINK.

The additional phenotypes and covariates are supported through an optional input
file (.pheno) for both sets of input files.

.. _iscript:

===========
SCRIPT file
===========

ONETOOL also support two different ways to run the program, through a command line and a script file (.script).

A script file includes the input file name(s) and all command-line options 
selected for a ONETOOL run.

.. code-block:: text

    $ onetool --script test.txt

