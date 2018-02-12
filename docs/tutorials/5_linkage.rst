.. _sec-linkage:

Linkage analysis in ONETOOL is composed of two parts, two-point model-based linkage 
and two-point model-free linkage accounting for linkage disequilibrium (LD).   

.. _tlodlink:

===========
Model-based
===========

LODLINK module performs model-based lod score calculations for two-point linkage between a main
trait (that follows Mendelian transmission and has either two or three types) and each of the 
variants in the input genotype file.  The output from SEGREG model (.typ file) is used as
input to specify the main trait model. 

LODLINK in S.A.G.E. is used which uses the genotype/phase elimination algorithms proposed 
by Lange and Boehnke (1983) and Lange and Goradia (1987), together with other enhancements, 
to perform relatively fast exact linkage calculations.

.. code-block:: text

   $ onetool --fam test_miss00.fam --vcf test_miss0.vcf --lodlink --typ test_segreg.typ

It generates 2 output files:

  - Summary output file (.lodlink.sum) - contains lod scores and results of linkage and linkage homogeneity tests. Results in this file are based on calculations done on the pedigree data file as a whole.
  - Detailed output file (.lodlink.det) - contains results on per family.


.. _tmerlin:

===========
Model-free
===========

ONETOOL uses Merlin for multipoint model-free linkage analysis accounting for linkage 
disequilibrium (LD), providing the Kong and Cox (1997) NPL statistics.

.. code-block:: text

   $ onetool --fam test_miss00.fam --vcf test_miss0.vcf --merlin --map test_miss0.map

It generates a output file with the extention .nonparametric.tbl.

.. note:: The same restriction applies to the size of family as in MERLIN. 

