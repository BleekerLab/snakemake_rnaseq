Welcome to Snakemake RNA-seq pipeline documentation
====================================================

Contents
========

This is the main documentation of the Snakemake RNA-seq pipeline. 
This pipeline takes single or paired-end mostly Illumina mRNA-seq reads and generates:    
1. A quality report (multiQC format),
2. A mapping summary, 
3. The gene *raw* read counts,
4. The gene *scaled* read counts

The gene raw counts can be directly used for differential gene expression with `DESeq2 <https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`. 

A dedicated tutorial, `available here <https://scienceparkstudygroup.github.io/rna-seq-lesson/06-differential-analysis/index.html>` can be followed to perform this analysis. 

.. toctree::
   :maxdepth: 2

   installation
   crunchomics
   license



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
