.. JIP GEM CNV PIPELINE documentation master file, created by
   sphinx-quickstart on Mon Apr 20 10:47:45 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

JIP GEM CNV PIPELINE  Quickstart
================================

JIP GEM CNV PIPELINE is a set of jip scripts that uses the **JIP pipeline system** to perform all pipeline steps, controlling each step dependency and output. The pipeline expects to find **JIP pipeline system** and **gemtools-cnv** installed on your system.

1) Download and install JIP PIPELINE system
-------------------------------------------

The JIP pipeline system is a python library and a set of command line utilities that allows you to create batch-process based computational pipeline that can be submitted and managed on a computer cluster or on your local machine.

Install the jip pipeline system following the instructions found in: `http://pyjip.readthedocs.org/en/latest/`_

.. _http://pyjip.readthedocs.org/en/latest/: http://pyjip.readthedocs.org/en/latest/


2) Download and install gemtools-cnv
------------------------------------

Download and install **gemtools-cnv** in your system. If that goes well, you will have a "gemtools-cnv" command line tool available.


3) Masking the reference
------------------------

Reference masking in gaps and repeats locations. It also creates an *index* for `bwa`_ mapper, *index* for the masked reference for `gem mapper`_ and a *configuration* file from the padded masked reference to perform copy number calls using `mrCaNaVaR`_. ::

    cnv_assembly_preparation.jip -f my_reference.fa -t 8 -r repeat_masker.bed tandem_repeats.bed -g gaps.bed -a chromInfo.txt -o out_dir

.. note::
   
    **Chromosome Length**

        *chromInfo.txt* must be a file with two fields per row separated by *tabulators*. 

        **First** field must be chromosome or conting **name** and **second** field its **length** in base pairs.
     


.. _gem mapper: http://algorithms.cnag.cat/wiki/The_GEM_library
.. _mrCaNaVaR: http://mrcanavar.sourceforge.net/
.. _bwa: http://bio-bwa.sourceforge.net/bwa.shtml 


4) Run configuration file
-------------------------

Creates a configuration json file for the cnv pipeline. This **JSON** file is used to manage the pipeline. **By default**, is performed the **long pipeline** (removing PCR duplicates) and assumes that input data is **paired end**. ::

    createConfigurationFile.py -bwa-reference my_reference.fa -gem-index my_reference_kmer_mask_mask-fasta.gem -reference-conf my_reference_pad36_padded-fasta.conf  --json-file my_pipeline_cfg.json -- submit --dry --show


These are the main parameters:

    .. glossary::

        --no-duplicates
            If specified, do not run remove duplicates steps. (**Fast Pipeline**)

        -se
            If specified, Input FASTQ data is treated as single end.

        --chunks
            Number of chunks to fragment FASTQ input data to later perform a mapping job per chunk. The more chunks, the more parallelization. Default 50.




4.1 Mandatory arguments
_______________________

    There are four arguments which are mandatory when creating the pipeline configuration.

    Mandatory arguments:

    .. glossary::

        -bwa-reference
            Path to the fasta index were are located the bwa index files.

        -gem-index
            Path to the gem index file of the fasta reference.

        -reference-conf 
            Path to the mrCaNaVar configuration of the reference file

        --json-file 
            Pipeline configuration file in *json* format.


    .. note::

        **createConfigurationFile.py**

        Check script parameters help to adjust the number of threads, job times, and other specific parameters per each step.




5) Run the pipeline
-------------------

Creation and submission of a set of jip jobs to perform the pipeline running according to a configuration file.  ::

    cnv_pipeline.jip  -c my_pipeline_cfg.json -i my_fastq_dir/*_1.fastq.gz -o ./MY_SAMPLE_DIR/ -n SAMPLE_NAME  -d SAMPLE_DESCRIPTION --gz -- submit -P development --dry --show


.. seealso::

    To know more about the CNV PIPELINE results check **GEMTOOLS CNV PIPELINE** documentation.




