3.2) Long CNV Pipeline
______________________

Long Copy Number pipeline version. 

Firstly, looks for PCR Duplicates using `bwa mem`_ for mapping reads against a *raw fasta reference* and `MarkDuplicates`_ to perform PCR artifacts removal.

Secondly, transforms the **BAM** alignment file, free of duplicates, to FASTQ and performs a mapping with `gem`_ to make a call of copy number.

.. _bwa mem: http://bio-bwa.sourceforge.net/bwa.shtml

.. _MarkDuplicates: http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates 

.. _gem: http://algorithms.cnag.cat/wiki/The_GEM_library

.. note::
    The process of looking for PCR Duplicates adds a significant cost on time and computing resources. You should consider this, if you
    already know that the percentage of duplicates is very low or 0 in your libraries then it is better to run the **fast cnv pipeline**.

.. warning::
    `PCR Duplicates`_ are artifacts generated at *library level*. Then, you must use as input data all FASTQ files that comes from the same 
    DNA library sample preparation. Duplicates reads can be found in different lanes and files. 


.. _PCR Duplicates: http://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/


3.2.0) Index reference
______________________

As `bwa mem`_ is going to be used as aligner we must index our raw fasta reference. Run *index* command only once for a new reference genome assembly ::

   bwa index my_raw_reference.fa

.. note:: 
    For removing duplicates it is not necessary to mask the reference as we must do for the copy number mapping. Indeed, we must map against the raw reference 
    to get a good estimation of the real duplicate percentage and to not ignore those duplicate reads involved in repeat regions. 

3.2.1) Basic Mapping
____________________

Basic mapping with standard `bwa mem`_ parameters. 

.. glossary::

    Input file detection
        If you do not specify ``--single`` to disable read pairing, it will look automatically for the second pair file in case you only specify one file. 
        For that to work, the second file has to end with .2 or _2, with the file extension .fastq or .txt (+ .gz for compressed files). 
    
    Picard Tools
        Basic Mapping is going to use Picard Tools for some of the steps to perform. You must provide the path were Picard Tools is located.    

::

    gemtools-cnv basic-mapping -f my_fastq_1.fastq.gz -bwa-ref my_raw_reference.fa -T 8 -sample-description description -picard-path /path/picard/tools/ -tmp-folder $tmp 

3.2.2) Remove Duplicates
________________________

Remove Duplicates looks for PCR artifacts in all input BAMS and outputs a new *BAM* alignment file free of duplicates. It also generates a HTML report with Remove Duplicates statistics.

.. warning::
    In order to detect all possible PCR duplicates you should use as input all reads for given library. You can specify a list of all *BAM* alignment files.

:: 

    gemtools-cnv remove-duplicates -f file_1.bam file_2.bam -picard-path /path/picard/tools/ -java-heap 25g -tmp-folder $tmp -sample-description description

3.2.3) BAM alignment file to chopped FASTQ
__________________________________________

Transforms a *BAM* alignment file to a set of *FASTQ* sequencing files. The idea is to *speed up* the mapping process creating a mapping job per each *FASTQ* generated. You should decide the number of chunks according to your computing resources.

.. glossary::

    BAM treatment
        The bam file will be treated as *paired end* unless you specify ``-single-end`` argument.
    
    Chopping
        Reads are chopped to 36 base pairs to not loose mapping information around masked regions. 
    
    Filtering
        By default, the first ten bases are removed because unexpected tendencies of nucletide distribution are usually found. 
        Nonetheless, this pattern should be confirmed through `fastqc`_ analysis. 

::

    gemtools-cnv bam-2-fastq -f my_file_rm-dups.rmdup.bam -split-times 20 -n sample_name --threads 8


.. _fastqc: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/


3.2.4) Alignment to the masked reference
----------------------------------------

Aligns reads against a *masked reference genome* using `gem mapper`_. Maps in *single en mode*, with the goal of generating a set of mappings to later perform a copy number calling based on Read Depth. ::

    gemtools-cnv cnv-mapping -f my_fastq_chop-reads.part-1.1.fq -i my_reference_kmer_mask_mask-fasta.gem -n my_mapping_name  --mappin-stats-json -T 8

.. _gem mapper: http://algorithms.cnag.cat/wiki/The_GEM_library

3.2.5) Copy Number Calling
--------------------------

Runs `mrCaNaVaR`_  to calculate read depth for a set of different kind of windows spanning all the genome. (Long Windows, Short Windows, Copy Windows). After a GC content correction a copy number call is performed.

It also generates a HTML Pipeline report with different statistics and plots. :: 

    gemtools-cnv cnv-call -d /mappings/map-sam/ --conf_file my_reference_pad36_padded-fasta.conf --gz --no-duplications -sample-description description -n my_name

.. _mrCaNaVaR: http://mrcanavar.sourceforge.net/













 



