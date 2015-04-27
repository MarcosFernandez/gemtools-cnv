.. GEMTOOLS CNV PIPELINE documentation master file, created by
   sphinx-quickstart on Fri Apr 10 11:15:38 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################################
GEMTOOLS CNV PIPELINE  Quickstart
#################################

****************************************
1) Download and install the gemtools-cnv
****************************************

If that goes well, you will have a "gemtools" command line tool available.
Please check the `GEMTools homepage <http://gemtools.github.io/>`_ for download
and installation instructions.

******************
2) Mask the genome
******************

In order to run cnv-pipeline the assembly reference must be prepared to
perform in a right way all the cnv pipeline steps. 

The fasta reference must be masked on those regions known to be simple repeats 
(From `Tandem Repeat Finder`_), 
repeats detected by  (From `Repeat Masker`_) and gaps. 

.. _Tandem Repeat Finder: http://en.wikipedia.org/wiki/Tandem_repeat

.. _Repeat Masker: http://en.wikipedia.org/wiki/Interspersed_repeat


Additionaly, the reference genome is fragmented in kmer windows and mapped against itself to
find overrepresented kmer windows that are going to be also masked. These masked regions are not interesting for the CNV study, its masking avoids bias
caused by small repeats sequences. 

To compute read depth, the reference used is also masked 36 bps around any gap or repeat. The rationale behind this
is to avoid deflation of read coverage around repeats or gaps where less reads would have been able to map. Once the pipeline is finished, the kmer masked
reference must be indexed with gem, and the padded reference must be processed by mrcanavar preparation step to perform the cnv calling.:: 

    gemtools cnv-assembly-preparation -f my_fasta_file.fa -t 8 --bed-regions reference_repeats.bed reference_gaps.bed --chr-len chroms.size -kmer-mappings-threshold 20

.. note::
    | All bed files (repeats and gaps) must be without header.
    | The fasta file reference must be raw, without any type of hard masking.
    
.. caution::
    *chroms.size* Files must have two fields separated by tabulator. First one chromosome name and second one chromsome length.

2.1) Outputs
============

    *gemtools cnv-assembly-preparation* creates a set of folders with statistical results, plots, html documentation and output files.
 
    +------------------------+------------------------------------------------------+
    | Folder                 | Contents                                             |    
    |                        |                                                      |
    +========================+======================================================+
    | mask-fasta             | Fasta masked files to be used in the mapping step.   |
    +------------------------+------------------------------------------------------+
    | padded-fasta           | Reference to be used by cnv calling step. This fasta |
    |                        | must be prepared according to 2.2 step.              |
    +------------------------+------------------------------------------------------+
    | html-report            | HTML document with Pipeline stats and plots.         |
    +------------------------+------------------------------------------------------+

2.2) Create configuration reference file
========================================

    Run this command only once for a new reference genome assembly (or versions of it with different repeat masking or different window sizes). 

    In this step, reference genome file (FASTA format) is processed calculating the coordinates and GC content of three window classes (LW, SW, and CW) 
    and saves this information in a "configuration" file. ::

        gemtools cnv-prep -f padded-fasta-reference.fa -g reference_gaps.bed

    **Output**: A file will be generated with *.conf* extension.

    Windows created:

    +------------------------+------------------------------------------------------+
    | Type of Window         | Definition                                           |    
    |                        |                                                      |
    +========================+======================================================+
    | LW. Long Windows       | 5.000 bp of non-masked characters. 1.000 bp slide of |
    |                        | any character.                                       |
    +------------------------+------------------------------------------------------+
    | SW. Short Windows      | 1.000 bp of non-masked characters. 1.000 bp slide of |
    |                        | any character.                                       |
    +------------------------+------------------------------------------------------+
    | CW. Copy Windows       | 1.000 bp of non-masked characters.                   |
    |                        |                                                      |
    +------------------------+------------------------------------------------------+

****************** 
3) Run the pieline
******************

There are two types of pipelines to get copy number calls. 

    The **fast one** (Fast pipeline) performs a mapping using `gem mapper`_ and finally a copy number calling using `mrCaNaVaR`_.
    
    :doc:`fast_pipeline`

    The **long one** needs much more computational time because it looks for PCR Duplicates. Firstly, it runs `bwa mem`_ with default parameters to apply                
    `MarkDuplicates`_ in order to detect PCR artifacts. The BAM alignment file is tranformed to a FASTQ file to run the `gem mapper`_ and call copy number regions
    using `mrCaNaVaR`_.
    
    :doc:`long_pipeline`


.. _gem mapper: http://algorithms.cnag.cat/wiki/The_GEM_library
.. _mrCaNaVaR: http://mrcanavar.sourceforge.net/
.. _bwa mem: http://bio-bwa.sourceforge.net/bwa.shtml
.. _MarkDuplicates: http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

*******************
4) Pipeline outputs
*******************

The differents steps performed by the pipeline outputs a collection of results sorted in a set of folders. A description of the contents per each folder is done in the next table:

+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+ 
| Folder                 | Contents                                             | Fast Pipeline | Long Pipeline | Temporal              |
|                        |                                                      |               |               | Removed once finished | 
+========================+======================================================+===============+===============+=======================+ 
| base-mapping           | bwa alignment file (bam format)                      | No            | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| base-mapping-report    | HTML web document with mapping statistics            | No            | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| rm-dups                | Remove Duplicates (bam format) file                  | No            | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| rm-dups-report         | HTML web document Duplicates statistics              | No            | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| bam-fastq              | Fastq files after BAM conversion                     | No            | Yes           | Yes                   |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| fragment-reads         | Fastq files divided in chunks                        | No            | Yes           | Yes                   |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| split-reads            | Raw reads divided in chunks                          | Yes           | No            | Yes                   |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| chop-reads             | 36 bp Chunks of Fastq reads                          | Yes           | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| cnv-map                | Mappings in GEM format (MAP)                         | Yes           | Yes           | Yes                   |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| map-sam                | Mappings in SAM format                               | Yes           | Yes           | No                    |  
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| map-stats              | JSON mapping stats files                             | Yes           | Yes           | No                    | 
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| mrcanavar              | Copy number outputs                                  | Yes           | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| cn-distribution        | Copy number distribution plots                       | Yes           | Yes           | No                    | 
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+
| html-doc               | Pipeline documentation in HTML web format            | Yes           | Yes           | No                    |
+------------------------+------------------------------------------------------+---------------+---------------+-----------------------+

.. hint::
    If you are just interested on copy number results and HTML documentations then you could remove mapping and chopping files in order to save
    disk space.

4.1) Copy Number Results
========================

In folder ``mrcanavar`` are located the copy number calls. The files generated are:

**CNV CALLS**

+------------------------------------------+------------------------------------------------------+
| File                                     | Contents                                             |
+==========================================+======================================================+
| mysample_mrcanavar.calls.copynumber.bed  | Main pipeline file result. Copy Number for non       |
|                                          | overlapping windows (Copy windows)                   |
+------------------------------------------+------------------------------------------------------+
| mysample_mrcanavar.calls.cw_norm.bed     | Read Depth normalized by GC Content and control      |
|                                          | windows notification for copy windows                |
+------------------------------------------+------------------------------------------------------+
| mysample_mrcanavar.calls.sw_norm.bed     | Read Depth normalized by GC Content and control      |
|                                          | windows notification for short windows               |
+------------------------------------------+------------------------------------------------------+
| mysample_mrcanavar.calls.lw_norm.bed     | Read Depth normalized by GC Content and control      |
|                                          | windows notification for long windows                |
+------------------------------------------+------------------------------------------------------+
| mysample_mrcanavar.calls.log             | Brief description about read depth average and       |
|                                          | standard deviation of copy, short and long windows   | 
+------------------------------------------+------------------------------------------------------+
| mysample_mrcanavar.calls                 | Pipeline general information                         |
+------------------------------------------+------------------------------------------------------+


**READ DEPTH**
 
+------------------------------------------+-------------------------------------------------------+
| File                                     | Contents                                              |
+==========================================+=======================================================+
| mysample_mrcanavar.depth.cw_norm.bed     | Read Depth Normalized by GC Content for copy windows  |
+------------------------------------------+-------------------------------------------------------+
| mysample_mrcanavar.depth.sw_norm.bed     | Read Depth Normalized by GC Content for short windows |
+------------------------------------------+-------------------------------------------------------+
| mysample_mrcanavar.depth.lw_norm.bed     | Read Depth Normalized by GC Content for long windows  |
+------------------------------------------+-------------------------------------------------------+
| mysample_mrcanavar.depth.cw.txt          | Raw Read Depth for copy windows                       |
+------------------------------------------+-------------------------------------------------------+
| mysample_mrcanavar.depth.sw.txt          | Raw Read Depth for short windows                      |
+------------------------------------------+-------------------------------------------------------+
| mysample_mrcanavar.depth.lw.txt          | Raw Read Depth for long windows                       |
+------------------------------------------+-------------------------------------------------------+     
| mysample_mrcanavar.depth                 | Raw Read Depth binary file                            |  
+------------------------------------------+-------------------------------------------------------+  
                 
4.2) Pipeline HTML report
=========================

Located in ``html-doc`` are found two document files, **html** web document and **json** text document. These files contains statistical values and plots of each pipeline step. Plot images are ``png`` files located in ``cn-distribution``. **HTML** web document have references to images located in ``cn-distribution`` folder.

**Long Pipeline outputs:**

    .. glossary::

        base-mapping-report
            Mapping statistics after **bwa mem** mapping.

        rm-dups-report
            Remove duplicates statistics after performing **MarkDuplicates** from `MarkDuplicates`_

******************
5) Important Notes 
******************

5.1) SAMTOOLS
=============

The pipeline expects to find samtools installed on the system. Try to get the latest samtools from their github repository 
(`https://github.com/samtools/samtools`_ – clone or download and call make to build it). 
       
The latest version is multi-threaded (i.e. samtools view --help will show a -@ paramter). 

Also, see if you have pigz installed in the system you try to run gemtools on. pigz is a parallel compressor and the pipeline makes use 
of it if it is available. It will speed up compression steps a lot!

5.2) BWA
========

The pipeline also expects to find bwa mapper installed on the system. Try to get the latest samtools from their github repository 
(`https://github.com/lh3/bwa`_ – clone or download and call make to build it). 


5.2) PICARD TOOLS
=================

The pipeline uses picard tools when performing the **long version**. Download the package from (`https://github.com/broadinstitute/picard`_) or 
follow `http://broadinstitute.github.io/picard/`_ instructions.


5.3) R Package
==============

R package must be installed in your system. Get last version from `http://www.r-project.org/`_ .


5.4) BEDTOOLS
=============

The pipeline uses BEDTOOLS package to perform some of the steps. Get last version from `https://code.google.com/p/bedtools/`_ and follow install instructions.




.. _https://github.com/samtools/samtools: https://github.com/samtools/samtools

.. _https://github.com/lh3/bwa: https://github.com/lh3/bwa

.. _https://github.com/broadinstitute/picard: https://github.com/broadinstitute/picard

.. _http://broadinstitute.github.io/picard/: http://broadinstitute.github.io/picard/

.. _http://www.r-project.org/: http://www.r-project.org/

.. _https://code.google.com/p/bedtools/: https://code.google.com/p/bedtools/



    





