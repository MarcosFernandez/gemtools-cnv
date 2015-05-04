#!/usr/bin/env python
"""Python wrapper around the GEM2 mapper that provides
ability to feed data into GEM and retrieve the mappings"""
import os
import sys
import logging
import subprocess

import tempfile
import files
from . import utils

import pkg_resources
import splits
import gem.filter as gemfilter
import gem.gemtools as gt
import gem
import gem.junctions
import gem.duplications
import __builtin__


LOG_NOTHING = 1
LOG_STDERR = 2
LOG_FORMAT = '%(asctime)-15s %(levelname)s: %(message)s'
# add custom log level
LOG_GEMTOOLS = logging.WARNING
logging.addLevelName(LOG_GEMTOOLS, "")
logging.basicConfig(format=LOG_FORMAT, level=logging.WARNING)
gemtools_logger = logging.getLogger("gemtools")
gemtools_logger.propagate = 0
gemtools_logger.setLevel(LOG_GEMTOOLS)

def log_gemtools(message, *args, **kws):
     gemtools_logger.log(LOG_GEMTOOLS, message, *args, **kws)

gemtools_logger.gt = log_gemtools

logging.gemtools = gemtools_logger

__parallel_samtools = None

class GemtoolsFormatter(logging.Formatter):
    info_fmt = "%(message)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        format_orig = self._fmt
        if record.levelno == LOG_GEMTOOLS:
            self._fmt = GemtoolsFormatter.info_fmt
        result = logging.Formatter.format(self, record)
        self._fmt = format_orig
        return result

gemtools_formatter = GemtoolsFormatter('%(levelname)s: %(message)s')
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
console.setFormatter(gemtools_formatter)
gemtools_logger.addHandler(console)
logging.gemtools.level = logging.WARNING

# default logger configuration
log_output = LOG_NOTHING


default_splice_consensus = [("GT", "AG")]
extended_splice_consensus = [("GT", "AG"),
    ("GC", "AG"),
    ("ATATC", "A."),
    ("GTATC", "AT")]


#default filter
default_filter = "same-chromosome,same-strand"

## use the bundled executables
use_bundled_executables = True

## max mappings to replace mapping counts for + and ! summaries
_max_mappings = 999999999

## filter to work around GT-32 and #006 in gem-2-gem
__awk_filter = ["awk", "-F", "\t", '{if($4 == "*" || $4 == "-"){print $1"\t"$2"\t"$3"\t0\t"$5}else{if($4 == "!" || $4 == "+"){print $1"\t"$2"\t"$3"\t' + str(_max_mappings) + '\t"$5}else{print}}}']
## filter to put before teh pairaligner if qualities are ignored
## current version has a bug where it does not expect to see the quality column
__awk_pair_quality_fix = ["awk", "-F", "\t", '{print $1"\t"$2"\t"$4"\t"$5}']
## The current version of gem does not work without qualities properly.
## We workaround by inserting dummy qualities I
__awk_gem_2_sam_quality_fix = [
    'awk', '-F', '\t',
    '{x=gensub(/\w/, "I", "g", $2); print $1"\t"$2"\t"x"\t"$4"\t"$5}'
]


class execs_dict(dict):
    """Helper dictionary that resolves bundled binaries
    based on the configuration. We check first for GEM_PATH
    environment variable. If its set, and points to a directory with
    the executable, the path to that executable is returned.
    Next, we check the use_bundled_executable flag. If that is true(default)
    the path to the bundled executable is returned.
    If nothing is found, the plain executable name is returned and we
    assume it can be found in PATH

    """
    def __getitem__(self, item):
        # check if there is an environment variable set
        # to specify the path to the GEM executables
           
        base_dir = os.getenv("GEM_PATH", None)
        if base_dir is not None:
            file = "%s/%s" % (base_dir, item)
            if os.path.exists(file):
                logging.debug("Using binary from GEM_PATH : %s" % file)
                return file

        if use_bundled_executables and pkg_resources.resource_exists("gem", "gembinaries/%s" % item):
            f = pkg_resources.resource_filename("gem", "gembinaries/%s" % item)
            logging.debug("Using bundled binary : %s" % f)
            return f
        # try to find from static distribution
        if use_bundled_executables and len(sys.argv) > 0:
            try:
                base = os.path.split(os.path.abspath(sys.argv[0]))[0]
                binary = base + "/" + item
                if os.path.exists(binary):
                    logging.debug("Using bundled binary : %s" % binary)
                    return binary
            except Exception:
                pass

        logging.debug("Using binary from PATH: %s" % item)
        return dict.__getitem__(self, item)

## paths to the executables
executables = execs_dict({
    "gem-indexer": "gem-indexer",
    "gem-mapper": "gem-mapper",
    "gem-2-gem": "gem-2-gem",
    "gem-2-sam": "gem-2-sam",
    "samtools": "samtools",
    "gem-info": "gem-info",
    "gem-retriever": "gem-retriever",
    "gem-rna-tools": "gem-rna-tools",
    "gt.filter": "gt.filter",
    "gt.map.2.sam": "gt.map.2.sam",
    "gt.mapset": "gt.mapset",
    "gt.gtfcount": "gt.gtfcount",
    "gt.stats": "gt.stats",
    "bam2fastq": "bam2fastq",
    "FastqSplitter": "FastqSplitter",
    "kmermaker": "kmermaker",
    "gt-scorereads": "gt.scorereads",
    "mrcanavar": "mrcanavar",
    "copyNumberDistribution.R": "copyNumberDistribution.R",
    "makeHistogramCounts.R": "makeHistogramCounts.R",
    "makeBedIntervalsAssembly": "makeBedIntervalsAssembly",
    "kmerCount": "kmerCount",
    "galculator": "galculator"
    })


def loglevel(level):
    """Simple way to set the current log level globally for the root logger.
    Accepts either 'debug','info','warning', 'error'

    Log levels debug also ensures executable output is written to stderr

    level -- one of debug, info, warn, error
    """
    global log_output
    numeric_level = level
    if isinstance(level, basestring):
        numeric_level = getattr(logging, level.upper(), None)

    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logging.basicConfig(level=numeric_level)
    logging.getLogger().setLevel(numeric_level)
    logging.gemtools.level = numeric_level


# cleanup functions
def _cleanup_on_shutdown():
    gem.utils.terminate_processes()
    gem.files._cleanup()

import atexit
atexit.register(_cleanup_on_shutdown)

def _prepare_index_parameter(index, gem_suffix=True):
    """Prepares the index file and checks that the index
    exists. The function throws a IOError if the index file
    can not be found.

    index      -- the path to the index file
    gem_suffix -- if true, the function ensures that the index ends in .gem,
                  otherwise, it ensures that the .gem suffix is removed.

    """
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    file_name = index

    if not file_name.endswith(".gem"):
        file_name = file_name + ".gem"

    if not os.path.exists(file_name):
        raise IOError("Index file not found : %s" % file_name)

    if gem_suffix:
        if not index.endswith(".gem"):
            index = index + ".gem"
    else:
        if index.endswith(".gem"):
            index = index[:-4]
    return index

def _prepare_splice_consensus_parameter(splice_consensus):
    """Convert the splice consensus tuple to
    valid gem parameter input.

    If the given splice_consensus is None, the
    default splice consensus is used

    splice_consensus -- if string, it is just passed on, if its a list of tuples like ("A", "G") it
                        is translated into gem representation
    """
    if splice_consensus is None:
        splice_consensus = default_splice_consensus
    if isinstance(splice_consensus, basestring):
        splice_cons = splice_consensus
    else:
        ## translate the splice consensus tupel structure
        splice_cons = ",".join(['%s+%s' % (x[0], x[1]) for x in splice_consensus])
    return splice_cons


def _prepare_quality_parameter(quality, input=None):
    """Prepare and returnn the quality parameter for gem runs. If the
    input is None, qualities are disabled, if it is a string and
    valid parameter it is returned as is. Otherwise it as to be 33
    or 64 and will be translated to a valid gem parameter

    quality -- the quality offset, 33|64, None to ignore or a string that
                is either offset-33|offset-64|ignore
    input   -- optional input that can be checkd for a quality parameter. This
               currently works only for gemtools.TemplateIterator or gemtools.InputFile
               instances
    """
    if quality is None and isinstance(input, (gt.InputFile)):
        quality = input.quality

    ## check quality
    if quality is not None and quality in ["offset-33", "offset-64", "ignore"]:
        return quality
    if quality is not None and quality not in ["none", "ignore"]:
        i = int(quality)
        if i not in [33, 64]:
            raise ValueError("%s is not a valid quality value, try None, 33 or 64" % (str(quality)))
        quality = "offset-%d" % int(quality)
    else:
        quality = 'ignore'

    return quality

def _prepare_output(process, output=None, quality=None, bam=False):
    """Creates a new gem.gemtools.Inputfile from the given process.
    If output is specified, the function blocks and waits for the process to finish
    successfully before the InputFile is created on the specified output.

    Otherwise, a stream based InputFile is created using the process stdout
    stream.

    If quality is specified it is passed on to the input file.

    process -- the process
    output  -- output file name or None when process stdout should be used
    quality -- optional quality that is passed to the InputFile
    """
    # wrap process list
    if not isinstance(process, (list, tuple)):
        process = [process]

    if output is not None:
        # we are writing to a file
        # wait for the process to finish
        if process is not None:
            for k, p in enumerate(process):
                if p is not None and p.wait() != 0:
                    raise gem.utils.ProcessError("Execution failed!")
                # close streams next process in
                if (k + 1) < len(process):
                    next_process = process[k + 1]
                    if next_process is not None and \
                       next_process.stdin is not None:
                        next_process.stdin.close()
        logging.debug("Opening output file %s" % (output))
        if output.endswith(".bam") or bam:
            return gt.InputFile(gem.files.open_bam(output), quality=quality, process=process[0])
        return gt.InputFile(output, quality=quality, process=process[0])
    else:
        logging.debug("Opening output stream")
        if bam:
            return gt.InputFile(gem.files.open_bam(process[0].stdout), quality=quality, process=process[0])
        ## running in async mode, return iterator on
        ## the output stream
        return gt.InputFile(process[0].stdout, quality=quality, process=process[0])

def _prepare_bwa_reference_parameter(reference):
    """Prepares the bwa reference file and checks that the index
    exists. The function throws a IOError if the reference file
    can not be found.

    reference  -- the path to the index file
    
    """
    if reference is None:
        raise ValueError("No valid BWA reference specified!")
    if not isinstance(reference, basestring):
        raise ValueError("BWA reference must be a string")
    
    bwa_index = reference + ".sa"

    if not os.path.exists(reference):
        raise IOError("Reference fasta file not found : %s" % reference)

    if not os.path.exists(bwa_index):
        raise IOError("Index file not found : %s" % bwa_index)

    return reference
    
def _prepare_sam_directory_parameter(samfile):
   """ Return directory where is located mapping file.
   Checks if the sam files files exists. The function throws a IOError if the
   samfile can not be found
   
   samfile - sam mapping file from which calculate read depth and copy number
   
   """
   if not os.path.exists(samfile):
        raise IOError("Sam file not found : %s" % samfile)

   return os.path.dirname(samfile)
   
def _prepare_copywindows_bed(mr_canavar_dir):
    """ Returns copy windows bed file 
    
    mr_canavar_dir -- mrCanavar Director results
    """
    callsDir = os.path.dirname(mr_canavar_dir)
    
    copyWindowsBed = gem.utils.runFindCommand(callsDir,"*.calls.cw_norm.bed").rstrip('\n')
    
    if not os.path.exists(copyWindowsBed):
        raise IOError("copy windows file not found in: %s" % mr_canavar_dir)
   
    return copyWindowsBed
    
def _prepare_copyNumber_bed(mr_canavar_dir):
    """ Returns copy number bed file 
    
    mr_canavar_dir -- mrCanavar Director results
    """
    callsDir = os.path.dirname(mr_canavar_dir)
    
    copyNumberBed = gem.utils.runFindCommand(callsDir,"*.calls.copynumber.bed").rstrip('\n')
    
    if not os.path.exists(copyNumberBed):
        raise IOError("copy number bed file not found in: %s" % mr_canavar_dir)
   
    return copyNumberBed
    
def _prepare_calls_log(mr_canavar_dir):
    """ Returns calls log file 
    
    mr_canavar_dir -- mrCanavar Director results
    """
    callsDir = os.path.dirname(mr_canavar_dir)
    
    copyNumberBed = gem.utils.runFindCommand(callsDir,"*.calls.log").rstrip('\n')
    
    if not os.path.exists(copyNumberBed):
        raise IOError("mrcanvar logs file not found in: %s" % mr_canavar_dir)
   
    return copyNumberBed


def validate_executables():
    """Validate the gem executables and
    print the paths to the executables in use
    """
    for exe, path in executables.items():
        path = executables[exe]
        exe_path = utils.which(executables[exe])
        found = exe_path is not None
        if found:
            print >> sys.stderr, "Executable '%s' (%s) : %s" % (exe, path, exe_path)
        else:
            print >> sys.stderr, "Executable '%s' (%s) : Unknown" % (exe, path)
            
def is_gem_index(fasta_file):
    """Checks for the existance of the file reference, return the absolute path to the index reference """
    gem_index = ""
    
    if fasta_file.endswith(".fa"):
       gem_index += fasta_file.replace(".fa",".gem")
    else:
       gem_index += fasta_file.replace(".fasta",".gem") 
        
    if os.path.exists(gem_index):
        return os.path.abspath(gem_index)
    return ""

def mapper(input, index, output=None,
           mismatches=0.04,
           delta=0,
           quality=33,
           quality_threshold=26,
           max_decoded_matches=20,
           min_decoded_strata=1,
           min_matched_bases=0.80,
           max_big_indel_length=15,
           max_edit_distance=0.20,
           mismatch_alphabet="ACGT",
           trim=None,
           unique_mapping=False,
           threads=1,
           extra=None,
           key_file=None,
           force_min_decoded_strata=False,
           compress=False
           ):
    """Start the GEM mapper on the given input.
    If input is a file handle, it is assumed to
    provide fastq entries. If input is a string,
    it is checked for its extension. In case of a
    .map file, the input is converted from gem format
    to fastq and passed to the mapper.

    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- A ReadIterator with the input
    output -- output file name
    index -- valid GEM2 index
    mismatches - number or % mismatches, default=0.04
    delta -- strata after best <number> (default=0)
    quality -- one of 'ignore'|'offset-33'|'offset-64' defaults to offset-33
    quality_threshold <number> -- (default=26, that is e<=2e-3)
    max_edit_distance -- max edit distance, 0.20 per default
    max_decoded_matches -- maximum decoded matches, defaults to 20
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 1
    min_matched_bases -- minimum number (or %) of matched bases, defaults to 0.80
    trim -- tuple or list that specifies left and right trimmings
    extra -- list of additional parameters added to gem mapper call
    """

    ## prepare inputs
    index = _prepare_index_parameter(index)
    quality = _prepare_quality_parameter(quality, input)

    # if delta >= min_decoded_strata and not force_min_decoded_strata:
    #     logging.warning("Changing min-decoded-strata from %s to %s to cope with delta of %s" % (
    #         str(min_decoded_strata), str(delta + 1), str(delta)))
    #     min_decoded_strata = delta + 1
    if compress and output is None:
        logging.warning("Disabeling stream compression")
        compress = False

    if compress and not output.endswith(".gz"):
        output += ".gz"

    ## prepare the input
    pa = [executables['gem-mapper'], '-I', index,
          '-q', quality,
          '-m', str(mismatches),
          '-s', str(delta),
          '--max-decoded-matches', str(max_decoded_matches),
          '--min-decoded-strata', str(min_decoded_strata),
          '--min-matched-bases', str(min_matched_bases),
          '--gem-quality-threshold', str(quality_threshold),
          '--max-big-indel-length', str(max_big_indel_length),
          '--mismatch-alphabet', mismatch_alphabet,
          '-T', str(threads)
    ]

    if unique_mapping:
        pa.append("--unique-mapping")

    if max_edit_distance > 0:
        pa.append("-e")
        pa.append("%s" % str(max_edit_distance))

    ## extend with additional parameters
    _extend_parameters(pa, extra)

    trim_c = [executables['gem-2-gem'], '-c', '-T', str(threads)]
    if trim is not None:
        ## check type
        if not isinstance(trim, (list, tuple)) or len(trim) != 2:
            raise ValueError("Trim parameter has to be a list or a tuple of size 2")
        input = gemfilter.trim(input, trim[0], trim[1], append_label=True)

    # workaround for GT-32 - filter away the !
    # build list of tools
    tools = [pa]
    if unique_mapping:
        tools.append(__awk_filter)

    if trim is not None:
        tools.append(trim_c)

    # convert to genome coordinates if mapping to transcriptome
    if key_file is not None:
        convert_to_genome = [executables['gem-rna-tools'],
                             'transcriptome-2-genome',
                             '-k', key_file,
                             '--threads', str(max(1, threads / 2))
                             ]
        tools.append(convert_to_genome)

    if compress:
        gzip = _compressor(threads=max(1, threads / 2))
        tools.append(gzip)

    raw = False
    if isinstance(input, gt.InputFile) and input.raw_sequence_stream():
        raw = False
        pa.append("-i")
        pa.append(input.filename)
        input = None

    ## run the mapper
    process = utils.run_tools(tools, input=input, output=output, name="GEM-Mapper", raw=raw)
    return _prepare_output(process, output=output, quality=quality)


def bwaMapper(input, reference, output=None,
           threads=1,
           compress=False,
           name="",
           sort_memory="768M"
           ):
    """Start the BWA mapper on the given input.
    If input is a file handle, it is assumed to
    provide fastq entries. If input is a string,
    it is checked for its extension. 

    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- A ReadIterator with the input
    output -- output file name
    reference -- bwa reference
    name -- namfe of the file to be mapped
    """
    ## prepare inputs
    reference = _prepare_bwa_reference_parameter(reference)



    ## prepare the input
    bwaList = ['bwa','mem',
          '-t', str(threads),
          '-M',
          reference,
          input[0]
    ]
 
    
    if len(input) > 1:
        bwaList.append(input[1])      
    
    # build list of tools
    tools = [bwaList]
    
    #sam to bam 
    sam2bam_p = _check_samtools("view", threads=threads, extend=["-S", "-b", "-h"])
   
    sam2bam_p.append('-')
    tools.append(sam2bam_p)
    
    
    #sort bam
    bamToSort = ['samtools', 'sort',
                '-@',str(threads),
                '-T','deleteme.' + name,
                '-O','bam',
    ]

    tools.append(bamToSort) 
    
    process = utils.run_tools(tools, input=None, output=output, name="BWA-MAP")
    if process.wait() != 0:
        raise ValueError("Error while running BWA MEM Mapping")
    
def runBamSummaryMetrics(input=None,output=None,picard_tools_path=None,java_heap="25g",tmp_folder="/tmp",reference=None):
    """ Perform Bam Summary Metric using PicardTools
    input -- BAM file to build stats
    output -- txt metrics file ouput
    picard_tools_path -- Path were is found .jar PicardTools applications
    java_heap -- Java Memory Heap
    tmp_folder -- Temporary folder
    """
    
     ## prepare the input
    pa = ['java','-Xmx' + java_heap,
          '-Djava.io.tmpdir=' + tmp_folder,
          '-jar',picard_tools_path + '/CollectAlignmentSummaryMetrics.jar',
          'REFERENCE_SEQUENCE=' + reference,
          'METRIC_ACCUMULATION_LEVEL=ALL_READS',
          "OUTPUT=" + output,
          "VALIDATION_STRINGENCY=SILENT",
          "INPUT="+input
         ]
       
    # build list of tools
    tools = [pa]  

    process = utils.run_tools(tools, input=None, output=None, name="run-bam-summary-metrics")
    if process.wait() != 0:
        raise ValueError("Error while running Bam Summary Metrics")
		

def mergeBams(input=None,output=None):
    """ Merge Bam files 
    input -- List of bam files
    output -- Merged bam 
    """
    
    #merge command
    merge = ['samtools','merge',output]
    for i in input:
        merge.append(i)
        
    #build list of tools
    tools = [merge]
    
    process = utils.run_tools(tools, input=None, output=None, name="merge-bams")
    if process.wait() != 0:
        raise ValueError("Error while running merge bam")
    
 
def markDuplicates(input,picard_tools_path=None,java_heap="25g",tmp_folder="/tmp",output=None):
    """ Remove PCR Artefacts from a sorted mapping bam file 
    input -- A ReadIterator with the input
    picard_tools_path -- Path were is found .jar PicardTools applications
    java_heap -- Java Memory Heap
    tmp_folder -- Temporary folder 
    output -- output file name    
    """
    
    ## prepare the input
    pa = ['java','-Xmx' + java_heap,
          '-Djava.io.tmpdir=' + tmp_folder,
          '-jar',picard_tools_path + '/MarkDuplicates.jar',
          'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
          'MAX_RECORDS_IN_RAM=1500000',
          'METRICS_FILE=' + output[0],
          'REMOVE_DUPLICATES=true',
          'ASSUME_SORTED=true',
          'VALIDATION_STRINGENCY=SILENT',
          'CREATE_INDEX=true',
          'INPUT=' + input,
          'COMPRESSION_LEVEL=9',
          'OUTPUT=' + output[1]
         ]
       
    # build list of tools
    tools = [pa]   
    
    process = utils.run_tools(tools, input=input, output=output, name="PCR-DUP")
    return _prepare_output(process, output=output[1], quality=33, bam=True)
          
def bamToFastq(input,output=None):
    """ Translates a BAM mapping file to a FASTQ one
    input -- A ReadIterator with input
    output -- output file name
    """
    
     ## prepare the input
    pa = [executables['bam2fastq'],
          '-o',output,
          input]
          
    # build list of tools
    tools = [pa]
    
    process = utils.run_tools(tools, input=None, output=None, name="BAM-FASTQ")
    return _prepare_output(process, output=output)
                 
def fragmentFastq(inputs,output=None,split_times=100,gz=False,prefix="",threads=1):
    """ Chops a FASTQ file to 36 bp
    inputs -- list of fastq files to be processed,  one in case of SE, two for the case of PE
    output -- list of output fastq files 
              n = number of splits
              SE: fq.part-1,...,fq.part-i,...,fq.part-n,    
              PE: 1.fq.part-1,2.fq.part-1,...,1.fq.part-i,2.fq.part-i,...,1.fq.part-n,2.fq.part-n, 
    threads           - number of threads
    """
    is_single = True
    # PE or SE    
    if len(inputs) == 2:
        is_single = False

    if is_single:
        splitter = [executables['FastqSplitter'],"--n-parts",str(split_times),"--prefix",prefix]
        
        if gz:
            splitter.append("--gz")
        
        splitter.append(inputs[0])
        
        tools = [splitter]        
        processOne = utils.run_tools(tools, input=None, output=None, name="split-fastq")
        
        if processOne.wait() != 0:
            raise ValueError("Error while  Fastq Splitter")
            
    else:
        splitterOne = [executables['FastqSplitter'],"--n-parts",str(split_times),"--prefix",prefix + ".1"]

        if gz:
            splitterOne.append("--gz")
            
        splitterOne.append(inputs[0])

        splitterTwo = [executables['FastqSplitter'],"--n-parts",str(split_times),"--prefix",prefix + ".2"]
        
        if gz:
            splitterTwo.append("--gz")

        splitterTwo.append(inputs[1])        
                
        toolsOne = [splitterOne]
        toolsTwo = [splitterTwo]
       
        if threads > 1:
            processOne = utils.run_tools(toolsOne, input=None, output=None, name="split-1-fastq")
            processTwo = utils.run_tools(toolsTwo, input=None, output=None, name="split-2-fastq")       
       
            if processOne.wait() != 0:
                raise ValueError("Error while Fastq Splitter Pair One")   

            if processTwo.wait() != 0:
                raise ValueError("Error while Fastq Splitter Pair Two")
        else:
            processOne = utils.run_tools(toolsOne, input=None, output=None, name="split-1-fastq")
            if processOne.wait() != 0:
                raise ValueError("Error while Fastq Splitter Pair One") 
                
            processTwo = utils.run_tools(toolsTwo, input=None, output=None, name="split-2-fastq")
            if processTwo.wait() != 0:
                raise ValueError("Error while Fastq Splitter Pair Two")


def chopFQ(inputs=None,output=None,first_position=10,
           split_times=100,kmer_len=36,slide=36,threads=1):
    """ Chop FASTQ file 
        inputs -- list of input fastq files 
              SE: fq.part-1,...,fq.part-i,...,fq.part-n,    
              PE: 1.fq.part-1,2.fq.part-1,...,1.fq.part-i,2.fq.part-i,...,1.fq.part-n,2.fq.part-n,    
        output -- list of output fastq files 
              SE: chop.part-1.fq,...,chop.part-i.fq,...,chop.part-n.fq    
              PE: chop.part-1.1.fq,chop.part-1.2.fq,...,chop.part-i.1.fq,chop.part-i.2.fq,...,chop.part-n.1.fq,chop.part-n.2.fq    
        first_position    - First position in the read to start the chopping process
        split_times       - Number of chunks
        kmer_len          - kmer length
        slide             - sliding, equalto to kmer_len for non overlapping windows
        threads           - number of threads
    """

     #List of task
    listInputOutput = []
    for (fragment,outputFile) in zip(inputs,output):
        listInputOutput.append([fragment,outputFile])
        
    #Scheduler of task
    listExecution = []    
    while len(listInputOutput) != 0:
        parallelCommands = []                
        for i in range(threads):
            if len(listInputOutput) != 0:
                parallelCommands.append(listInputOutput.pop())
        listExecution.append(parallelCommands)

    #Task launching
    for parallelTasks in listExecution:
        #Task MAKE BED INTERVAL
        processIntervals = []        
        for tasks in parallelTasks:
            chopping = [executables['kmermaker'],
                    '-i', tasks[0],
                    '-k', str(kmer_len),
                    '-s', str(slide)
                    ]      

            if first_position > 0:
                chopping.append('-f')
                chopping.append(str(first_position))
         
            tools = [chopping]
         
            processIntervals.append(utils.run_tools(tools, input=None, output=tasks[1], name="chop-fastq"))
            
            
        for process in processIntervals:
            if process.wait() != 0:
                raise ValueError("Error while running kmermaker, chopping fastq!!")
       



def cnvMapper(input,index,output=None,
            mismatch_alphabet="ACGT",
            mismatches=2,
            max_edit_distance=4,       
            strata_after_best=2,        
            max_decoded_matches=20,
            min_decoded_strata=1,
            threads_number=8   
        ):  
    """Start the GEM mapper on the given input.
    Mapping for the CNV analysis.
    
    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- A ReadIterator with the input
    output -- output file name
    index -- valid GEM2 index
    mismatch_alphabet -- Specifies the set of characters which are valid replacements in case of mismatch, default='ACGT'
    mismatches -- number or % mismatches, default=2
    max_edit_distance -- max edit distance, default=4 
    strata_after_best -- How many strata should be explored after the best one, default=2
    max_decoded_matches -- maximum decoded matches, defaults to 20
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 1
    
    threads_number -- Number of threads, default=8
    """ 
    
    ## prepare inputs
    index = _prepare_index_parameter(index)
        
    ## prepare the input
    pa = [executables['gem-mapper'], '-I', index,
          '-q', 'ignore',
          '--mismatch-alphabet', mismatch_alphabet,          
          '-m', str(mismatches),
          '-e', str(max_edit_distance),
          '-s', str(strata_after_best),
          '-d', str(max_decoded_matches),
          '-D', str(min_decoded_strata),
           '-T', str(threads_number)
    ]
    
    # build list of tools
    tools = [pa]
    
    ##compress mappings
    gzip = _compressor(threads=max(1, int(threads_number) / 2))
    tools.append(gzip)    
    
    pa.append("-i")
    pa.append(input)

    ## run the mapper
    process = utils.run_tools(tools, input=None, output=output, name="CNV-Mapper")
    return _prepare_output(process, output=output)

def mapToSam(input,index,output=None,name=None,threads=8,sort_memory="768M"):
    """Transforms MAP to SAM 

    input -- A ReadIterator with the input
    index -- valid GEM2 index
    output -- output file name
    name -- basic name ouput
    threads -- threads default 8    
    """ 
    ## prepare the input
    pa = [executables['gt-scorereads'],
           '--i1', input,
           '--gem-index', index,
           '-q','offset-33',
           '--output-format', 'SAM',
           '-t',str(threads)
    ]
        
    # build list of tools
    tools = [pa]
    
    #sam to bam 
    samToBam = ['samtools', 'view',
                '-Shb',
                '-F','4',
                '-@',str(threads),
                '-'
    ]
    
    tools.append(samToBam)
    
    #sort bam
    bamToSort = ['samtools', 'sort',
                '-@',str(threads),
                '-T','deleteme.' + name,
                '-O','sam',
    ]

    tools.append(bamToSort) 
   
    ##compress mappings
    gzip = _compressor(threads=max(1, int(threads) / 2))
    tools.append(gzip) 
   
    ## run the mapper
    process = utils.run_tools(tools, input=None, output=output, name="map-sam")
    return _prepare_output(process, output=output)
    
def mrCanavarReadDepth(input,conf_file = None, depth_file = None, gzipped = False):
    """mrCanavar compute Read Depth

    input -- A ReadIterator with the input
    conf_file -- Configuration file path
    output -- output file list
    depth_file -- depth file to store the results
    gzipped -- If true mappings are gzipped
    """        
    
    ## prepare input compute read depth
    readDepth = [executables['mrcanavar'],
          '--read','-conf',conf_file,
          '-samdir', input,
          '-depth', depth_file
    ]
    
    #GZIP mapping
    if gzipped:
        readDepth.append("--gz")
        
    # build list of tools
    tools = [readDepth]
    process = utils.run_tools(tools, input=None, output=None, name="cnv-read.depth")

    ## run mrcanavar
    if process.wait() != 0:
        raise ValueError("Mr Canavar Read Depth failed!")
        
    
    #process = utils.run_tools(tools, input=None, output=output, name="cnv-read.depth")
    #return _prepare_output(process, output=output)

def mrCanavarCalls(input,conf_file = None, calls_file = None):
    """mrCanavar perform copy number calls

    input -- A ReadIterator with the input
    conf_file -- Configuration file path
    output -- output file list
    """        
    
    ## prepare the input copy number call
    copyNumber = [executables['mrcanavar'],
          '--call','-conf',conf_file,
          '-depth', input,
          '-o', calls_file
    ]    
    
    tools = [copyNumber]

    ## run the mapper
    process = utils.run_tools(tools, input=None, output=calls_file, name="cnv-call")
    if process.wait() != 0:
        raise ValueError("Mr Canavar Read Call failed!")
    
def copyNumberDistribution(input,output = None, sampleName = None):
    """copyNumberDistribution Copy Number Analysis 

    input -- A ReadIterator with the input
    cn_RData_file -- Copy Number RData output   
    cut_offs_file -- Cut offs file R output
    plot_file -- Control Regions Density plot output file
    control_regions_distribution_file -- Control Regions Distribution output file
    """
    ## prepare inputs                               
    copyWindowsBed = _prepare_copywindows_bed(input)
    copyNumberBed = _prepare_copyNumber_bed(input)   
    
    cn_RData_file = output[0]
    cut_offs_file = output[1]
    plot_file = output[2]
    control_regions_distribution_file = output[3]
    
    ## prepare input R script
    cnRScript = [executables['copyNumberDistribution.R'],
                 '-sample=' + sampleName,
                 '-cw=' + copyWindowsBed,
                 '-cn=' + copyNumberBed,
                 '-rData=' + cn_RData_file,
                 '-cutOffs=' + cut_offs_file,
                 '-plotFile=' + plot_file,
                 '-distribution=' + control_regions_distribution_file
    ]
    
    # build list of tools
    tools = [cnRScript]
    
    process = utils.run_tools(tools, input=None, output=None, name="cn-distribution")
    if process.wait() != 0:
        raise ValueError("Copy number distribution analysis failed!")
    return output

def callDuplications(input, output = None, sampleName = None,bed_repeat_regions = None, 
                     bed_gaps_coordinates = None, copyNumberFile = None):
    """Call of Duplications
    
    input -- A ReadIterator with the input
    output -- output file list
    bed_repeat_regions -- Bed file with repeat regions coordinates (include TRF, Simple Repateas amd kmer masking regions
    bed_gaps_coordinates -- Bed file with gaps coordinates
    """

    duplication = gem.duplications.Duplications(sampleName,bedRepeatRegions=bed_repeat_regions,gapsBedCoord=bed_gaps_coordinates,
                                                cutOffFile=input[1][1], pathMrCanavar=input[0][0],outputlist=output)

    #First Method
    duplication.runMethod1Duplications()
    
    #Second Method
    duplication.runMethod2Duplications()
        

def transcript_mapper(input, indices, key_files, output=None,
           mismatches=0.04,
           delta=0,
           quality=33,
           quality_threshold=26,
           max_decoded_matches=100,
           min_decoded_strata=1,
           min_matched_bases=0.80,
           max_big_indel_length=15,
           max_edit_distance=0.20,
           mismatch_alphabet="ACGT",
           trim=None,
           threads=1,
           extra=None,
           ):
    """Start the GEM mapper on the given input.
    If input is a file handle, it is assumed to
    provide fastq entries. If input is a string,
    it is checked for its extension. In case of a
    .map file, the input is converted from gem format
    to fastq and passed to the mapper.

    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- A ReadIterator with the input
    output -- output file name
    index -- valid GEM2 index
    mismatches - number or % mismatches, default=0.04
    delta -- strata after best <number> (default=0)
    quality -- one of 'ignore'|'offset-33'|'offset-64' defaults to offset-33
    quality_threshold <number> -- (default=26, that is e<=2e-3)
    max_edit_distance -- max edit distance, 0.20 per default
    max_decoded_matches -- maximum decoded matches, defaults to 20
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 1
    min_matched_bases -- minimum number (or %) of matched bases, defaults to 0.80
    trim -- tuple or list that specifies left and right trimmings
    extra -- list of additional parameters added to gem mapper call
    """

    if not isinstance(indices, (list, tuple)):
        indices = [indices]

    if not isinstance(key_files, (list, tuple)):
        key_files = [key_files]

    outputs = []
    output_files = []
    for i, index in enumerate(indices):
        output_file = output
        if len(indices) > 1:
            (fifo, output_file) = tempfile.mkstemp(suffix=".map", prefix="transcript_mapping_output", dir=".")
            os.close(fifo)

        output_files.append(output_file)
        outputs.append(mapper(input.clone(), index,
           key_file=key_files[i],
           output=output_file,
           mismatches=mismatches,
           delta=delta,
           quality=quality,
           quality_threshold=quality_threshold,
           max_decoded_matches=max_decoded_matches,
           min_decoded_strata=min_decoded_strata,
           min_matched_bases=min_matched_bases,
           max_big_indel_length=max_big_indel_length,
           max_edit_distance=max_edit_distance,
           mismatch_alphabet=mismatch_alphabet,
           trim=trim,
           threads=threads,
           extra=extra,
           force_min_decoded_strata=True
           )
        )
    if len(indices) > 1:
        merged = merge(outputs[0], outputs[1:], paired=False, threads=threads,
                       same_content=True, output=output)
        if(output is not None):
            for f in output_files:
                os.remove(f)
            return _prepare_output(None, output=output, quality=quality)
        return merged
    else:
        return _prepare_output(None, output=output, quality=quality)


def splitmapper(input,
                index,
                output=None,
                mismatches=0.04,
                splice_consensus=extended_splice_consensus,
                filter=default_filter,
                refinement_step_size=2,
                min_split_size=15,
                matches_threshold=100,
                strata_after_first=1,
                mismatch_alphabet="ACGT",
                quality=33,
                trim=None,
                filter_splitmaps=True,
                post_validate=True,
                threads=1,
                extra=None):
    """Start the GEM split mapper on the given input.
    If input is a file handle, it is assumed to
    provide fastq entries. If input is a string,
    it is checked for its extension. In case of a
    .map file, the input is converted from gem format
    to fastq and passed to the mapper.

    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- string with the input file or a file handle or a generator
    output -- output file name or file handle
    index -- valid GEM2 index
    """

    ## check the index
    index = _prepare_index_parameter(index, gem_suffix=True)
    if quality is None and isinstance(input, files.ReadIterator):
        quality = input.quality
    quality = _prepare_quality_parameter(quality)
    splice_cons = _prepare_splice_consensus_parameter(splice_consensus)

    pa = [executables['gem-rna-tools'],
          'split-mapper',
          '-I', index,
          '-q', quality,
          '-m', str(mismatches),
          '--min-split-size', str(min_split_size),
          '--refinement-step-size', str(refinement_step_size),
          '--matches-threshold', str(matches_threshold),
          '-s', str(strata_after_first),
          '--mismatch-alphabet', mismatch_alphabet,
          '-T', str(threads)
    ]
    min_threads = int(round(max(1, threads / 2)))

    if filter is not None:
        pa.append("-f")
        pa.append(filter)
    if splice_cons is not None:
        pa.append("-c")
        pa.append(splice_cons)

    ## extend with additional parameters
    _extend_parameters(pa, extra)
    trim_c = [executables['gem-2-gem'], '-c', '-T', str(min_threads)]
    if trim is not None:
        ## check type
        if not isinstance(trim, (list, tuple)) or len(trim) != 2:
            raise ValueError("Trim parameter has to be a list or a tuple of size 2")
        input = gemfilter.trim(input, trim[0], trim[1], append_label=True)

    tools = [pa]
    if filter_splitmaps:
        tools.append(__awk_filter)
    if trim is not None:
        tools.append(trim_c)

    ## run the mapper
    process = None
    original_output = output
    if post_validate:
        output = None

    raw = False
    if isinstance(input, gt.InputFile) and input.raw_sequence_stream():
        raw = False
        pa.append("-i")
        pa.append(input.filename)
        input = None

    process = utils.run_tools(tools, input=input, output=output, name="GEM-Split-Mapper", raw=raw)
    splitmap_out = _prepare_output(process, output=output, quality=quality)

    if post_validate:
        return validate(splitmap_out, index, original_output, threads=threads)

    return splitmap_out


def extract_junctions(input,
                      index,
                      filter="ordered,non-zero-distance",
                      mismatches=0.04,
                      refinement_step_size=2,
                      min_split_size=15,
                      matches_threshold=75,
                      splice_consensus=extended_splice_consensus,
                      strata_after_first=1,
                      quality=33,
                      threads=1,
                      merge_with=None,
                      min_split=4,
                      max_split=2500000,
                      coverage=0,
                      max_junction_matches=5,
                      tmpdir=None,
                      annotation=None,
                      extra=None):
    ## run the splitmapper
    splitmap = splitmapper(input,
        index,
        output=None,
        filter=filter,
        mismatches=mismatches,
        refinement_step_size=refinement_step_size,
        min_split_size=min_split_size,
        matches_threshold=matches_threshold,
        splice_consensus=splice_consensus,
        quality=quality,
        strata_after_first=strata_after_first,
        filter_splitmaps=False,
        post_validate=False,
        threads=threads,
        extra=extra)

    annotation_junctions = None
    if annotation is not None:
        annotation_junctions = annotation
        if isinstance(annotation, basestring):
            annotation_junctions = gem.junctions.from_gtf(annotation)

    denovo_junctions = splits.extract_denovo_junctions(
        splitmap.raw_stream(),  # pass the raw stream
        minsplit=min_split,
        maxsplit=max_split,
        coverage=coverage,
        sites=merge_with,
        max_junction_matches=max_junction_matches,
        process=splitmap.process,
        threads=max(1, threads / 2),
        annotation_junctions=annotation_junctions
    )
    return denovo_junctions


def pairalign(input, index, output=None,
              quality=33,
              quality_threshold=26,
              max_decoded_matches=20,
              min_decoded_strata=1,
              min_insert_size=0,
              max_insert_size=1000,
              max_edit_distance=0.30,
              min_matched_bases=0.80,
              max_extendable_matches=0,
              max_matches_per_extension=0,
              unique_pairing=False,
              map_both_ends=False,
              filter_max_matches=0,
              threads=1,
              compress=False,
              extra=None):
    ## check the index
    index = _prepare_index_parameter(index)
    quality = _prepare_quality_parameter(quality, input)
    if compress and output is None:
        logging.warning("Disabeling stream compression")
        compress = False

    if compress and not output.endswith(".gz"):
        output += ".gz"

    pa = [executables['gem-mapper'],
          '-p',
          '-I', index,
          '-q', quality,
          '--gem-quality-threshold', str(quality_threshold),
          '--max-decoded-matches', str(max_decoded_matches),
          '--min-decoded-strata', str(min_decoded_strata),
          '--min-insert-size', str(min_insert_size),
          '--max-insert-size', str(max_insert_size),
          '-E', str(max_edit_distance),
          '--min-matched-bases', str(min_matched_bases),
          '--max-extendable-matches', str(max_extendable_matches),
          '--max-extensions-per-match', str(max_matches_per_extension),
          '-T', str(threads)
    ]

    ## extend with additional parameters
    _extend_parameters(pa, extra)

    if unique_pairing:
        pa.append("--unique-pairing")
    if map_both_ends:
        pa.append("--map-both-ends")

    # if qualities are ignored, make sure we remove the quality column
    # otherwise the pairaligner will crash
    tools = []
    if quality in ['ignore', 'none']:
        tools.append(__awk_pair_quality_fix)
    tools.append(pa)

    filter_pa = [executables["gt.filter"], "-t", str(threads), "-p"]
    if filter_max_matches > 0:
        filter_pa.extend(["--max-output-matches", str(filter_max_matches)])
    tools.append(filter_pa)
    if compress:
        gzip = _compressor(threads=max(1, threads / 2))
        tools.append(gzip)

    raw = False
    if isinstance(input, gt.InputFile):
        raw = True

    ## run the mapper and trim away all the unused stuff from the ids
    process = utils.run_tools(tools, input=input, output=output, name="GEM-Pair-align", write_map=True, clean_id=True, append_extra=False, raw=raw)
    return _prepare_output(process, output=output, quality=quality)


def realign(input,
            index,
            output=None,
            threads=1, ):
    index = _prepare_index_parameter(index, gem_suffix=True)
    validate_p = [executables['gem-2-gem'],
                  '-I', index,
                  '-r',
                  '-T', str(threads)
    ]
    process = utils.run_tool(validate_p, input=input, output=output, name="GEM-Realign", write_map=True)
    return _prepare_output(process, output=output)


def validate(input,
             index,
             output=None,
             validate_score=None,  # "-s,-b,-i"
             validate_filter=None,  # "1,2,25"
             threads=1, ):
    index = _prepare_index_parameter(index, gem_suffix=True)
    validate_p = [executables['gem-2-gem'],
                  '-I', index,
                  '-v', '-r',
                  '-T', str(max(threads, 1))
    ]
    if validate_score is not None:
        validate_p.extend(["-s", validate_score])
    if validate_filter is not None:
        validate_p.extend(['-f', validate_filter])

    process = utils.run_tool(validate_p, input=input, output=output, name="GEM-Validate", write_map=True, raw=True)
    return _prepare_output(process, output=output)


def score(input,
          index,
          output=None,
          scoring="+U,+u,-s,-t,+1,-i,-a",
          filter=None,  # "1,2,25"
          quality=None,
          compress=False,
          threads=1,
          raw=False,
          remove_existing=False):
    """Score the input. In addition, you can specify a tuple with (<score_strata_to_keep>,<max_strata_distance>,<max_alignments>) to
    filter the result further.
    """
    if compress and output is None:
        logging.warning("Disabeling stream compression")
        compress = False

    if compress and not output.endswith(".gz"):
        output += ".gz"

    quality = _prepare_quality_parameter(quality)
    if quality in ['none', 'ignore']:
        quality = 'offset-33'
    index = _prepare_index_parameter(index, gem_suffix=True)
    score_p = [executables['gem-2-gem'],
               '-I', index,
               '-q', quality,
               '-s', scoring,
               '-T', str(threads)
    ]

    if filter is not None:
        score_p.append("-f")
        ff = filter
        if not isinstance(filter, basestring):
            ff = ",".join([str(f) for f in filter])
        score_p.append(ff)

    if raw or isinstance(input, gt.InputFile):
        raw = True
        if isinstance(input, gt.InputFile) and remove_existing:
            input.remove_scores = True
            raw = False
        #input = input.raw_stream()

    tools = [score_p]

    if compress:
        gzip = _compressor(threads=threads)
        tools.append(gzip)

    process = utils.run_tools(tools, input=input, output=output, name="GEM-Score", write_map=True, raw=raw)
    return _prepare_output(process, output=output)


def gem2sam(input, index=None, output=None,
            single_end=False, compact=False, threads=1,
            quality=None, check_ids=True, add_length=True, consensus=None,
            exclude_header=False, calc_xs=True, raw=False):

    if index is not None:
        index = _prepare_index_parameter(index, gem_suffix=True)

    gem_2_sam_p = [executables['gem-2-sam'], '-T', str(threads)]
    if index is not None:
        gem_2_sam_p.extend(['-I', index])
        if not exclude_header:
            gem_2_sam_p.append("-l")

    tools = []
    quality = _prepare_quality_parameter(quality, input)
    if quality is not None and not quality == "ignore":
        gem_2_sam_p.extend(["-q", quality])
    else:
        # apply fix
        gem_2_sam_p.extend(["-q", 'offset-33'])
        tools.append(__awk_gem_2_sam_quality_fix)
        raw = True

    if consensus is not None and calc_xs:
        gem_2_sam_p.extend([
            '-s', _prepare_splice_consensus_parameter(consensus)
        ])

    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")

    tools.append(gem_2_sam_p)
    
    # GT-25 transform id's
    process = utils.run_tools(tools, input=input, output=output,
                              name="GEM-2-sam", write_map=True,
                              clean_id=not single_end,
                              append_extra=False, raw=raw)
    return _prepare_output(process, output=output, quality=quality)


def _check_samtools(command, threads=1, extend=None):
    """Return based on threads and existence of parallel samtools"""
    global __parallel_samtools
    if threads == 1:
        p = [executables["samtools"], command]
    else:
        if __parallel_samtools is None:
            (s, e) = subprocess.Popen([executables["samtools"], "view"], stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
            f = __builtin__.filter(lambda x: x.startswith("-@"), [l.strip() for l in e.split("\n")])
            __parallel_samtools = len(f) > 0

        if __parallel_samtools:
            p = [executables["samtools"], command, "-@", str(threads)]
        else:
            p = [executables["samtools"], command]
    if extend is not None:
        p.extend(extend)
    return p

def __is_parallel_samtools():
    if __parallel_samtools is None:
        _check_samtools("view", threads=2)
    return __parallel_samtools

def sam2bam(input, output=None, sorted=False, tmpdir=None, mapq=None, threads=1, sort_memory="768M"):
    sam2bam_p = _check_samtools("view", threads=threads, extend=["-S", "-b"])
    if mapq is not None and int(mapq) > 0:
        sam2bam_p.append("-q")
        sam2bam_p.append(str(mapq))
    sam2bam_p.append('-')

    tools = [sam2bam_p]
    out_name = output
    if sorted:
        if not __is_parallel_samtools():
            # check the memory paramters
            try:
                m = int(sort_memory)
                if m < 128 * 1024 * 128:  # ugly but we assume you give it at least 128 mb
                    sort_memory = 768 * 1024 * 1024
                    if m < 1024 * 32:
                        sort_memory = m * 1024 * 1024
            except Exception, e:
                # convert to default byte
                sort_memory = 768 * 1024 * 1024

        bam_sort = _check_samtools("sort", threads=threads, extend=["-m", str(sort_memory), "-o", "-"])
        suffix = ""
        if output is not None:
            suffix = "-" + os.path.basename(output)
        tmpfile = tempfile.NamedTemporaryFile(prefix="sort", suffix=suffix)
        tmpfile.close()
        if os.path.exists(tmpfile.name):
            os.remove(tmpfile.name)
        out_name = os.path.basename(tmpfile.name)
        bam_sort.append(out_name)
        tools.append(bam_sort)

    process = utils.run_tools(tools, input=input, output=output, name="SAM-2-BAM", raw=True)
    return _prepare_output(process, output=output, quality=33, bam=True)


def bamIndex(input, output=None):
    in_name = input
    if not isinstance(input, basestring):
        in_name = input.filename
    bam_p = ['samtools', 'index', in_name]
    if output is None:
        output = in_name + ".bai"
    bam_p.append(output)
    tools = [bam_p]
    process = utils.run_tools(tools, input=None, output=None, name="BAM-Index")
    if process.wait() != 0:
        raise ValueError("BAM indexing failed!")
    return output


def compute_transcriptome(max_read_length, index, junctions, substract=None, output_name=None):
    """Compute the transcriptome based on a set of junctions. You can optionally specify
    a *substract* junction set. In that case only junctions not in substract are passed
    to compute the transcriptome.
    The function returns a tuple of a .fa file with the transcriptome genome and a
    .keys file with the translation table.

    max_read_length -- the maximum read length
    index -- path to the gem index
    junctions -- path to the junctions file
    substract -- additional juntions that are not taken into account and substracted from the main junctions
    """
    if output_name is None:
        output_name = os.path.abspath(junctions)
    transcriptome_p = [
        executables['gem-rna-tools'],
        'compute-transcriptome',
        '-l', str(max_read_length),
        '-I', index,
        '-j', junctions,
        '-o', output_name
    ]
    if substract is not None:
        transcriptome_p.extend(['-J', substract])

    process = utils.run_tools([transcriptome_p], input=None, output=None, name="compute-transcriptome")
    if process.wait() != 0:
        raise ValueError("Error while computing transcriptome")

    return (os.path.abspath("%s.fa" % junctions), os.path.abspath("%s.keys" % junctions))


def index(input, output, content="dna", threads=1):
    """Run the gem-indexer on the given input. Input has to be the path
    to a single fasta file that contains the genome to be indexed.
    Output should be the path to the target index file. Note that
    the gem index has to end in .gem and the prefix is added if necessary and
    the returned path will always be the correct path to the index.

    The method checks for the existence of the target index file
    and will NOT override but exit silently without recreating the index!

    Returns the absolute path to the resulting index file
    """
    indexer_p = [
        executables['gem-indexer'],
        '-T', str(threads),
        '--content-type', content.lower()
    ]

    if isinstance(input, basestring):
        if not os.path.exists(input):
            raise ValueError("Indexer input file %s not found" % input)
        indexer_p.extend(["-i", input])
    else:
        raise ValueError("The indexer wrapper can not handle the input %s, pass a file or a list of files" % input)

    existing = output
    if existing[-4:] != ".gem": existing = "%s.gem" % existing
    if os.path.exists(existing):
        logging.warning("Index %s already exists, skipping indexing" % existing)
        return os.path.abspath(existing)

    # indexer takes the prefix
    if output[-4:] == ".gem":
        output = output[:-4]
    indexer_p.extend(['-o', output])

    # the indexer need the other indexer tools in PATH
    path = "%s:%s" % (os.path.dirname(executables['gem-indexer']), os.getenv("PATH"))

    process = utils.run_tools([indexer_p], name="gem-indexer", env={"PATH": path})
    if process.wait() != 0:
        raise ValueError("Error while executing the gem-indexer")
    return os.path.abspath("%s.gem" % output)


def prep(input, output, gaps=None,lw_size=None,lw_slide=None,sw_size=None,\
         sw_slide=None,cw_size=None,pseudoa=None):
    """Run PREP step of mrCanavar
       input - FASTA reference file
       output - Congirutation output file
       gaps - Gap Bed file annotation
       lw_size - Long window span size
       lw_slide - Long window slide size
       sw_size - Short window span size
       sw_slide - Short window slide size
       cw_size - Copy number window size
       pseudoa - Coordinates for pseudoautosomal regions in the reference genome in BED format
    """
    canavar = [
        executables['mrcanavar'],'--prep',
        '-fasta',input,
        '-gaps',gaps,
        '-conf',output,
        '-lw_size',str(lw_size),
        '-lw_slide',str(lw_slide),
        '-sw_size',str(sw_size),
        '-sw_slide',str(sw_slide),
        '-cw_size',str(cw_size)
    ]
    
    if pseudoa:
        canavar.append('-pseudoa',pseudoa)
        
   
    process = utils.run_tools([canavar], name="gem-prep",output=input + ".log")
    if process.wait() != 0:
        raise ValueError("Error while executing the canavar PREP")
	 
	



def gtfcounts(inputs, annotation, output=None, json_output=None, threads=1,
              counts=None, weight=True, multimaps=False, exon_threshold=0,
              paired=False, coverage=False):
    """Run the count stats. This returns the gtf count stats as dictionary"""
    p = [
        executables['gt.gtfcount'],
        '--threads', str(threads),
        '-f', 'both',
        '-a', annotation
    ]
    if paired:
        p.append("-p")

    if counts is not None:
        p.extend(["-g", counts])
        if weight:
            p.append("-w")
        if multimaps:
            p.append("-m")
        if exon_threshold > 0:
            p.extend(["-e", str(exon_threshold)])
    if coverage:
        p.append('-c')

    from subprocess import PIPE
    process = utils.run_tools([p], name="gtfcounts", input=inputs,
                              output=output, raw=True, write_map=True,
                              logfile=PIPE)

    if json_output is not None:
        json_output = open(json_output, 'w')

    lines = []
    for line in process.processes[0].process.stderr:
        lines.append(line.strip())
        if json_output is not None:
            json_output.write(line)

    if json_output is not None:
        json_output.close()

    if process.wait() != 0:
        raise ValueError("Error while running gt.gtfcounts")
    import json
    return json.loads("\n".join(lines))


def stats(inputs, output=None, json_output=None, threads=1, paired=False, raw=True):
    """Run the count gt.stats with all options enabled. The function returns the parsed
    json stats.
    """
    p = [
        executables['gt.stats'],
        '-a',
        '-t', str(threads),
        '-f', 'both'
    ]
    if paired:
        p.append("-p")
        

    from subprocess import PIPE
    process = utils.run_tools([p], name="stats", input=inputs,
                              output=output, raw=raw, write_map=True,
                              force_debug=True, logfile=PIPE)
    if json_output is not None:
        json_output = open(json_output, 'w')

    lines = []
    for line in process.processes[0].process.stderr:
        lines.append(line.strip())
        if json_output is not None:
            json_output.write(line)

    if json_output is not None:
        json_output.close()

    if process.wait() != 0:
        raise ValueError("Error while running gt.stats")
    import json
    return json.loads("\n".join(lines))

def statsCnv(inputs, output=None, json_output=None, threads=1, paired=False, raw=True):
    """Run the count gt.stats with all options enabled. The function returns the parsed
    json stats.
    """
    p = [
        executables['gt.stats'],
        '-a',
        '-i',inputs,
        '-t', str(threads),
        '-f', 'both'
    ]
    if paired:
        p.append("-p")
        

    from subprocess import PIPE
    process = utils.run_tools([p], name="stats", input=None,
                              output=output, raw=raw, write_map=True,
                              force_debug=True, logfile=PIPE)
    if json_output is not None:
        json_output = open(json_output, 'w')

    lines = []
    for line in process.processes[0].process.stderr:
        lines.append(line.strip())
        if json_output is not None:
            json_output.write(line)

    if json_output is not None:
        json_output.close()

    if process.wait() != 0:
        raise ValueError("Error while running gt.stats")
    import json
    return json.loads("\n".join(lines))

def _compressor(threads=1):
    """Returns compressor configuration
    for compressing streams"""
    pigz = gem.utils.which("pigz")
    if threads == 1 or pigz is None:
        return ["gzip", "-"]
    else:
        return [pigz, "-p", str(threads), "-"]


def merge(master, slaves, output=None, paired=False, same_content=False,
          compress=False, threads=1):
    """Merge the content of the master with the content
    of the salve(s).

    """
    merge_out = subprocess.PIPE
    if output is not None:
        if output == sys.stdout:
            merge_out = output
            output = None
        else:
            if not compress:
                merge_out = output
                if isinstance(output, basestring):
                    merge_out = open(output, 'wb')
            else:
                p = subprocess.Popen(_compressor(threads=max(1, threads / 2)),
                                     stdout=open(output, 'wb'),
                                     stdin=subprocess.PIPE, close_fds=True)
                merge_out = p.stdin

    # create tmpdir for the fifos
    tmpdir = tempfile.mkdtemp()
    gem.files.delete_on_exit.append(tmpdir)
    current_master = master
    current_process = None
    for i, slave in enumerate(slaves):
        current_output = None
        if i == (len(slaves) - 1):
            # last one
            current_output = merge_out
        current_process = _merge_two(current_master, slave, current_output,
                                     tmpdir, i, paired, same_content, threads)
        current_master = current_process.stdout
    return _prepare_output(current_process, output=output)


def _merge_two(master, slave, output, tmpdir, count,
               paired=False, same_content=False, threads=1):
    """Helper function to merge to files. Return the merging
    process.
    """
    pa = [executables['gt.mapset'], '-C', 'merge-map', '-t', str(threads)]
    if paired:
        pa.append('-p')
    if same_content:
        pa.append('-s')


    inmaster = master
    if isinstance(inmaster, basestring):  # from string
        inmaster = gem.files.open_file(inmaster)
    elif hasattr(inmaster, 'stdout'):  # from process
        inmaster = inmaster.stdout
    elif isinstance(inmaster, gt.InputFile):  # from gt input file
        inmaster = inmaster.raw_stream()

    inslave = slave
    fifo_in = None
    fifo_out = None
    if not isinstance(inslave, basestring):
        # create a fifo and
        filename = os.path.join(tmpdir, "%d" % count)
        os.mkfifo(filename)
        if not hasattr(inslave, 'stdout'):
            inslave = inslave.raw_stream()
        fifo_in = inslave
        inslave = filename

    if output is None:
        output = subprocess.PIPE
    elif isinstance(output, basestring):
        output = open(output, 'wb')

    pa.append('--i1')
    pa.append(inslave)

    p = subprocess.Popen(pa, stdin=inmaster, stdout=output)
    if fifo_in is not None:
        fifo_out = open(filename, 'wb')
        subprocess.Popen(['cat'], stdin=fifo_in, stdout=fifo_out)
    return p


def maskFastaFromBed(input=None,output=None,regions=None):
    """Mask a fasta file according to a set of bed file
    regions

    input -- A raw fasta input file
    output -- Bed Regions output and Masked output in a list
    regions -- list of files to be masked
    """
    #Get all bed regions
    allRegions = ["cat"]
    for bedFile in regions:
        allRegions.append(bedFile)
        
    #build list of tools
    tools = [allRegions]
    
    #Get just coordinates fields
    coord = ["awk", "-F", "\t", '{print $1"\t"$2"\t"$3}']
    
    tools.append(coord)

    #Sort bed regions
    sortCoord = ["sortBed"]
    
    tools.append(sortCoord)
    
    #Merge Regions
    mergeCoord = ["mergeBed"]
    
    tools.append(mergeCoord)
    
    processOne = utils.run_tools(tools, input=None, output=output[0], name="bed-basic-mask")

    if processOne.wait() != 0:
        raise ValueError("Error while running bed-basic-mask")
    
    #maskFastaFromBed
    maskFasta = ["maskFastaFromBed",
                 "-fi",input,
                 "-fo",output[1],
                 "-bed",output[0]
    ]
    
    tools = [maskFasta]
    
    processTwo = utils.run_tools(tools, input=None, output=None, name="mask-fasta-bed")

    if processTwo.wait() != 0:
        raise ValueError("Error while running mask-fasta-bed")
    return _prepare_output(processTwo, output=None)


def buildAssemblyKmers(input=None,output=None,reference=None,threads=1):
    """ Create assembly fragements of 36 bp
    input - dictionary of chromsome name and its lengths
    output - list of output files per chromosome and length 
    reference - reference to store the data    
    threads - number of threads
    """
    
    #Get Path Out
    pathOut = os.path.dirname(output[0]) 

    #List of task
    listChrLenOut = []
    for (chrom,outputFa) in zip(sorted(input.keys()),output):
        listChrLenOut.append([chrom,input[chrom],outputFa])

    #Scheduler of task
    listExecution = []    
    while len(listChrLenOut) != 0:
        parallelCommands = []                
        for i in range(threads):
            if len(listChrLenOut) != 0:
                parallelCommands.append(listChrLenOut.pop())
        listExecution.append(parallelCommands)

    #Task launching
    for parallelTasks in listExecution:
        #Task MAKE BED INTERVAL
        processIntervals = []        
        for tasks in parallelTasks:
            chrOutBed = pathOut + "/" + tasks[0]  + ".bed"
            makeBedInterval = [executables['makeBedIntervalsAssembly'],
                           "-chrName",tasks[0],
                           "-chrLen",tasks[1],
                           "-outFile",chrOutBed
            ]
    
            tools = [makeBedInterval]
            processIntervals.append(utils.run_tools(tools, input=None, output=None, name="ASS-KMERS-BED"))
            
        for process in processIntervals:
            if process.wait() != 0:
                raise ValueError("Error while running make bed intervals")
        #Task Fasta from bed
        processFastaBed = []
        for tasks in parallelTasks:
            chrOutBed = pathOut + "/" + tasks[0]  + ".bed"
            fastaBed = ["fastaFromBed",
                    "-fi",reference,
                    "-bed",chrOutBed,
                    "-fo", tasks[2]
            ]
            tools = [fastaBed]
            processFastaBed.append(utils.run_tools(tools, input=None, output=None, name="ASS-KMERS-FASTA"))

        for process in processFastaBed:
            if process.wait() != 0:
                raise ValueError("Error while running fasta from bed")

        
        
def mapKmerWindows(input=None,output=None,reference=None,threads=1):
    """ Map each chromosome fragment against itself
    input - list of files with the sequences to be mapped
    output - list of output files chromosome  
    reference - fasta reference to map the set of reads """
    
    #Check index
    indx = is_gem_index(reference)
    if indx == "":
        file_name = ""
        if reference.endswith(".fa"):
            file_name += reference.replace(".fa","")
        else:
            file_name += reference.replace(".fasta","") 
        indx += index(reference, file_name, content="dna", threads=threads)
        

    ## Mapping
    for (kmer,maps) in zip(input,output):
        kmerMap = [executables['gem-mapper'], '-I', indx,
          '-q', 'ignore',
          '-m', '0',
          '-e', '0',
          '-d', '0',
          '-D', '0',
          '-T', str(threads),
          '-i', kmer
        ]
        tools = [kmerMap]
        
        ## run the mapper
        process = utils.run_tools(tools, input=None, output=maps, name="kmer-map")
        if process.wait() != 0:
            raise ValueError("Error while running map kmer windows")
            
def kmerCounts(input=None,output=None):
    """ Calls kmerCounts script, from a gem map file gets those reads without N on its sequence
        and more than two hits 
        input - list of input map files
        output - list of bed output results """
    
    for (mapFile,bedFile) in zip(input,output):
        counter = [executables['kmerCount'],
                   mapFile,
                   bedFile
                   ]
                   
        tools = [counter]
        
        ## Run kmer Counting 
        process = utils.run_tools(tools, input=None, output=None, name="kmer-Counting")
        if process.wait() != 0:
            raise ValueError("Error while running kmer counts")
    
def distributionPlot(input=None,output=None,hits=None):
    """Builds a distribution plot and creates a files of kmer regions over a given thershold 
    of reads
    input - list of kmer count files per chromosome
    output - list of output files (global mapping kmer counts, accumulative distribution plot, global kmer overrepresentation)
    hits - threshold of reads to output overrepresentated reads """

    mergeKmerCounts = ["cat"]    
    for kmerCountFile in input:
        mergeKmerCounts.append(kmerCountFile)
        
    tools = [mergeKmerCounts]
    
    #Run make histogram accumulative distribution
    processOne = utils.run_tools(tools, input=None, output=output[0], name="distribution-plot-merging")
    if processOne.wait() != 0:
        raise ValueError("Error while running distribution plot, merging step")
    
    buildHistogram = [executables['makeHistogramCounts.R'],
                      "--countFile=" + output[0], 
                      "--pngPlot=" + output[1], 
                      "--histoPlot=" + output[2],
                      "--kmerOverThreshold=" + output[3], 
                      "--threshold=" + str(hits) ]
                      
    tools = [buildHistogram]
    
    #Run make histogram accumulative distribution
    processTwo = utils.run_tools(tools, input=None, output=None, name="distribution-plot")
    if processTwo.wait() != 0:
        raise ValueError("Error while running distribution plot")
        
        
def maskReference(basic_reference=None,kmer_over=None,bed_regions=None,kmer_masked_fasta=None,threads=1):
    """Mask Basic Reference with over represented kmers
    basic_reference = basic masked reference    
    kmer_over = bed file with overrepresented kmers
    bed_regions = (output) over represented regions
    kmer_masked_fasta = (output) fasta reference file kmer masked
    """
    #Correction coordinates    
    correction = ["awk", "-F", "\t", '{start=$2-1; if(start < 0){start=$2;} end=$2+36; print $1"\t"start"\t"end}']
    correction.append(kmer_over)
    tools = [correction]        
        
    #coordinates sorting    
    #sort = ["sortBed"]
    #tools.append(sort)    
    
    #Merge Regions
    #merge = ["mergeBed"]
    #tools.append(merge)
    
    #Run Create merging bed
    processOne = utils.run_tools(tools, input=None, output=bed_regions, name="mask-reference")
    if processOne.wait() != 0:
        raise ValueError("Error while running masking reference merging regions")
        
    masking = ["maskFastaFromBed",
               "-fi", basic_reference, 
               "-bed", bed_regions,
               "-fo", kmer_masked_fasta
               ]
               
    tools = [masking]    
    
    #Run Masking Process
    processTwo = utils.run_tools(tools, input=None, output=None, name="mask-reference")
    if processTwo.wait() != 0:
        raise ValueError("Error while running masking reference maskFastaFromBed")
        
    #Check index
    indx = is_gem_index(kmer_masked_fasta)
    if indx == "":
        file_name = ""
        if kmer_masked_fasta.endswith(".fa"):
            file_name += kmer_masked_fasta.replace(".fa","")
        else:
            file_name += kmer_masked_fasta.replace(".fasta","") 
        indx += index(kmer_masked_fasta, file_name, content="dna", threads=threads)
        
    

def pad36Fasta(basic_masked_regions=None,kmer_regions=None,masked_fasta=None,chr_len=None,pad_regions=None,chrom_info_bed=None,pad_fasta=None):
    """Pad 36 bp over repeats, gaps and kmer regions
    basic_masked_regions - list of bed files for the basic regions masked    
    kmer_regions - kmer regions bed file
    masked_fasta - Kmer Masked Reference
    chr_len - dictionary of chromsome lengths
    pad_regions - (output) regions padded with 36 bp
    chrom_info_bed - chrom length file in bed format
    pad_fasta - (output) path to the fasta to be paddd masked
    """
    #Create bed chrom info
    with open(chrom_info_bed, 'w') as bedFile:
        for chrom in sorted(chr_len.keys()):
            bedFile.write(chrom + '\t0\t' + chr_len[chrom] + '\n')    
    
    #Get all maskeable regions
    merge = ["cat",
             kmer_regions]
    for bedFile in basic_masked_regions:
        merge.append(bedFile)
    tools = [merge]    
    
    #Bed Coordinates
    bedCoord = ["awk",'{print $1"\t"$2"\t"$3}']     
    tools.append(bedCoord)   
    
    #Sort them
    #sort = ["sortBed"]
    #tools.append(sort) 
    
    #Merge Regions
    #merge = ["mergeBed"]
    #tools.append(merge)
    
    #Padding
    padding = ["awk", '{start=$2; if(start-36>0){start=start-36;} end=$3+36; print $1"\t"start"\t"end}']
    tools.append(padding)  
    
    #Intersection chromosome limits
    intersectionLen = ["intersectBed",
                     "-a","stdin",
                     "-b",chrom_info_bed]  
    tools.append(intersectionLen)
   
    #Run Create merging bed
    processTwo = utils.run_tools(tools, input=None, output=pad_regions, name="pad36-fasta")
    if processTwo.wait() != 0:
        raise ValueError("Error while running pad 36 masking regions")
      
    #Masking Pad 36            
    masking = ["maskFastaFromBed",
               "-fi", masked_fasta, 
               "-bed", pad_regions,
               "-fo", pad_fasta
               ]
               
    tools = [masking]    
    
    #Run Masking Process
    processThree = utils.run_tools(tools, input=None, output=None, name="pad36-fasta-mask-reference")
    if processThree.wait() != 0:
        raise ValueError("Error while running masking reference maskFastaFromBed pad 36")

def galculator(input=None,output=None):
    """ Calculates Nuleotide distribution over a fasta file 
    input - Path to the fasta file to be processed    
    output - File to store galculator results
    """
    galculator = [executables["galculator"],input]  
    tools = [galculator]

    #Run Create merging bed
    processTwo = utils.run_tools(tools, input=None, output=output, name="galculator")
    if processTwo.wait() != 0:
        raise ValueError("Error while running galculator")
        

def _is_i3_compliant(stream):
    """Reads lines from the input stream and scans "flags" lines
    and returns true if the flags are compatible with the GEM
    i3 bundle. This is usually filled with the content of /proc/cpuinfo
    to determine the current systems capabilities

    The input is a stream that must provide a readline method"""
    i3_flags = set(["popcnt", "ssse3", "sse4_1", "sse4_2"])
    cpu_flags = set([])
    for line in iter(stream.readline, ''):
        line = line.rstrip()
        if line.startswith("flags"):
            for e in line.split(":")[1].strip().split(" "):
                cpu_flags.add(e)
    return i3_flags.issubset(cpu_flags)


def _extend_parameters(pa, extra=None):
    """Extend parameter array pa
    with extra parameters. If extra is a string, it is
    split by space and parameters are added separatly,
    otherwise it is assumed that extra is a list and
    it is appended to the parameter array as is.
    """
    if extra is not None:
        if isinstance(extra, (basestring,)):
            pa.extend(extra.split(" "))
        else:
            pa.extend(extra)
