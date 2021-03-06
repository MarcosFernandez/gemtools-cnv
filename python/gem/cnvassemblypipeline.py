#!/usr/bin/env
"""Pipeline utilities"""
import gem
import json
import logging
import os
import errno
import signal
import traceback
import re

from gem.utils import Timer
import gem.gemtools as gt
import gem.filter
import gem.utils
import gem.assemblyReport
from gem.pipeline import *
from types import *

class BasicMaskingStep(PipelineStep):
    """Basic Masking step"""
    
    def files(self):
        """Return the output files generated by this step.
        By default one .map output file is generated
        """
        if self._files is None:
            self._files = []

            bedRegions = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_basic_mask",file_suffix="bed")
            basicMask = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_basic_mask",file_suffix=self.file_suffix)
            
            self._files.append(bedRegions)
            self._files.append(basicMask)
                
        return self._files

    def run(self):
        cfg = self.configuration
        gem.maskFastaFromBed(input=cfg["raw_fasta"],output=self.files(),regions=cfg["bed_files"])

class GetAssemblyKmers(PipelineStep):
    """Fastq Fragmentation """
    def files(self):
        """Return the output files generated by this step.
        By default one .map output file is generated
        """        
        if self._files is None:
            self._files = []
            for chrom in sorted(self.pipeline.chrLen.keys()):
                chromFragments = "_" + chrom + "_kmers_36"
                file = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix=chromFragments,file_suffix=self.file_suffix)
                self._files.append(file)
            
        return self._files
            
        
    def run(self):
        cfg = self.configuration
        gem.buildAssemblyKmers(input=self.pipeline.chrLen,output=self.files(),reference=self.getAllFiles()[0][1],threads=int(cfg["threads"])) 
        
    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on basic masking step!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0]
        
class MapKmerWindows(PipelineStep):
    """Map kmer fragments to the masked reference"""
    def files(self):
        """Return the output files generated by this step.
        By default one .map output file is generated
        """        
        if self._files is None:
            self._files = []
            for chrom in sorted(self.pipeline.chrLen.keys()):
                chromFragments = "_" + chrom + "_kmers_36"
                fileKmer = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix=chromFragments,file_suffix=self.file_suffix)
                 
                self._files.append(fileKmer)
            
        return self._files
    
    def run(self):
        cfg = self.configuration
        gem.mapKmerWindows(input=self.getAllFiles()[1],output=self.files(),reference=self.getAllFiles()[0][1],threads=cfg["threads"])
    
    
    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on mapping Get Assembly Kmers!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0]
        
class CountKmerMappings(PipelineStep):
    """Once al fasta are mapped then is reported those kmer regions that were mapped at least more than 2 times"""
    def files(self):
        """Return the output files generated by this step."""
        if self._files is None:
            self._files = []
            for chrom in sorted(self.pipeline.chrLen.keys()):
                chromFragments = "_" + chrom + "_kmers_36"
                fileKmer = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix=chromFragments,file_suffix=self.file_suffix)
                 
                self._files.append(fileKmer)
            
        return self._files

    def run(self):
        gem.kmerCounts(input=self.getAllFiles()[0],output=self.files())   
    
    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on mapping Map Kmer Windows!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0]      
        
        
class DistributionPlot(PipelineStep):
    """Build accumulative distribution """
    def files(self):
        """Return the output files generated by this step."""
        if self._files is None:
            self._files = []
            
            rawCounts = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_raw_counts",file_suffix="bed")            
            plot = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_distribution",file_suffix="png")
            histoPlot = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_histogram",file_suffix="png")
            overThreshold = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_overrepresented",file_suffix="bed")
                 
            self._files.append(rawCounts)
            self._files.append(plot)
            self._files.append(histoPlot)
            self._files.append(overThreshold)
            
        return self._files       
    
    def run(self):
        cfg = self.configuration
        gem.distributionPlot(input=self.getAllFiles()[0],output=self.files(),hits=cfg["kmer_mappings_threshold"])
        
    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on mapping Map Kmer Windows!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0]  
    
class MaskFastaFromBed(PipelineStep):
    """Mask Fasta from Bed"""
    def files(self):
        """Return the output files generated by this step."""
        if self._files is None:
            self._files = []
            
            bedRegions = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_regions_to_mask",file_suffix="bed")
            kmerMaskedFasta = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_kmer_mask",file_suffix="fa")
            
            self._files.append(bedRegions)
            self._files.append(kmerMaskedFasta)

        return self._files

    def run(self):
        cfg = self.configuration
        gem.maskReference(basic_reference=self.getAllFiles()[0][1],kmer_over=self.getAllFiles()[1][3],bed_regions=self.files()[0],kmer_masked_fasta=self.files()[1],threads=cfg["threads"])

    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on distribution plot!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0]   

class CreatePaddedReference(PipelineStep):
    """Create padded reference masking 36 bp on flanking regions"""
    def files(self):
        """Return the output files generated by this step."""
        if self._files is None:
            self._files = []
            
            bedRegions = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_pad36",file_suffix="bed")
            bedChrLen = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_chr_len",file_suffix="bed")
            paddedFasta = self.pipeline.create_file_name(self.name, sub_directory=self.name, name_suffix="_pad36",file_suffix="fa")
            
            self._files.append(bedRegions)
            self._files.append(bedChrLen)
            self._files.append(paddedFasta)

        return self._files
        
    def run(self):
        """Run padded reference masking"""
        gem.pad36Fasta(basic_masked_regions=self.pipeline.bed_files,kmer_regions=self.getAllFiles()[0][0],masked_fasta=self.getAllFiles()[0][1],
                       chr_len=self.pipeline.chrLen,pad_regions=self.files()[0],
                       chrom_info_bed=self.files()[1],pad_fasta=self.files()[2])
        

    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on mask fasta from bed!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0] 
        
    
class DocumentationStep(PipelineStep):
    """Documentation Step"""
    
    def files(self):
        """Return the output files generated by this step. """
        if self._files is None:
            self._files = []
            
            html_file = self.pipeline.create_file_name(self.name, sub_directory=self.name,file_suffix="html")
            json_file = self.pipeline.create_file_name(self.name, sub_directory=self.name,file_suffix="json")    
            galculator_original = self.pipeline.create_file_name(self.name, sub_directory=self.name,file_suffix="original.txt")
            galculator_kmer = self.pipeline.create_file_name(self.name, sub_directory=self.name,file_suffix="kmer.txt")
            galculator_pad36 = self.pipeline.create_file_name(self.name, sub_directory=self.name,file_suffix="pad.txt")           
            
            self._files.append(html_file)
            self._files.append(json_file)
            self._files.append(galculator_original)
            self._files.append(galculator_kmer)
            self._files.append(galculator_pad36)
            
        return self._files
            
    def run(self):
        filesOut = self.files()
        html_doc = filesOut[0]        
        json_doc = filesOut[1]
        galculator_original = filesOut[2]        
        galculator_kmer = filesOut[3]
        galculator_pad36 = filesOut[4]
        
        inputsList = self.getAllFiles()
        plotPng = inputsList[0][1]
        histoPlot = inputsList[0][2]
        kmerRegions = inputsList[1][0]
        kmerFasta = inputsList[1][1]
        paddedRegions = inputsList[2][0]
        paddedFasta = inputsList[2][2]
        
        gem.galculator(input=self.pipeline.raw_fasta,output=galculator_original)
        gem.galculator(input=kmerFasta,output=galculator_kmer)
        gem.galculator(input=paddedFasta,output=galculator_pad36)        
        
        gem.assemblyReport.create_report(html_file=html_doc,json_file=json_doc,galculator_original_fasta=galculator_original,
                                         list_regions_mask=self.pipeline.bed_files,kmer_regions_mask=kmerRegions, accumulative_plot=plotPng,
                                         histogram_plot=histoPlot,galculator_kmer_mask_fasta=galculator_kmer,galculator_kmer_pad_fasta=galculator_pad36,
                                         regions_padded=paddedRegions)
                
    def getAllFiles(self):
        """Return the output of all
        dependencies"""
        if not self.dependencies:
            raise PipelineError("This step depends on Padded Reference!")
        return [self.pipeline.steps[i].files() for i in self.dependencies if i >= 0] 





        
class CopyNumberAssemblyPreparationPipeline(object):
    """General mapping pipeline class."""

    def __init__(self, args=None):
        self.steps = []  # pipeline steps
        self.run_steps = []  # steps to run
        
        self.output_dir = None  # Output directory
        self.name = None # target name

        self.raw_fasta = None #raw fasta assembly reference
        self.bed_files = [] #set of bed files for regions that must be masked
        self.chromsomeLengthFile = None #chromsome length file
        self.chrLen = {}
        self.mrcanavar = None #Path to the mrcanvar binary
        self.kmer_mappings_threshold = 20 #Number of mappings to considered a kmer regions as overrepresented

        self.threads = 1  # number of threads 
        self.write_config = None  # write configuration
        self.dry = False  # only dry run
        self.sort_memory = "768M"  # samtools sort memory
        self.force = False  # force computation of all steps
        
        self.remove_temp = True  # remove temporary
        
        if args is not None:
            # initialize from arguments
            # load configuration
            try:
                if args.load_configuration is not None:
                    self.load(args.load_configuration)
            except AttributeError:
                pass
            ## update parameter
            self.update(vars(args))
            ## initialize pipeline and check values
            self.initialize()
        
    def update(self, configuration):
        """Update configuration from given map

        configuration -- the input configuration
        """
        for k, v in configuration.items():
            try:
                if v is not None:
                    setattr(self, k, v)
            except AttributeError:
                pass

    def __update_dict(self, target, source):
        if source is None:
            return
        for k, v in source.items():
            #if v is not None:
            target[k] = v
            

    def basicMasking(self,name, description="Mask known regions, repeats and gaps", dependencies=None,configuration=None,final=False):
        """Basic Reference Masking"""
        step = BasicMaskingStep(name, final=final, dependencies=dependencies, description=description, file_suffix="fa")
        config = dotdict()
        
        config.raw_fasta = self.raw_fasta
        config.bed_files = self.bed_files
        
        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
        
        
    def getAssemblyKmers(self,name,dependencies=None,description="Create kmer windows from the assembly reference",configuration=None,final=False):  
        """Fragmentation of the assembly in kmer of 36 bp to be mapped against itself"""
        step = GetAssemblyKmers(name, final=final, dependencies=dependencies, description=description, file_suffix="fa")        
        config = dotdict()
        
        config.threads = self.threads 
        
        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
    
    def mapKmerWindows(self,name,dependencies=None,description="Kmer Mapping",configuration=None,final=True):
        """Map kmer fragmentation against the own assembly"""
        step = MapKmerWindows(name, final=final, dependencies=dependencies, description=description, file_suffix="map")
        config = dotdict()
        
        config.threads = self.threads
        
        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
    
    def countKmerMappings(self,name,dependencies=None,description="Count kmers mappings",configuration=None,final=False):
        """Count kmer mappings"""
        step = CountKmerMappings(name, final=final, dependencies=dependencies, description=description, file_suffix="bed")
        config = dotdict()
        
        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
        
    
    def distributionPlot(self,name,dependencies=None, description="Accumulative distribution plot",configuration=None,final=False):
        step = DistributionPlot(name, final=final, dependencies=dependencies, description=description)
        config = dotdict()
         
        config.kmer_mappings_threshold = self.kmer_mappings_threshold    
        
        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
    
    def maskFastaFromBed(self,name,dependencies=None,description="Mask Reference Fasta",configuration=None,final=False):
        """Mask fasta reference file """
        step = MaskFastaFromBed(name, final=final, dependencies=dependencies, description=description)
        config = dotdict()

        config.threads = self.threads

        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
    
    def createPaddedReference(self,name,dependencies=None,description="Create Padded reference fasta",configuration=None,final=False):
        """Mask with 36 bps flanking regions of repeats, kmers, simple repaeats and gaps"""
        step = CreatePaddedReference(name, final=final, dependencies=dependencies, description=description)
        config = dotdict()

        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
    
    def createReport(self,name, dependencies=None,description="Create Html Report",configuration=None,final=False):
        """Creates an html and json report """
        step = DocumentationStep(name, final=final, dependencies=dependencies, description=description)
        config = dotdict()

        if configuration is not None:
            self.__update_dict(config, configuration)
            
        step.prepare(len(self.steps), self, config)
        self.steps.append(step)
        return step.id
    
    def open_input(self,pair_end_files = None):
        return gem.files.open(self.input[0])
        
    def open_step(self, id, raw=False):
        """Open the original input files"""
        return self.steps[id].open(raw=raw)
        
    def load(self, file):
        """Load pipeline configuration from file"""
        if file is None or not os.path.exists(file):
            raise PipelineError("Configuration file not found: %s" % file)
        fd = open(file, "r")
        logging.gemtools.info("Loading configuraiton from %s", file)
        data = json.load(fd)
        for k, v in data.items():
            if hasattr(self, k):
                setattr(self, k, v)
        fd.close()
        
    def write_pipeline(self, file_name):
        """Write the pipeline and its configuration to a file
        based on the name
        """

        json_container = dict(self.__dict__)
        # skip the steps here, we convert them manually
        del json_container["steps"]
        del json_container["run_steps"]
        del json_container["write_config"]
        # remove non default values
        default_pipeline = CnvAssemblyPreparation()
        default_pipeline.initialize(silent=True)
        for k, v in json_container.items():
            if hasattr(default_pipeline, k) and getattr(default_pipeline, k) == v:
                del json_container[k]

        # json_container['piepline_steps'] = json_steps

        fd = open(file_name, "w")
        json.dump(json_container, fd, indent=2, sort_keys=True)
        fd.close()
        logging.gemtools.gt("Configuration saved to %s\n", file_name)
    
    def log_parameter(self):
        """Print selected parameters"""
        printer = logging.gemtools.gt
        run_step = len(self.run_steps) > 0

        printer("------------ Input Parameter ------------")
        printer("Raw Fasta        : %s", self.raw_fasta)
        printer("Bed Regions      : %s", self.bed_files)
        printer("")
        printer("MrCaNaVar path   : %s", self.mrcanavar)
        printer("kmer Maps    Lim.: %s", self.kmer_mappings_threshold)
        printer("")
        printer("Threads : %s", self.threads)  
        printer("")

        if not run_step:
            printer("------------ Pipeline Steps  ------------")
            for i, s in enumerate(self.steps):
                printer("%-2d - %20s : %s", i, s.name, s.description)
        else:
            printer("------------ Single Step execution  ------------")

        for i, s in enumerate(self.steps):
            if run_step and s.id not in self.run_steps:
                continue
            printer("")
            if len(s.dependencies) > 0:
                printer("------------ [ID:{0:-3} -- '{1}'] [Depends On: {2}] ------------".format(i, s.name, ", ".join([self.steps[j].name for j in s.dependencies])))
            else:
                printer("------------ [ID:{0:-3} -- '{1}'] ------------".format(i, s.name))

            for k, v in s.configuration.items():
                printer("%25s : %s", k, str(v))

            for i, f in enumerate(s.files()):
                if i == 0:
                    printer("%25s : %s", "Temporary Outputs", not s.final)
                    printer("%25s : %s", "Outputs", f)
                else:
                    printer("%25s : %s", "", f)
        printer("--------------------------------------------------------------")
        printer("")        

    def initialize(self, silent=False):
        # check general parameter
        errors = []
       
        if self.raw_fasta is None:
            errors.append("No input fasta file specified")

            
        if self.name is None and self.raw_fasta is not None:
            # get name from input files
            name = os.path.basename(self.raw_fasta)
            if name.endswith(".gz"):
                name = name[:-3]
            idx = name.rfind(".")
            if idx > 0:
                self.name = name[:idx]

        if not os.path.exists(self.raw_fasta):
            errors.append("Fasta reference file not found: %s" % (self.raw_fasta))
        else:
            self.raw_fasta = os.path.abspath(self.raw_fasta)
            
        for bedFile in self.bed_files:
            if not os.path.exists(bedFile):
                errors.append("Bed file not found: %s" % (bedFile))

        if self.output_dir is None:
            self.output_dir = os.getcwd()

        self.output_dir = os.path.abspath(self.output_dir)

        if self.threads <= 0:
            self.threads = 1

        if self.chromsomeLengthFile is None:
            errors.append("No chromosome length file specified")
        else: 
            if not os.path.exists(self.chromsomeLengthFile):
                errors.append("Chromosome length file not found: %s" % (self.chromsomeLengthFile))
            else:
                self.chromsomeLengthFile = os.path.abspath(self.chromsomeLengthFile)
                
                with open(self.chromsomeLengthFile, "r") as chrLen:
                    for line in chrLen:
                        fields = re.split('\s+', line.rstrip('\n'))
                        if len(fields) > 1:
                            if str(fields[1]).isdigit():
                                self.chrLen [fields[0]] = fields[1]
                            else:
                                errors.append("Field two %s in %s is not a numeric one" % (fields[1],self.chromsomeLengthFile))


        if not silent and len(errors) > 0 and self.write_config is None:
            raise PipelineError("Failed to initialize neccessary parameters:\n\n%s" % ("\n".join(errors)))
        if self.write_config is not None:
            # log configuration errors
            logging.gemtools.warning("---------------------------------------------")
            logging.gemtools.warning("Writing configuration")
            logging.gemtools.warning("")
            logging.gemtools.warning("Note that some of the parameters are missing:\n")
            for e in errors:
                logging.gemtools.warning("\t" + str(e))
            logging.gemtools.warning("---------------------------------------------") 
            
            
    def run(self):
        run_step = len(self.run_steps) > 0

        if self.write_config is not None:
            self.write_pipeline(self.write_config)
            return
        if self.dry:
            # check and print states
            print "Checking Job states"
            print "-----------------------------------------------------"
            for step in self.steps:
                print "Step %3s:%25s :: %s" % (str(step.id), step.name,
                                         "Done" if step.is_done() else "Compute")
            return

        error = False

        all_done = True
        final_files = []
        # check final steps if we are not running a set of steps
        if not run_step and not self.force:
            for step in self.steps:
                if step.final:
                    final_files.extend(step.files())
                    all_done = all_done & step.is_done()
            if all_done:
                logging.gemtools.warning("The following files already exist. Nothing to be run!\n\n%s\n" % ("\n".join(final_files)))
                return

        time = Timer()
        times = {}

        if run_step:
            # sort by id
            self.run_steps.sort()
        ids = [s.id for s in self.steps]
        if run_step:
            ids = self.run_steps

        step = None

        # register signal handler to catch
        # interruptions and perform cleanup
         # register cleanup signal handler
        def cleanup_in_signal(signal, frame):
            logging.gemtools.warning("Job step canceled, forcing cleanup!")
            step.cleanup(force=True)

        signal.signal(signal.SIGINT, cleanup_in_signal)
        signal.signal(signal.SIGQUIT, cleanup_in_signal)
        signal.signal(signal.SIGHUP, cleanup_in_signal)
        signal.signal(signal.SIGTERM, cleanup_in_signal)

        for step_id in ids:
            step = self.steps[step_id]

            if run_step:
                # check dependencies are done
                for d in step.dependencies:
                    if not self.steps[d].is_done():
                        logging.gemtools.error("Step dependency is not completed : %s", self.steps[d].name)
                        error = True
                        break
            if run_step or self.force or not step.is_done():
                logging.gemtools.gt("Running step: %s" % step.name)
                t = Timer()

                if not os.path.exists(self.output_dir):
                    # make sure we create the ouput folder
                    logging.gemtools.warn("Creating output folder %s", self.output_dir)
                    try:
                        os.makedirs(self.output_dir)
                    except OSError as exc: # Python >2.5
                        if exc.errno == errno.EEXIST and os.path.isdir(path):
                            pass
                        else:
                            logging.gemtools.error("unable to create output folder %s", self.output_dir)
                            error = True
                            break

                try:
                    step.run()
                except KeyboardInterrupt:
                    logging.gemtools.warning("Job step canceled, forcing cleanup!")
                    error = True
                    step.cleanup(force=True)
                    break
                except PipelineError, e:
                    logging.gemtools.error("Error while executing step %s : %s" % (step.name, str(e)))
                    logging.gemtools.warning("Cleaning up after failed step : %s", step.name)
                    step.cleanup(force=True)
                    error = True
                    break
                except gem.utils.ProcessError, e:
                    logging.gemtools.error("Error while executing step %s : %s" % (step.name, str(e)))
                    logging.gemtools.warning("Cleaning up after failed step : %s", step.name)
                    step.cleanup(force=True)
                    error = True
                    break
                except Exception, e:
                    traceback.print_exc()
                    logging.gemtools.error("Error while executing step %s : %s" % (step.name, str(e)))
                    logging.gemtools.warning("Cleaning up after failed step : %s", step.name)
                    step.cleanup(force=True)
                    error = True
                    break
                finally:
                    t.stop(step.name + " completed in %s", loglevel=None)
                    times[step.id] = t.end
                    if not error:
                        logging.gemtools.gt("Step %s finished in : %s", step.name, t.end)
                    else:
                        logging.gemtools.gt("Step %s failed after : %s", step.name, t.end)
            else:
                logging.gemtools.warning("Skipping step %s, output already exists" % (step.name))

        # do celanup if not in error state
        if not error:
            logging.gemtools.debug("Cleanup after run")
            for step in self.steps:
                step.cleanup()

            time.stop("Completed in %s", loglevel=None)
            logging.gemtools.gt("Step Times")
            logging.gemtools.gt("-------------------------------------")
            for s in self.steps:
                if s.id in times:
                    logging.gemtools.gt("{0:>25} : {1}".format(s.name, times[s.id]))
                else:
                    logging.gemtools.gt("{0:>25} : skipped".format(s.name))
            logging.gemtools.gt("-------------------------------------")
            logging.gemtools.gt("Pipeline run finshed in %s", time.end)
            
    def create_file_name(self,suffix,sub_directory=None, name_suffix=None, file_suffix="map", final=False):
        """Create a result file name"""
        file = ""
        
        if final:
            suffix = None
        if name_suffix is None:
            name_suffix = ""
            
        out_path = self.output_dir
        if sub_directory is not None:
            out_path = out_path + "/" + sub_directory
            
        #Checking output path    
        if not os.path.exists(out_path):
            os.makedirs(out_path)
             
        if suffix is not None and len(suffix) > 0:
            file = "%s/%s%s_%s.%s" % (out_path, self.name, name_suffix, suffix, file_suffix)
        else:
            file = "%s/%s%s.%s" % (out_path, self.name, name_suffix, file_suffix)
           
        return file
    
    def register_parameter(self, parser):
        """Register all parameters with the given
        arparse parser"""
        self.register_general(parser)
        self.register_overRepresentedKmer(parser)
        self.register_canavarPreparation(parser)
        self.register_execution(parser)
        
    def register_general(self,parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Input')
        ## general pipeline paramters
        input_group.add_argument('-f', '--fasta', dest="raw_fasta", metavar="unmasked_reference.fa",
            help='''Unmasked genome reference. It is not mandatory to have repeats, simple repeats or gaps masked.''')
        input_group.add_argument('--bed-regions', dest="bed_files", nargs="+", 
                                 help='''Set of bed files for regions to be masked. Can be one or more bed files. It is expected to have at
                                 least gaps, repeats from Repeat Masker and simple repeats from Tandem Repeat finder.''')

        input_group.add_argument('--chr-len', dest="chromsomeLengthFile", metavar="chrLenFile",
            help='''File containing information about contig name and its length. First column must be the name and the second the length.''')     
        
        input_group.add_argument('-t', '--threads', dest="threads", metavar="NUM_THREADS", 
                                 help='Number of threads. Default to %d' % self.threads)
     
    def register_overRepresentedKmer(self,parser):
        """Register over represented kmer parameters with the given
        argparse parser

        parser -- the argparse parser
        """                
        over_represented_group = parser.add_argument_group('Kmer Windows Over Representation')
        over_represented_group.add_argument('-kmer-mappings-threshold', dest="kmer_mappings_threshold", metavar="NUM_MAPPINGS", 
                                 help='Number of mappings to consider a kmer regions as overrepresented. Default to %d' % self.kmer_mappings_threshold)
        
    def register_canavarPreparation(self,parser):
        """Register canvar preparation with the given
        argparse parser

        parser -- the argparse parser
        """                
        over_represented_group = parser.add_argument_group('mrCaNaVar Assembly Preparation Step')
        over_represented_group.add_argument('-mrcanavar-path', dest="mrcanavar", metavar="MRCANAVAR_PATH", 
                                        help='''Path to the mrcanavar binary''')    
    def register_execution(self,parser):
        """Register the execution mapping parameters with the
        given arparse parser

        parser -- the argparse parser
        """
        execution_group = parser.add_argument_group('Execution')
        execution_group.add_argument('--save', dest="write_config", nargs="?", const=None, help="Write the given configuration to disk")
        execution_group.add_argument('--dry', dest="dry", action="store_true", default=None, help="Print and write configuration but do not start the pipeline")
        execution_group.add_argument('--load', dest="load_configuration", default=None, metavar="cfg", help="Load pipeline configuration from file")
        execution_group.add_argument('--run', dest="run_steps", type=int, default=None, nargs="+", metavar="cfg", help="Run given pipeline steps idenfified by the step id")
        execution_group.add_argument('--force', dest="force", default=None, action="store_true", help="Force running all steps and skip checking for completed steps")
       
