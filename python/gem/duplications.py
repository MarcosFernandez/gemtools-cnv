#!/usr/bin/env python
"""Management of the two methods for calling Duplications"""
import os
import gem
import re

from . import utils

class Duplications(object):

    """General parameters class"""
    def __init__(self, sample, bedRepeatRegions=None, gapsBedCoord=None,cutOffFile=None,pathMrCanavar=None,outputlist=None):
        self.sample = sample        
        self.bedRepeatRegions = bedRepeatRegions
        self.gapsBedCoord = gapsBedCoord
        self.cutOffFile = cutOffFile
        self.mrcanavar_path = self.getCanavarPath(pathMrCanavar)
        self.copy_number_file = self.getCopyNumber()
        self.mrcanavar_log = self.getMrCanavarLog()
        self.swCall = self.getShortWindowsPath()
        self.lwCall = self.getLongWindowsPath()
        self.cwCall = self.getCopyNumberWindowsPath()
   
        #Process output file
        self.duplication_file = None
        self.duplication_without_gaps_file = None
        self.filteredSwCall = None
        self.filteredLwCall = None
        self.filteredCwCall = None
        self.wssdPicked = None
        self.wssdMerged = None
        self.wssdMerged10K = None
        self.wssdMerged10KNoGaps = None 
        
        self.prepare_outputs(outputlist)
        
        
        ## paths to the executables
        self.executables = {
             "wssd_picker": "wssd_picker.pl",
             "twoPartOvpMgsrt": "twoPartOvpMgsrt", 
             "MergeCoord" : "MergeCoord"
        }

    def getCanavarPath(self,pathMrCanavar):
        """Return directory name of the file poited by pathMrCanavar"""
        return os.path.dirname(pathMrCanavar)
        
    def getCopyNumber(self):
        """Return the path to copy number calls"""
        return gem.utils.runFindCommand(self.mrcanavar_path,"*.calls.copynumber.bed").rstrip('\n')   
        
    def getMrCanavarLog(self):
        """Return the path to mr canavar logs"""
        return gem.utils.runFindCommand(self.mrcanavar_path,"*.calls.log").rstrip('\n')
        
    def getLongWindowsPath(self):
        """Return long windows path"""
        return gem.utils.runFindCommand(self.mrcanavar_path,"*.calls.lw_norm.bed").rstrip('\n')   
        
    def getShortWindowsPath(self):
        """Return short windows path"""
        return gem.utils.runFindCommand(self.mrcanavar_path,"*.calls.sw_norm.bed").rstrip('\n')

    def getCopyNumberWindowsPath(self):
        """Return copy number windows path"""
        return gem.utils.runFindCommand(self.mrcanavar_path,"*.calls.cw_norm.bed").rstrip('\n')   
        
    def prepare_outputs(self, ouputlist):
        """Outputs prepare"""        
        self.duplication_file = ouputlist[0]
        self.duplication_without_gaps_file = ouputlist[1]
        self.filteredSwCall = ouputlist[2]
        self.filteredLwCall = ouputlist[3]
        self.filteredCwCall = ouputlist[4]
        self.wssdPicked = ouputlist[5]
        self.wssdMerged = ouputlist[6]
        self.wssdMerged10K = ouputlist[7]
        self.wssdMerged10KNoGaps = ouputlist[8]      
       
    def getSamplesCuttoff(self):
        """Reads a Cutoff file and stores in a dictionary object (hash) each sample name and its cutoff value 
        Returns dictionary instance object """
        sampleCutoffs = {}
        with open(self.cutOffFile, "r") as cutFile:
            for lineSample in cutFile:
                fields = re.split('\s+', lineSample.rstrip('\n'))
                sampleCutoffs[fields[0]] =  fields[5]

        return sampleCutoffs
        
    def runMethod1Duplications(self):
        """ Method 1 for Calling Duplications.
        1. Selects 1-Kbps windows with copy number exceeding the sample-specific cutoff (but lower than 100 copies), merging adjacent windows into single regions.
        2. Selects regions with at least five windows and larger than 10 Kbp
        3. Retains duplications whose at least 85% of their size does not overlap with repeats
        
        duplicationFile -- duplications calls output file     
        duplicationWoGapsFile -- duplications calls output file without coordinates
        """
        gain = self.getSamplesCuttoff()[self.sample]
        first=["sed","1,2d",self.copy_number_file]
        tools = [first]
        threshold = ['awk','{if('+ gain +'<=$5 && $5<100) print $1"\t"$2"\t"$3"\t"$5}']
        tools.append(threshold)
        mergeFile = ["mergeBed","-i","stdin","-n"]
        tools.append(mergeFile)
        coverage = ['awk','{if($4>=5 && $3-$2>10000) print $0}']
        tools.append(coverage)
        repeats = ["intersectBed","-wo","-a","stdin","-b",self.bedRepeatRegions]
        tools.append(repeats)
        groups = ["groupBy","-i","stdin","-g","1,2,3,4","-c","8","-o","sum"]
        tools.append(groups)
        
        filter = ['awk','-F','\t','{if($5 < 0.85*($3-$2)){ print $1"\t"$2"\t"$3}}']
        tools.append(filter)
       
        
        process = utils.run_tools(tools, input=None, output=self.duplication_file, name="m1-calls-1")
        if process.wait() != 0:
            raise ValueError("M1 Call Duplications analysis failed!")
        
        second = ["subtractBed","-a",self.duplication_file,
                  "-b",self.gapsBedCoord]
                  
        tools = [second]
        
        process = utils.run_tools(tools, input=None, output=self.duplication_without_gaps_file, name="m1-calls-2")
        if process.wait() != 0:
            raise ValueError("M1 Call Duplications analysis failed!")
                             
    def getAvgStDev(self,typeWindow):  
        """ CW Get Average Read Depth and StDev for 1kbp non repetitive sequence non overlapping 
            LW Get Average Read Depth and StDev for 5kbp non repetitive sequence overlapping 1kbp of real sequence
            SW Get Average Read Depth and StDev for 1kbp non repetitive sequence overlapping 1kbp of real sequence
        """
        meanAvg = []
        with open(self.mrcanavar_log, "r") as logFile:
            for lineLog in logFile:
                fields = lineLog.rstrip('\n').split()
                if len(fields) > 0:
                    if fields[0] == typeWindow:
                        meanAvg.append(float(fields[4][:-1]))
                        meanAvg.append(float(fields[7][:-1]))
                        return meanAvg

    def callsChromFiltering(self,fileInput, fileOutput):
        """ Removes sexual and mitocondrial chromosomes """
        #Open output file descriptor
        with open(fileOutput, 'w') as goalFile:
            #Process input file
            with open(fileInput, "r") as originalFile:
                for lineOriginal in originalFile:
                    fields = lineOriginal.rstrip('\n').split('\t')
                    if len(fields) > 0:
                        if fields[0] != "chrX" and fields[0] != "chrY" and fields[0] != "chrM" \
                           and fields[0] != "chrY" and fields[0] != "#CHROM":
                            goalFile.write(lineOriginal)   
                    

    def runMethod2Duplications(self):
        cwAvgDev = self.getAvgStDev("CW")   
        lwAvgDev = self.getAvgStDev("LW")
        swAvgDev = self.getAvgStDev("SW")
        #Chromosome filtering mrCanavar Calls
        self.callsChromFiltering(self.swCall, self.filteredSwCall)
        self.callsChromFiltering(self.lwCall, self.filteredLwCall)
        self.callsChromFiltering(self.cwCall, self.filteredCwCall)
        #WSSD
        thr5 = lwAvgDev[0] + (4 * lwAvgDev[1])
        thr1 = swAvgDev[0] + (4 * swAvgDev[1]) 
        #Pick those 7 windows of which 6 should have more than the average + 4std of the distribution or read depth in the control regions
        commandPicker = [self.executables['wssd_picker'],
                         '-f',self.filteredLwCall,
                         '-w',"7",
                         '-s',"6",
                         '-b',"4",
                         '-k',self.filteredSwCall,
                         '-n',"5",
                         '-i',"1",
                         '-c', str(thr5),
                         '-t', str(thr1),
                         '-o', self.wssdPicked
        ]                         
                         
                         
        tools = [commandPicker]
        
        process = utils.run_tools(tools, input=None, output=None, name="m2-calls-1")
        if process.wait() != 0:
            raise ValueError("M2 Call Duplications analysis failed!")
            
        #MERGE COORDINATES
        commandMerge = [self.executables['MergeCoord'],
                         '-f1',self.wssdPicked,
                         '-h1',
                         '-outF',self.wssdMerged
        ]        

        tools = [commandMerge]        
        
        process = utils.run_tools(tools, input=None, output=None, name="m2-calls-2")
        if process.wait() != 0:
            raise ValueError("M2 Call Duplications analysis failed!")
        
        
        #GET 10K size regions
        tenKiloBpRegions = ['awk', '-F', '\t', '{if($3 - $2 >= 10000){print $0}}',self.wssdMerged]
        tools = [tenKiloBpRegions]

        process = utils.run_tools(tools, input=None, output=self.wssdMerged10K, name="m2-calls-3")
        if process.wait() != 0:
            raise ValueError("M2 Call Duplications analysis failed!")
       
        #REMOVE GAPS
        commandGaps = [self.executables['twoPartOvpMgsrt'],
                         '-i',self.wssdMerged10K,
                         '-f',
                         '-j',self.gapsBedCoord,
                         '-t',
                         '-L',
                         '-o',self.wssdMerged10KNoGaps
        ]            
        
        tools = [commandGaps]
        
        process = utils.run_tools(tools, input=None, output=None, name="m2-calls-4")
        if process.wait() != 0:
            raise ValueError("M2 Call Duplications analysis failed!")
        