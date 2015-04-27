#!/usr/bin/env python

import json
import os
import re
import numpy

from gem.cnvReport import BaseDocumentHtml,BaseStatistical

class BuildHtmlAssembly(BaseDocumentHtml):
    """Base Mapping step"""
    
    def titleDocument(self):
        '''Title Web Document'''
        self.vContent.append("  <H1 id=\"title\"> <U> ASSEMBLY MASKING FOR CNV PIPELINE REPORT </U> </H1>\n")
    
    
    def linkersTable(self):
        ''' Links table section '''
        self.vContent.append("<a id=\"linkers\"></a>\n")

        self.vContent.append("<table id=\"linksTable-b\" >\n")
        self.vContent.append("    <tbody>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#originalReference\">Original Reference</a> </td>\n")
        self.vContent.append("        </tr>\n") 
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#basicMaskingRegions\">Basic Masked Regions</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#kmerMaskingRegions\">Kmer Masking Regions</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#accumulativePlot\">Acummulative Distribution Plot</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"histogramPlot\">Histogram Plot</a> </td>\n")
        self.vContent.append("        </tr>\n")        
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#kmerMaskedReference\">Kmer Masked Reference</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#kmerMaskedReferencePad36\">Kmer Masked Reference Pad 36</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#regionMaskedPad36\">Regions Masked 36 bps around any gap/repeat</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("    </tbody>\n")
        self.vContent.append("</table>\n")

        self.vContent.append(" <br> \n")    
    
    def sectionTitle(self, typeSection):
        ''' Creates a linker anchor for the section to be open, it also set section name'''
	
        linker = ""
        sectionName = ""
        
        if typeSection == "originalReference":
            linker += "id=\"originalReference\""
            sectionName += "Original Reference"
        elif typeSection == "basicMaskingRegions":
            linker += "id=\"basicMaskingRegions\""
            sectionName += "Basic Masked Regions"
        elif typeSection == "kmerMaskingRegions":
            linker += "id=\"kmerMaskingRegions\""
            sectionName += "Kmer Masking Regions"
        elif typeSection == "accumulativePlot":
            linker += "id=\"accumulativePlot\""
            sectionName += "Acummulative Distribution Plot"
        elif typeSection == "histogramPlot":
            linker += "id=\"histogramPlot\""
            sectionName += "Histogram Plot"
        elif typeSection == "kmerMaskedReference":
            linker += "id=\"kmerMaskedReference\""
            sectionName += "Kmer Masked Reference"
        elif typeSection == "kmerMaskedReferencePad36":
            linker += "id=\"kmerMaskedReferencePad36\""
            sectionName += "Kmer Masked Reference Pad 36"
        elif typeSection == "regionMaskedPad36":
            linker += "id=\"regionMaskedPad36\""
            sectionName += "Regions Masked 36 bps around any gap/repeat"
        
        self.vContent.append("<a "+ linker +"></a>\n")
        self.vContent.append("<H1 id=\"section\"> " + sectionName + " </H1>\n")
    

class FastaComposition(BaseStatistical):
    '''Class for manage mrCanvar stats'''

    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Name")        
        self.fieldsList.append("Contigs")
        self.fieldsList.append("Total (base pairs)")
        self.fieldsList.append("A")
        self.fieldsList.append("C")
        self.fieldsList.append("G")
        self.fieldsList.append("T")
        self.fieldsList.append("N")
        
    def parseInput(self):
        '''Parse input source'''
        self.fieldValue['Name'] = os.path.basename(self.source_file)
        with open(self.source_file, "r") as statsFile:
            for line in statsFile:
                if line.find("|") != -1:
                    vFields = re.split('\s+',line.rstrip('\n'))
                    self.fieldValue['contigs'] = vFields[0]
                    self.fieldValue['total'] = long(vFields[1])
                    self.fieldValue['A'] = long(vFields[4])
                    self.fieldValue['C'] = long(vFields[6])
                    self.fieldValue['G'] = long(vFields[8])
                    self.fieldValue['T'] = long(vFields[10])
                    self.fieldValue['N'] = long(vFields[12])

    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append(self.fieldValue['Name'])
        self.valuesList.append(self.fieldValue['contigs'])
        self.valuesList.append("{:,}".format(self.fieldValue['total'])) 
        self.valuesList.append("{:,}".format(self.fieldValue['A'])) 
        self.valuesList.append("{:,}".format(self.fieldValue['C'])) 
        self.valuesList.append("{:,}".format(self.fieldValue['G']))
        self.valuesList.append("{:,}".format(self.fieldValue['T']))
        self.valuesList.append("{:,}".format(self.fieldValue['N']))        
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["FastaComposition_"+os.path.basename(self.source_file)] = self.fieldValue 
        return self.allValues
        
class BedStats(BaseStatistical):
    ''' BED file Statistics '''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("")        
        self.fieldsList.append("Regions")
        self.fieldsList.append("Base Pairs")
        
    def parseInput(self):
        '''Parse input source'''        
        totalLines = 0
        totalBases = 0
        with open(self.source_file, "r") as bedData:
            for line in bedData:
                totalLines = totalLines + 1
                vFields = line.rstrip('\n').split("\t")
                totalBases = totalBases + (long(vFields[2]) - long(vFields[1]))

        self.fieldValue['name'] = os.path.basename(self.source_file)
        self.fieldValue['regions'] = totalLines
        self.fieldValue['bases'] = totalBases
 
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append(self.fieldValue["name"])
        self.valuesList.append("{:,}".format(self.fieldValue["regions"]))
        self.valuesList.append("{:,}".format(self.fieldValue["bases"]))
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues [os.path.basename(self.source_file)] = self.fieldValue 
        return self.allValues
 
class AccumulativeDistributionPlot(BaseStatistical):
    ''' Accumulative Distribution Plot '''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Accumulative distribution plot")
    
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append(self.source_file)  
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["CumulativeDistribution"] = self.source_file 
        return self.allValues 
        
def create_report(html_file,json_file,galculator_original_fasta,list_regions_mask,kmer_regions_mask, accumulative_plot,
                  histogram_plot,galculator_kmer_mask_fasta, galculator_kmer_pad_fasta,regions_padded):
    """ Generate HTML Report from mapping stat, mrCanavar log, control regions distribution, control regions plot and cutoffs file 
    
    Parameters
    ----------
    html_file: html file to store the document 
    json_file: json file to store the document
    galculator_original_fasta: output file from galculator for the original reference
    list_regions_mask: list of bed files for basic masking
    kmer_regions_mask: kmer masked regions bed file
    accumulative_plot: accumulative plot png image
    galculator_kmer_mask_fasta: output file from galculator for the kmer masked reference
    galculator_kmer_pad_fasta: output file from galculator for the padded masked reference
    regions_padded: padded regions masked bed
    """
    
    #Original Fasta Reference Nucleotide Composition
    originalFasta = FastaComposition(source_file=galculator_original_fasta,is_json = False)

    #Basic Bed Regions
    basicRegions = []
    for bedFile in list_regions_mask:
        basicRegions.append(BedStats(source_file=bedFile,is_json = False))

    #Kmer Masking (regions,bases)
    kmerMasking = BedStats(source_file=kmer_regions_mask,is_json = False)  
    
    #Accumulative Plot
    pngPlot = AccumulativeDistributionPlot(source_file=accumulative_plot,is_json = False) 
    
    #Histogram Plot
    histogramPlot = AccumulativeDistributionPlot(source_file=histogram_plot,is_json = False)
    
    #Kmer Masking Compositon
    kmerMaskingComposition = FastaComposition(source_file=galculator_kmer_mask_fasta,is_json = False)
    
    #Kmer Masking Padded Composition
    kmerMaskingPadded = FastaComposition(source_file=galculator_kmer_pad_fasta,is_json = False)
    
    #Regions Kmer Masked Pad 36
    pad36Masked = BedStats(source_file=kmer_regions_mask,is_json = False)  
    
    #1. HTML REPORT CONSTRUCTION
    #1.1 Header HTML
    vHtmlContent = []
    htmlManager = BuildHtmlAssembly(vHtmlContent)
    htmlManager.addHtmlReportHeader()
    
    #1.2 Section Original Fasta
    htmlManager.addHtmlNewSection(originalFasta,"originalReference","Original Reference Nucleotide Composition","blue")
    
    #1.3 Section Basic masked Regions
    htmlManager.addHtmlNewSection(basicRegions,"basicMaskingRegions","Basic Masked Regions","green",is_stack=True)
    
    #1.4 Section Kmer Masking Regions
    htmlManager.addHtmlNewSection(kmerMasking,"kmerMaskingRegions","Kmer Masking Regions","blue")
    
    #1.5 Section Accumulative plot
    htmlManager.addHtmlNewSection(pngPlot,"accumulativePlot","Accumulative Distribution","green",is_image=True)
    
    #1.6 Section Histogram plot
    htmlManager.addHtmlNewSection(histogramPlot,"histogramPlot","Histogram Plot","blue",is_image=True)
        
    #1.7 Section Kmer Masked composition
    htmlManager.addHtmlNewSection(kmerMaskingComposition,"kmerMaskedReference","Kmer Masked Reference Nucleotide Composition","green")
    
    #1.8 Section Kmer Masked Padded Composition
    htmlManager.addHtmlNewSection(kmerMaskingPadded,"kmerMaskedReferencePad36","Pad 36 Masked Reference Nucleotide Composition","blue")
    
    #1.9 Section Pad 36 Masked Regions
    htmlManager.addHtmlNewSection(pad36Masked,"regionMaskedPad36","Pad 36 Masked Regions","green")
    
    #1.10 Close Html Report
    htmlManager.closeHtmlReport()
    
    #2. SAVE HTML DOCUMENT
    htmlManager.saveDocument(html_file)
    
    #3. CREATE CASCADE STYLE SHEET
    vCSScontent = []
    cssManager = BuildHtmlAssembly(vCSScontent)
    cssManager.buildStyleSheet()

    cssManager.saveDocument(os.path.dirname(os.path.abspath(html_file)) + "/style.css")
    
    #4.CREATE JSON
    jsonDataDocument = {}
    jsonDataDocument.update(originalFasta.getHeaderValues())
    
    for basicRegion in basicRegions:
        jsonDataDocument.update(basicRegion.getHeaderValues())
        
    jsonDataDocument.update(kmerMasking.getHeaderValues())
    jsonDataDocument.update(pngPlot.getHeaderValues())
    jsonDataDocument.update(histogramPlot.getHeaderValues())
    jsonDataDocument.update(kmerMaskingComposition.getHeaderValues())
    jsonDataDocument.update(kmerMaskingPadded.getHeaderValues())
    jsonDataDocument.update(pad36Masked.getHeaderValues())
   
    with open(json_file, 'w') as of:
        json.dump(jsonDataDocument, of, indent=2)
    
    