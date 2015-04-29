#!/usr/bin/env python

import json
import os
import re
import numpy

class BaseDocumentHtml(object):

    def __init__(self,vectorContents):
        self.vContent = vectorContents
        

    def addHtmlReportHeader(self):
        ''' Add HTML Report Header to a given vector'''
        self.vContent.append("<HTML>\n")

        self.vContent.append(" <HEAD>\n")
        self.vContent.append(" <STYLE TYPE=\"text/css\">\n")

        self.vContent.append("  <!--\n")
        self.vContent.append("   @import url(\"style.css\"); \n")
        self.vContent.append("  -->\n")

        self.vContent.append(" </STYLE>\n")
        self.vContent.append(" </HEAD>\n")

        self.vContent.append(" <BODY>\n")
        self.titleDocument()

        self.linkersTable()

    def closeHtmlReport(self):
        ''' CLOSE HTML header to a given vector object '''  
        self.vContent.append(" </BODY>\n")
        self.vContent.append("</HTML>\n")
	
    def titleDocument(self):
        '''To Be Defined in the Child Class'''
        pass
 
    def linkersTable(self):
        ''' To Be Difined in the Child Class '''
        pass
        
    def buildStyleSheet(self):
        ''' Add HTML header to a given vector '''
        self.vContent.append(" #title\n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 20px;\n")
        self.vContent.append("   text-align: center;\n")
        self.vContent.append("   color: #039;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #section\n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 18px;\n")
        self.vContent.append("   text-align: center;\n")
        self.vContent.append("   color: #039;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        #LINKS TABLE
        self.vContent.append("#linksTable-b\n")
        self.vContent.append("{\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("	font-size: 12px;\n")
        self.vContent.append("	background: #fff;\n")
        self.vContent.append("	margin: 5px;\n")
        self.vContent.append("	width: 200px;\n")
        self.vContent.append("	border-collapse: collapse;\n")
        self.vContent.append("	text-align: left;\n")
        self.vContent.append("}\n")
        self.vContent.append("#linksTable-b th\n")
        self.vContent.append("{\n")
        self.vContent.append("	font-size: 14px;\n")
        self.vContent.append("	font-weight: normal;\n")
        self.vContent.append("	color: #039;\n")
        self.vContent.append("	padding: 10px 8px;\n")
        self.vContent.append("	border-bottom: 2px solid #6678b1;\n")
        self.vContent.append("}\n")
        self.vContent.append("#linksTable-b td\n")
        self.vContent.append("{\n")
        self.vContent.append("	border-bottom: 1px solid #ccc;\n")
        self.vContent.append("	color: #669;\n")
        self.vContent.append("	padding: 6px 8px;\n")
        self.vContent.append("}\n")
        self.vContent.append("linksTable-b tbody tr:hover td\n")
        self.vContent.append("{\n")
        self.vContent.append("	color: #009;\n")
        self.vContent.append("}\n")
        #BLUE TABLE
        self.vContent.append(" #hor-zebra\n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 12px;\n")
        self.vContent.append("   margin: 0px;\n")
        self.vContent.append("   width: 100%;\n")
        self.vContent.append("   text-align: left;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #hor-zebra th\n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-size: 14px;\n")
        self.vContent.append("   font-weight: normal;\n")
        self.vContent.append("   padding: 10px 8px;\n")
        self.vContent.append("   color: #039;\n")
        self.vContent.append("   border-bottom: 2px solid #6678b1;\n")
        self.vContent.append("   border-right: 1px solid #6678b1; \n")
        self.vContent.append("	border-left: 1px solid #6678b1;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #hor-zebra td\n")
        self.vContent.append(" {\n")
        self.vContent.append("   padding: 8px;\n")
        self.vContent.append("   color: #669;\n")
        self.vContent.append("   border-right: 1px solid #6678b1; \n")
        self.vContent.append("	border-left: 1px solid #6678b1;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #hor-zebra .odd\n")
        self.vContent.append(" {\n")
        self.vContent.append("   background: #e8edff;\n")
        self.vContent.append("   border-right: 1px solid #6678b1; \n")
        self.vContent.append("	border-left: 1px solid #6678b1;\n")
        self.vContent.append(" }\n")
        #GREEN TABLE
        self.vContent.append(" #green \n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 12px;\n")
        self.vContent.append("   margin: 0px;\n")
        self.vContent.append("   width: 100%;\n")
        self.vContent.append("   text-align: left;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #green th\n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-size: 14px;\n")
        self.vContent.append("   font-weight: normal;\n")
        self.vContent.append("   padding: 10px 8px;\n")
        self.vContent.append("   color: #2b9900;\n")
        self.vContent.append("   border-bottom: 2px solid #66b16f;\n")
        self.vContent.append("   border-right: 1px solid #66b16f; \n")
        self.vContent.append("	border-left: 1px solid #66b16f;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #green td\n")
        self.vContent.append(" {\n")
        self.vContent.append("   padding: 8px;\n")
        self.vContent.append("   color: #578055;\n")
        self.vContent.append("   border-right: 1px solid #66b16f; \n")
        self.vContent.append("	border-left: 1px solid #66b16f;\n")
        self.vContent.append(" }\n")
        self.vContent.append(" #green .odd\n")
        self.vContent.append(" {\n")
        self.vContent.append("   background: #eaffe8;\n")
        self.vContent.append("   border-right: 1px solid #66b16f; \n")
        self.vContent.append("	border-left: 1px solid #66b16f;\n")
        self.vContent.append(" }\n")
        #LINKS
        self.vContent.append("a.link:link {font-family: \"Sans-Serif\";font-size: 12px;color: #039;text-decoration:none;}")
        self.vContent.append("a.link:visited {font-family: \"Sans-Serif\";font-size: 12px;color: #039;text-decoration:none;}")
        self.vContent.append("a.link:hover {font-family: \"Sans-Serif\";font-size: 12px;color: #039;text-decoration:underline;}")
        #P BLUE
        self.vContent.append(" #descriptionBlue \n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 14px;\n")
        self.vContent.append("   text-align: left;\n")
        self.vContent.append("   color: #039;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        #P GREEN
        self.vContent.append(" #descriptionGreen \n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 14px;\n")
        self.vContent.append("   text-align: left;\n")
        self.vContent.append("   color: #009900;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        #P SUBTITLE BLUE
        self.vContent.append(" #subtitleBlue \n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 16px;\n")
        self.vContent.append("   font-weight: bold;\n")
        self.vContent.append("   text-decoration:underline;\n")
        self.vContent.append("   text-align: left;\n")
        self.vContent.append("   color: #039;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")
        #P SUBTITLE GREEN
        self.vContent.append(" #subtitleGreen \n")
        self.vContent.append(" {\n")
        self.vContent.append("   font-family: \"Sans-Serif\";\n")
        self.vContent.append("   font-size: 16px;\n")
        self.vContent.append("   font-weight: bold;\n")
        self.vContent.append("   text-decoration:underline;\n")
        self.vContent.append("   text-align: left;\n")
        self.vContent.append("   color: #009900;\n")
        self.vContent.append("   border-collapse: collapse;\n")
        self.vContent.append(" }\n")

    def addParagraph(self,sContent,color,isSubtitle = False):
        '''Add Paragrah with content and given color '''
        if color == "green":
            if isSubtitle :
                self.vContent.append("  <P id=\"subtitleGreen\">\n")
            else :
                self.vContent.append("  <P id=\"descriptionGreen\">\n")
        else:
            if isSubtitle :    
                self.vContent.append("  <P id=\"subtitleBlue\">\n")
            else :
                self.vContent.append("  <P id=\"descriptionBlue\">\n")
 
        self.vContent.append(sContent)
        self.vContent.append("  </P>")

    def addTableHeader(self,vFields,color):
        '''Add table header to a given vector with vFields vector of names'''
        if color == "green":
            self.vContent.append("  <TABLE id=\"green\">\n")
        else:
            self.vContent.append("  <TABLE id=\"hor-zebra\">\n")

        self.vContent.append("   <TR>\n")
	
        for field in vFields:
            self.vContent.append("    <TH scope=\"col\">" + field + "</TH>\n")
	
        self.vContent.append("   </TR>\n")

    def addHtmlReportTableContent(self,vFields, isOdd, isImage = False):
        ''' Add Formated Table Contents from a vector of fields 
            isOdd: Boolean to indicate if the row is odd or pair 
            isImage: Boolean in case the content to add is an image'''
        #Check row pairing
        if isOdd == True: 
           self.vContent.append("   <TR class=\"odd\">\n")
        else:
           self.vContent.append("   <TR>\n")
	
        for field in vFields:
            if isImage == True:
                self.vContent.append("    <TD> <img src=\"" + field + "\" alt=\"" + field + "\" </TD>\n")                
            else:
                self.vContent.append("     <TD> " + field + " </TD>")
	
        self.vContent.append("   </TR>\n")

    def closeHtmlReportTableContent(self):
        ''' Close Formated Table '''
        self.vContent.append("  </TABLE>\n")
        self.vContent.append("\n")

    def sectionTitle(self, typeSection):
        ''' Creates a linker anchor for the section to be open, it also set section name
            To be defined in the Child Class        
        '''
        pass

    def linkToMainMenu(self):
        ''' Adds link to come back to main menu '''
        self.vContent.append("\n")
        self.vContent.append("<a class=\"link\" href=\"#linkers\"><CENTER><B>BACK</B></CENTER></a> <br>\n")
        self.vContent.append("\n")

   
    def createDescriptionTable(self,vDescription,color):
        ''' From a vector of string prints HTML paragraphs per each element '''
        for feature in vDescription:
            self.addParagraph(feature,color)
    
    def addDescriptionSectionTable(self,list_statements,color="blue"):
        ''' Add Description Section
        Arguments---
        list_statements : List of sentences to be printed
        color : Background color of the format table 
        '''
        self.sectionTitle("sampleDescription")
        self.addParagraph("","blue")
        self.addParagraph("","blue")
        self.createDescriptionTable(list_statements,color)
        self.addParagraph("","blue")
        self.addParagraph("","blue")
        self.addParagraph("","blue")
        self.addParagraph("","blue")        
            
    def addHtmlNewSection(self,data_container=None,section_title="",section_description="",color="",is_composition=False,is_stack=False,is_image=False):
        ''' Adds a new Html Docuemntation Section
        Arguments--
        data_container : Instance of class base on BaseStatistical, or a class composed by a set of instances based on a BaseSatistical class 
                         Contains vector of fields and values
        section_title : Title to create links to the section
        section_description : Description text of the section
        color : Muste be green or blue, background color of the data table
        is_composition : Boolean value to indication if the data container is a composition of diferent set of data
        is_stack : Boolean value to indicate if the data container is a stack (list) of elements of the same type
        is_image : True if the contents to add are images
        '''
        #1.Section Title         
        self.sectionTitle(section_title)
        #2.Description Area
        self.addParagraph(section_description,"blue")
        self.addParagraph("","blue")  
        #3.Table Creation
        if is_stack == True:
            self.addTableHeader(data_container[0].getListHeaders(),color) 
            is_odd = False
            for partial in data_container:
                 self.addHtmlReportTableContent(vFields = partial.getListValues(),isOdd=True,isImage=is_image)
                 is_odd = not is_odd
        elif is_composition == True:
            self.addTableHeader(data_container.getVectorOfFields(),color) 
            self.addHtmlReportTableContent(vFields = data_container.getGlobalStats(),isOdd=True,isImage=is_image)
            is_odd = False
            for partial in data_container.getListPartialStats():
                self.addHtmlReportTableContent(vFields=partial,isOdd=is_odd,isImage=is_image)
                is_odd = not is_odd
        else:
            self.addTableHeader(data_container.getListHeaders(),color)   
            self.addHtmlReportTableContent(vFields=data_container.getListValues(),isOdd=True,isImage=is_image)
                 
        self.closeHtmlReportTableContent()   
        self.linkToMainMenu()

    def saveDocument(self,file_name):
        '''Save the document to a given name file
        file_name: File path to store the document
        '''
        with open(file_name, 'w') as fileDocument:
            for line in self.vContent:
                fileDocument.write(line+ '\n')

    def createImagesTable(self,matrix,color):
       '''Creates a table of two columns where odd rows are titles and pair rows are images'''
       if color == "green":
           self.vContent.append("  <TABLE id=\"green\">\n")
       else:
           self.vContent.append("  <TABLE id=\"hor-zebra\">\n")

       isOdd = True
       for fields in matrix:
           if isOdd == True: 
               self.vContent.append("   <TR class=\"odd\">\n")
               for field in fields:
                   self.vContent.append("    <TH scope=\"col\">" + field + "</TH>\n")
           else:
               self.vContent.append("   <TR>\n")
               for field in fields:
                   self.vContent.append("    <TD> <img src=\"" + field + "\" alt=\"" + field + "\" </TD>\n")

           self.vContent.append("   </TR>\n")
           isOdd = not isOdd
        
       #CLOSE TABLE
       self.vContent.append("  </TABLE>\n")
       self.vContent.append("\n")


class BuildHtml(BaseDocumentHtml):

    def titleDocument(self):
        '''Title for the HTML Document'''
        self.vContent.append("  <H1 id=\"title\"> <U> CNV PIPELINE REPORT </U> </H1>\n")
    
    def linkersTable(self):
        ''' Links table section '''
        self.vContent.append("<a id=\"linkers\"></a>\n")

        self.vContent.append("<table id=\"linksTable-b\" >\n")
        self.vContent.append("    <tbody>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#sampleDescription\">Sample Description</a> </td>\n")
        self.vContent.append("        </tr>\n") 
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#mappingStats\">Mapping Stats</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#canavarOutput\">mrCaNaVar Output</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#cnDistribution\">Copy Number Distribution</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#cnPlots\">Copy Number Distribution Plot</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#normCnDistribution\">Normalized Copy Number Distribution</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#m1Duplications\">Method 1 Duplications</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#m2Duplications\">Method 2 Duplications</a> </td>\n")
        self.vContent.append("        </tr>\n")
        self.vContent.append("    </tbody>\n")
        self.vContent.append("</table>\n")

        self.vContent.append(" <br> \n")

   
    def sectionTitle(self, typeSection):
        ''' Creates a linker anchor for the section to be open, it also set section name'''
	
        linker = ""
        sectionName = ""
        
        if typeSection == "sampleDescription":
            linker += "id=\"sampleDescription\""
            sectionName += "Sample"
        elif typeSection == "mappingStats":
            linker += "id=\"mappingStats\""
            sectionName += "Mapping Stats"
        elif typeSection == "canavarOutput":
            linker += "id=\"canavarOutput\""
            sectionName += "mrCaNaVar Output"
        elif typeSection == "cnDistribution":
            linker += "id=\"cnDistribution\""
            sectionName += "Copy Number Distribution"
        elif typeSection == "cnPlots":
            linker += "id=\"cnPlots\""
            sectionName += "Copy Number Distribution Plot"
        elif typeSection == "normCnDistribution":
            linker += "id=\"normCnDistribution\""
            sectionName += "Normalized Copy Number Distribution"
        elif typeSection == "m1Duplications":
            linker += "id=\"m1Duplications\""
            sectionName += "Method 1 Duplications"
        elif typeSection == "m2Duplications":
            linker += "id=\"m2Duplications\""
            sectionName += "Method 2 Duplications"

        self.vContent.append("<a "+ linker +"></a>\n")
        self.vContent.append("<H1 id=\"section\"> " + sectionName + " </H1>\n")
        
class ReportBamHtml(BaseDocumentHtml):

    def titleDocument(self):
        '''Title for the HTML Document'''
        self.vContent.append("  <H1 id=\"title\"> <U> BAM ALIGNMENT REPORT </U> </H1>\n")
    
    def linkersTable(self):
        ''' Links table section '''
        self.vContent.append("<a id=\"linkers\"></a>\n")

        self.vContent.append("<table id=\"linksTable-b\" >\n")
        self.vContent.append("    <tbody>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#bamAlignment\">Bam Alignment</a> </td>\n")
        self.vContent.append("        </tr>\n") 
        self.vContent.append("    </tbody>\n")
        self.vContent.append("</table>\n")

        self.vContent.append(" <br> \n")

   
    def sectionTitle(self, typeSection):
        ''' Creates a linker anchor for the section to be open, it also set section name'''
	
        linker = ""
        sectionName = ""
        
        if typeSection == "bamAlignment":
            linker += "id=\"bamAlignment\""
            sectionName += "Bam Alignment Stats"
        
        self.vContent.append("<a "+ linker +"></a>\n")
        self.vContent.append("<H1 id=\"section\"> " + sectionName + " </H1>\n")

class ReportRmDupHtml(BaseDocumentHtml):

    def titleDocument(self):
        '''Title for the HTML Document'''
        self.vContent.append("  <H1 id=\"title\"> <U> REMOVE PCR DUPLICATES REPORT </U> </H1>\n")
    
    def linkersTable(self):
        ''' Links table section '''
        self.vContent.append("<a id=\"linkers\"></a>\n")

        self.vContent.append("<table id=\"linksTable-b\" >\n")
        self.vContent.append("    <tbody>\n")
        self.vContent.append("        <tr>\n")
        self.vContent.append("            <td> <a class=\"link\" href=\"#rmDup\">Remove PCR Duplicates</a> </td>\n")
        self.vContent.append("        </tr>\n") 
        self.vContent.append("    </tbody>\n")
        self.vContent.append("</table>\n")

        self.vContent.append(" <br> \n")

   
    def sectionTitle(self, typeSection):
        ''' Creates a linker anchor for the section to be open, it also set section name'''
	
        linker = ""
        sectionName = ""
        
        if typeSection == "rmDup":
            linker += "id=\"rmDup\""
            sectionName += "Remove PCR Duplicates"
        
        self.vContent.append("<a "+ linker +"></a>\n")
        self.vContent.append("<H1 id=\"section\"> " + sectionName + " </H1>\n")



class BaseStatistical(object):
    ''' General Statistical Methods and Steps '''
    def __init__(self, source_file=None, is_json = False):
        self.source_file = source_file
        self.is_json = is_json
        self.fieldValue = {}
        self.fieldsList = []
        self.valuesList = []
        self.allValues = {}
 
        self.parseInput()
        self.buildListHeaders()
        self.buildListValues()
        
        self.is_virgin = False #Control case were the instance is set to 0 for calculating total values    
        
    def parseInput(self):
        '''Implement this method to Parse input source'''
        pass
    
    def buildListHeaders(self):
        '''Implement this method to get a list of title fields'''
        pass
    
    def buildListValues(self):
        '''Implement this method to get a list of value fields'''
        pass
    
    def getListHeaders(self):
        '''Returns a sorted list of Headers'''
        return self.fieldsList
        
    def getListValues(self):
        '''Returns a sorted list of Values'''
        return self.valuesList
        
    def getHeaderValues(self):
        '''Implement this method to return a dictionary of fields and its values'''
        pass

class MappingJsonStats(BaseStatistical):
    '''Mapping Stats'''
                
    def parseInput(self):
        '''Parse input source'''
        if self.source_file == None:
            self.fieldValue["nameFiles"] = ""
            self.fieldValue["numberReads"] = 0
            self.fieldValue["numberTemplates"] = 0
            self.fieldValue["readLengthMin"] = 100
            self.fieldValue["readLengthAvg"] = 0
            self.fieldValue["readLengthMax"] = 0
            self.fieldValue["templatesMapped"] = 0
            self.fieldValue["readsMapped"] = 0
            self.fieldValue["perCenMapped"] = 0.0
            self.fieldValue["readsMappedLengthMin"] = 100
            self.fieldValue["readsMappedLengthAvg"] = 0
            self.fieldValue["readsMappedLengthMax"] = 0
            self.fieldValue["numBases"] = 0
            self.fieldValue["numBasesAligned"] = 0
            self.fieldValue["BasesA"] = 0
            self.fieldValue["BasesC"] = 0
            self.fieldValue["BasesG"] = 0
            self.fieldValue["BasesT"] = 0
            self.fieldValue["BasesN"] = 0
            self.is_virgin = True
        else:
            of = None
            if isinstance(self.source_file, basestring):
                of = open(self.source_file, 'r')

            data = json.load(of)
            of.close()         
          
            self.fieldValue["nameFiles"] = os.path.basename(self.source_file)
            self.fieldValue["numberReads"] = data["general"]["num_reads"]
            self.fieldValue["numberTemplates"] = data["general"]["num_templates"]
            self.fieldValue["readLengthMin"] = data["general"]["read_lenght_min"]
            self.fieldValue["readLengthAvg"] = data["general"]["read_lenght_avg"]
            self.fieldValue["readLengthMax"] = data["general"]["read_lenght_max"]
            self.fieldValue["templatesMapped"] = data["general"]["templates_mapped"]
            self.fieldValue["readsMapped"] = data["general"]["reads_mapped"]
            self.fieldValue["perCenMapped"] = (float(self.fieldValue["readsMapped"]) / self.fieldValue["numberReads"] ) * 100  
            self.fieldValue["readsMappedLengthMin"] = data["general"]["reads_mapped_length_min"]
            self.fieldValue["readsMappedLengthAvg"] = data["general"]["reads_mapped_length_avg"]
            self.fieldValue["readsMappedLengthMax"] = data["general"]["reads_mapped_length_max"]
            self.fieldValue["numBases"] = data["general"]["num_bases"]
            self.fieldValue["numBasesAligned"] = data["general"]["num_bases_aligned"]
            self.fieldValue["BasesA"] = data["general"]["bases_prop"]["A"]
            self.fieldValue["BasesC"] = data["general"]["bases_prop"]["C"]
            self.fieldValue["BasesG"] = data["general"]["bases_prop"]["G"]
            self.fieldValue["BasesT"] = data["general"]["bases_prop"]["T"]
            self.fieldValue["BasesN"] = data["general"]["bases_prop"]["N"]
            
            
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("")
        self.fieldsList.append("Number of Reads")
        self.fieldsList.append("Number of Templates")
        self.fieldsList.append("Read Len. Min.")
        self.fieldsList.append("Read Len. Avg.")
        self.fieldsList.append("Read Len. Max.")
        self.fieldsList.append("Templates Mapped")
        self.fieldsList.append("Reads Mapped")
        self.fieldsList.append("% Mapped")
        self.fieldsList.append("Reads Mapped Len. Min.")
        self.fieldsList.append("Reads Mapped Len. Avg.")
        self.fieldsList.append("Reads Mapped Len. Max.")
        self.fieldsList.append("Number of Bases")
        self.fieldsList.append("Number of Bases Aligned")
        self.fieldsList.append("A")
        self.fieldsList.append("C")
        self.fieldsList.append("G")
        self.fieldsList.append("T")
        self.fieldsList.append("N")
        
    def buildListValues(self):
        '''Returns a sorted list of Values'''
        self.valuesList.append(str(self.fieldValue["nameFiles"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["numberReads"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["numberTemplates"]))
        self.valuesList.append(str(self.fieldValue["readLengthMin"]))
        self.valuesList.append(str(self.fieldValue["readLengthAvg"]))
        self.valuesList.append(str(self.fieldValue["readLengthMax"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["templatesMapped"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["readsMapped"]))
        self.valuesList.append("{:.2f}".format(self.fieldValue["perCenMapped"]))
        self.valuesList.append(str(self.fieldValue["readsMappedLengthMin"]))
        self.valuesList.append(str(self.fieldValue["readsMappedLengthAvg"]))
        self.valuesList.append(str(self.fieldValue["readsMappedLengthMax"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["numBases"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["numBasesAligned"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["BasesA"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["BasesC"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["BasesG"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["BasesT"]))
        self.valuesList.append('{:,}'.format(self.fieldValue["BasesN"]))
   
            
    def __iadd__(self, other):
        ''' Overload iadd operator to sum two Mapping Json Stats Objects '''
        self.fieldValue["numberReads"] += other.fieldValue["numberReads"]
        self.fieldValue["numberTemplates"] += other.fieldValue["numberTemplates"]

        if  self.is_virgin:
            self.fieldValue["readLengthMin"] = other.fieldValue["readLengthMin"]
            self.fieldValue["readLengthAvg"] = other.fieldValue["readLengthAvg"]
            self.fieldValue["readLengthMax"] = other.fieldValue["readLengthMax"]
            self.is_virgin = False
        else:
            self.fieldValue["readLengthMin"] = min ([self.fieldValue["readLengthMin"], other.fieldValue["readLengthMin"]])
            self.fieldValue["readLengthAvg"] = numpy.mean([self.fieldValue["readLengthAvg"], other.fieldValue["readLengthAvg"]])
            self.fieldValue["readLengthMax"] = max ([self.fieldValue["readLengthMax"],other.fieldValue["readLengthMax"]])

        self.fieldValue["templatesMapped"] += other.fieldValue["templatesMapped"]
        self.fieldValue["readsMapped"] += other.fieldValue["readsMapped"]
        self.fieldValue["perCenMapped"] = (float(self.fieldValue["readsMapped"]) / self.fieldValue["numberReads"] ) * 100

        if  self.is_virgin:
            self.fieldValue["readsMappedLengthMin"] = other.fieldValue["readsMappedLengthMin"]
            self.fieldValue["readsMappedLengthAvg"] = other.fieldValue["readsMappedLengthAvg"]
            self.fieldValue["readsMappedLengthMax"] = other.fieldValue["readsMappedLengthMax"]
            self.is_virgin = False
        else:
            self.fieldValue["readsMappedLengthMin"] = min ([self.fieldValue["readsMappedLengthMin"], other.fieldValue["readsMappedLengthMin"]])
            self.fieldValue["readsMappedLengthAvg"] = numpy.mean([self.fieldValue["readsMappedLengthAvg"], other.fieldValue["readsMappedLengthAvg"]])
            self.fieldValue["readsMappedLengthMax"] = max ([self.fieldValue["readsMappedLengthMax"],other.fieldValue["readsMappedLengthMax"]])
        
        self.fieldValue["numBases"] += other.fieldValue["numBases"]
        self.fieldValue["numBasesAligned"] += other.fieldValue["numBasesAligned"]
        self.fieldValue["BasesA"] += other.fieldValue["BasesA"]
        self.fieldValue["BasesC"] += other.fieldValue["BasesC"]
        self.fieldValue["BasesG"] += other.fieldValue["BasesG"]
        self.fieldValue["BasesT"] += other.fieldValue["BasesT"]
        self.fieldValue["BasesN"] += other.fieldValue["BasesN"]
        
        self.valuesList = []
        self.buildListValues()
        
        return self
        


class GlobalMappingJsonStats(object):
    '''Mapping Stats for all mapped files'''
    def __init__(self,listFilesJson = None):
        '''Add new MappingJsonStats object for each Json file'''
        self.mapStatsList = []
        for jsonFile in listFilesJson:
            self.mapStatsList.append(MappingJsonStats(jsonFile,True))
 
        self.totalStats = MappingJsonStats()
        self.buildGlobalStats()
                         
    def getVectorOfFields(self):
        '''Returns a list of fields name to be printed as HTML table headers '''
        return self.mapStatsList[0].getListHeaders()
        
    def buildGlobalStats(self):
        ''' Returns a vector of global stats '''
        for stats in self.mapStatsList:
           self.totalStats += stats
        
    def getGlobalStats(self):
        ''' Returns a vector of global stats '''
        return self.totalStats.getListValues()
            
    def getListPartialStats(self):
        ''' Returns a list of stats map instances'''
        partialStats = []
        for stats in self.mapStatsList:
           partialStats.append(stats.getListValues())
        
        return partialStats
        
    def getAllDataValues(self):
        '''Return a dictionary of all distionaries generated per each mapping json file'''
        allDataValues = {}
        globalStats = {}
        globalStats["total"] = self.totalStats.fieldValue
        allDataValues.update(globalStats)
        
        for stats in self.mapStatsList:
            currentStats = {}
            currentStats [stats.fieldValue["nameFiles"]] = stats.fieldValue
            allDataValues.update(currentStats)
        
        mappingStats = {}
        mappingStats ["MappingStats"] = allDataValues
        return mappingStats

class PcrDuplicatesStats(BaseStatistical):
    '''Class for manage RM Duplicates Picard Tools metrics stats'''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Name")
        self.fieldsList.append("Unpaired Reads Examined")
        self.fieldsList.append("Read Pairs Examined")
        self.fieldsList.append("Unmapped Reads")
        self.fieldsList.append("Unpaired Read Duplicates")
        self.fieldsList.append("Read Pair Duplicates")
        self.fieldsList.append("Read Pair Optical Duplicates")
        self.fieldsList.append("Duplication Coefficient [0-1]")
        self.fieldsList.append("Estimated Library Size")
        
    def parseInput(self):
        '''Parse input source'''
        if self.source_file == None:
            self.fieldValue['Name'] = ""
            self.fieldValue['UnpairedReadsExamined'] = 0
            self.fieldValue['ReadsPairsExamined'] = 0
            self.fieldValue['UnmappedReads'] = 0
            self.fieldValue['UnpairedReadDuplicates'] = 0
            self.fieldValue['ReadPairDuplicates'] = 0
            self.fieldValue['ReadPairOpticalDuplicates'] = 0
            self.fieldValue['DuplicationCoefficient'] = 0
            self.fieldValue['EstimatedLibrarySize'] = 0
            self.is_virgin = True
        else:
            processLine = False        
            with open(self.source_file, "r") as fileLog:
                for line in fileLog:
                    vFields = re.split('\t',line.rstrip('\n'))
                    if processLine == False and len(vFields) > 1:
                        if vFields[0] == "LIBRARY":
                            processLine = True
                    elif processLine == True:
                        self.fieldValue['Name'] = os.path.basename(self.source_file)
                        self.fieldValue['UnpairedReadsExamined'] = long(vFields[1])
                        self.fieldValue['ReadsPairsExamined'] = long(vFields[2])
                        self.fieldValue['UnmappedReads'] = long(vFields[3])
                        self.fieldValue['UnpairedReadDuplicates'] = long(vFields[4])
                        self.fieldValue['ReadPairDuplicates'] = long(vFields[5])
                        self.fieldValue['ReadPairOpticalDuplicates'] = long(vFields[6])
                        self.fieldValue['DuplicationCoefficient'] = float(vFields[7])
                        if vFields[8] != "":
                            self.fieldValue['EstimatedLibrarySize'] = long(vFields[8])
                        else:
                            self.fieldValue['EstimatedLibrarySize'] = ""
                        break
                    
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append(self.fieldValue['Name'])
        self.valuesList.append("{:,}".format(self.fieldValue['UnpairedReadsExamined']))
        self.valuesList.append("{:,}".format(self.fieldValue['ReadsPairsExamined']))
        self.valuesList.append("{:,}".format(self.fieldValue['UnmappedReads']))
        self.valuesList.append("{:,}".format(self.fieldValue['UnpairedReadDuplicates']))
        self.valuesList.append("{:,}".format(self.fieldValue['ReadPairDuplicates']))
        self.valuesList.append("{:,}".format(self.fieldValue['ReadPairOpticalDuplicates']))
        self.valuesList.append("{:.5f}".format(self.fieldValue['DuplicationCoefficient']))
        if  self.fieldValue["EstimatedLibrarySize"] != "": 
            self.valuesList.append("{:,}".format(self.fieldValue['EstimatedLibrarySize']))
        else:
            self.valuesList.append("")
        
    def __iadd__(self, other):
        ''' Overload iadd operator to sum two PcrDuplicatesStats Objects '''
        self.fieldValue["UnpairedReadsExamined"] += other.fieldValue["UnpairedReadsExamined"]
        self.fieldValue["ReadsPairsExamined"] += other.fieldValue["ReadsPairsExamined"]
        self.fieldValue["UnmappedReads"] += other.fieldValue["UnmappedReads"]
        self.fieldValue["UnpairedReadDuplicates"] += other.fieldValue["UnpairedReadDuplicates"]
        self.fieldValue["ReadPairDuplicates"] += other.fieldValue["ReadPairDuplicates"]
        self.fieldValue["ReadPairOpticalDuplicates"] += other.fieldValue["ReadPairOpticalDuplicates"]
        if self.is_virgin == True:
            self.fieldValue["DuplicationCoefficient"] = other.fieldValue["DuplicationCoefficient"]
            self.is_virgin = False
        else:
            self.fieldValue["DuplicationCoefficient"] = numpy.mean([self.fieldValue["DuplicationCoefficient"], other.fieldValue["DuplicationCoefficient"]])
        if  other.fieldValue["EstimatedLibrarySize"] != "":       
            self.fieldValue["EstimatedLibrarySize"] += other.fieldValue["EstimatedLibrarySize"]
        else:
            self.fieldValue["EstimatedLibrarySize"] = ""
        
        self.valuesList = []
        self.buildListValues()
        
        return self
    
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["rmDuplicatesStats"] = self.fieldValue 
        return self.allValues

class GlobalDuplicatesStats(object):
    '''Mangement of Duplicates stats, per each mapped file'''
   
    def __init__(self,listFilesMetrics = None):
        '''Add new PcrDuplicatesStats object for each metric file'''
        self.dupStatsList = []
        
        for metricsFile in listFilesMetrics:
            self.dupStatsList.append(PcrDuplicatesStats(metricsFile,True))
 
        self.totalStats = PcrDuplicatesStats()
        self.buildGlobalStats()
        
        
    def getVectorOfFields(self):
        '''Returns a list of fields name to be printed as HTML table headers '''
        return self.dupStatsList[0].getListHeaders()
        
    def buildGlobalStats(self):
        ''' Returns a vector of global stats '''
        for stats in self.dupStatsList:
            self.totalStats += stats
        
    def getGlobalStats(self):
        ''' Returns a vector of global stats '''
        return self.totalStats.getListValues()
            
    def getListPartialStats(self):
        ''' Returns a list of stats map instances'''
        partialStats = []
        for stats in self.dupStatsList:
           partialStats.append(stats.getListValues())
        
        return partialStats
        
    def getAllDataValues(self):
        '''Return a dictionary of all distionaries generated per each mapping json file'''
        allDataValues = {}
        globalStats = {}
        globalStats["total"] = self.totalStats.fieldValue
        allDataValues.update(globalStats)
        
        for stats in self.dupStatsList:
            currentStats = {}
            currentStats [stats.fieldValue["Name"]] = stats.fieldValue
            allDataValues.update(currentStats)
        
        mappingStats = {}
        mappingStats ["DuplicatesStats"] = allDataValues
        return mappingStats
        
class BamMappingStats(BaseStatistical):
    '''Class for manage Bam Mapping Stats'''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Sample")
        self.fieldsList.append("Total Reads")
        self.fieldsList.append("Aligned Reads")
        self.fieldsList.append("% Alignment")
        self.fieldsList.append("Mismatch Rate")
        self.fieldsList.append("InDel Rate")
        self.fieldsList.append("Strand Balance")
        self.fieldsList.append("Aligned Bases")
        
    def parseInput(self):
        '''Parse input source '''
        self.fieldValue["sample"] = os.path.basename(self.source_file)
        with open(self.source_file, "r") as fileLog:
            for line in fileLog:
                vFields = re.split('\s+',line.rstrip('\n'))
                
                if len(vFields) > 18:
                    if vFields[0] == "PAIR" or vFields[0] == "UNPAIRED":                        
                        #total_read
                        self.fieldValue['total_reads'] = long(vFields[1])
                        #Aligned reads
                        self.fieldValue['aligned_reads'] = long(vFields[5])
                        #% Aligned
                        self.fieldValue["percen_aligned"] = (float(self.fieldValue["aligned_reads"]) / self.fieldValue["total_reads"] ) * 100
                        #Mismatch rate
                        self.fieldValue['mismatch_rate'] = vFields[12]
                        #Indel Rate
                        self.fieldValue['indel_rate'] = vFields[14]
                        #Strand Balance
                        self.fieldValue['strand_balance'] = vFields[19]
                        #Aligned Bases
                        self.fieldValue['aligned_bases'] = long(vFields[7])
                        
                        
    def buildListValues(self):
        '''Returns a sorted list of Values'''
        #Sample
        self.valuesList.append(self.fieldValue["sample"])
        #Total Reads
        self.valuesList.append('{:,}'.format(self.fieldValue["total_reads"]))
        #Aligned Reads
        self.valuesList.append('{:,}'.format(self.fieldValue["aligned_reads"]))
        #% Alignment
        self.valuesList.append("{:.2f}".format(self.fieldValue["percen_aligned"]))
        #Mismatch Rate
        self.valuesList.append(self.fieldValue["mismatch_rate"])
        #InDel Rate
        self.valuesList.append(self.fieldValue["indel_rate"])
        #Strand Balance
        self.valuesList.append(self.fieldValue["strand_balance"])
        #Aligned Bases        
        self.valuesList.append('{:,}'.format(self.fieldValue["aligned_bases"]))
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["bamMappingStats"] = self.fieldValue 
        return self.allValues
        
 
class MrCanavarStats(BaseStatistical):
    '''Class for manage mrCanvar stats'''

    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("LW Read Depth Avg.")
        self.fieldsList.append("LW Read Depth StDev.")
        self.fieldsList.append("SW Read Depth Avg.")
        self.fieldsList.append("SW Read Depth StDev.")
        self.fieldsList.append("CW Read Depth Avg.")
        self.fieldsList.append("CW Read Depth StDev.")
        
    def parseInput(self):
        '''Parse input source'''
        with open(self.source_file, "r") as fileLog:
            for line in fileLog:
                vFields = re.split('\s+',line.rstrip('\n'))
                if line.startswith('CW Average Read Depth:'):
                    self.fieldValue['CW'] = [float(vFields [4][:-1]), float(vFields [7])]
                elif line.startswith('LW Average Read Depth:'):
                    self.fieldValue['LW'] = [float(vFields [4][:-1]), float(vFields [7])]
                elif line.startswith('SW Average Read Depth:'):
                    self.fieldValue['SW'] = [float(vFields [4][:-1]), float(vFields [7])]
        
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append("{:,.2f}".format(self.fieldValue["LW"][0]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["LW"][1]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["SW"][0]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["SW"][1]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["CW"][0]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["CW"][1]))
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["mrcanavar"] = self.fieldValue 
        return self.allValues
 
class CopyNumberDistributionControlRegions(BaseStatistical):   
    ''' Statistical Values for copy number distribution in Control Regions '''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Mean")
        self.fieldsList.append("Median")
        self.fieldsList.append("Min")
        self.fieldsList.append("Max")
        self.fieldsList.append("Standard Dev.")
        
    def parseInput(self):
        '''Parse input source'''
        readLine = 0
        with open(self.source_file, "r") as fileResults:
            for line in fileResults:
                readLine = readLine + 1
                if readLine == 2:
                    vFields = re.split('\s+',line.rstrip('\n'))   
                    self.fieldValue['Mean'] = float(vFields[2])
                    self.fieldValue['Median'] = float(vFields[3])
                    self.fieldValue['Min'] = float(vFields[1])
                    self.fieldValue['Max'] = float(vFields[5])
                    self.fieldValue['StDev'] = float(vFields[4]) 
            
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Mean"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Median"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Min"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Max"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["StDev"]))
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["cnDistribution"] = self.fieldValue 
        return self.allValues        
        
class ControlRegionsPlot(BaseStatistical):
    ''' Plot of the Control Regions Distribution '''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Control Regions Distribution")
    
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append(self.source_file)  
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["cnDistributionPlot"] = self.source_file 
        return self.allValues        

class NormalizedCopyNumberDistribtuion(BaseStatistical):
    ''' Stats normalized of copy number '''
 
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Mean")
        self.fieldsList.append("Standard Deviation (Excluding windows with 1% most extreme copy number)")
        self.fieldsList.append("Windows Excluded")
        self.fieldsList.append("Coefficient windows exluded")
        self.fieldsList.append("Gain CutOff")
        self.fieldsList.append("Loss CutOff")
        
    def parseInput(self):
        '''Parse input source'''
        readLine = 0
        with open(self.source_file, "r") as fileResults:
            for line in fileResults:
                readLine = readLine + 1
                if readLine == 2:
                    vFields = re.split('\s+',line.rstrip('\n'))   
                    self.fieldValue['Mean'] = float(vFields[1])
                    self.fieldValue['StDev'] = float(vFields[2])
                    self.fieldValue['Windows'] = int(vFields[3])
                    self.fieldValue['Coeff'] = float(vFields[4])
                    self.fieldValue['Gain'] = float(vFields[5]) 
                    self.fieldValue['Loss'] = float(vFields[6])
                
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Mean"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["StDev"]))
        self.valuesList.append("{:,}".format(self.fieldValue["Windows"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Coeff"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Gain"]))
        self.valuesList.append("{:,.2f}".format(self.fieldValue["Loss"]))
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["normalizedCnDistribution"] = self.fieldValue 
        return self.allValues
              
class BedDuplicationsStats(BaseStatistical):
    ''' Duplications BED file '''
    
    def buildListHeaders(self):
        '''Implementation of list of headers construction'''
        self.fieldsList.append("Windows Called")
        self.fieldsList.append("Bases Called")
        
    def parseInput(self):
        '''Parse input source'''        
        totalLines = 0
        totalBases = 0
        with open(self.source_file, "r") as bedData:
            for line in bedData:
                totalLines = totalLines + 1
                vFields = line.rstrip('\n').split("\t")
                totalBases = totalBases + (long(vFields[2]) - long(vFields[1]))

        self.fieldValue['Windows'] = totalLines
        self.fieldValue['Bases'] = totalBases
 
    def buildListValues(self):
        '''Get a list of value fields'''
        self.valuesList.append("{:,}".format(self.fieldValue["Windows"]))
        self.valuesList.append("{:,}".format(self.fieldValue["Bases"]))
        
    def getHeaderValues(self):
        '''Return a dictionary of fields and its values'''
        self.allValues ["Duplications_"+os.path.basename(self.source_file)] = self.fieldValue 
        return self.allValues
 
 
def create_report(html_file,json_file,mapping_stats_files, mrcanavar_log_file, 
                  control_regions_distribution, control_regions_plot, 
                  cutoffs_file, duplications_M1, duplications_M2,sample_description=None):
    """ Generate HTML Report from mapping stat, mrCanavar log, control regions distribution, control regions plot and cutoffs file 
    
    Parameters
    ----------
    html_file: html file to store the document 
    json_file: json file to store the document
    mapping_stats_files: json mapping stats file
    mrcanavar_log_file: mrcanavar log file
    control_regions_distribution: control regions distribution file
    control_regions_plot: control regions plot file
    cutoffs_file: cutoffs normalization in control regions file
    duplications_M1: bed file of duplications of Method One
    duplications_M2: bed file of duplications of Method Two
    sample_description = Sample Description
    """    
    
    #Global Mapping Stats    
    mappingStatistics = GlobalMappingJsonStats(listFilesJson = mapping_stats_files)
    #MrCanavar stats
    canavarStatistics = MrCanavarStats(source_file=mrcanavar_log_file,is_json = False)
    #Copy Number Distribution on Control Regions stats
    cnControlRegions = CopyNumberDistributionControlRegions(source_file=control_regions_distribution,is_json = False)
    #Control Regions Plot
    cnControlRegionsPlot = ControlRegionsPlot(source_file=control_regions_plot,is_json = False)
    #Normalized Copy Number Distribution
    normCNControlRegions = NormalizedCopyNumberDistribtuion(source_file=cutoffs_file,is_json = False)
   
    if duplications_M1 is not None:
        #Method 1 Duplications
        m1Duplications = BedDuplicationsStats(source_file=duplications_M1,is_json = False)
        
    if duplications_M2 is not None:     
        #Method 2 Duplications
        m2Duplications = BedDuplicationsStats(source_file=duplications_M2,is_json = False)
    
    #1. HTML REPORT CONSTRUCTION
    #1.1 Header HTML
    vHtmlContent = []
    htmlManager = BuildHtml(vHtmlContent)
    htmlManager.addHtmlReportHeader()
    
    #1.1 Sample Description
    if sample_description != None:
        htmlManager.addDescriptionSectionTable([sample_description])

    #1.3 Section Mapping Output
    htmlManager.addHtmlNewSection(mappingStatistics,"mappingStats","GEM Mapping stats","green",is_composition=True)   
    
    #1.4 Section MrCanavar output
    htmlManager.addHtmlNewSection(canavarStatistics,"canavarOutput","MrCanavar Output","blue")
    
    #1.5 Section Copy Number Distribution in Control Regions 
    htmlManager.addHtmlNewSection(cnControlRegions,"cnDistribution","Copy Number Distribution in Control Regions","green")
    
    #1.6 Section Copy Number Distribution Plots
    htmlManager.addHtmlNewSection(data_container=cnControlRegionsPlot,section_title="cnPlots",
                                  section_description="Copy Number Distribution Plots",color="blue",is_image=True)
    
    #1.7 Section Normalized Distribution in Control Regions    
    htmlManager.addHtmlNewSection(normCNControlRegions,"normCnDistribution","Normalized Copy Number Distribution in Control Regions","green")

    if duplications_M1 is not None:
        #1.8 Method One for calling Duplications
        htmlManager.addHtmlNewSection(m1Duplications,"m1Duplications","Method 1: Copy Number non overlapping windows of 1000bp, +3 St.Dev, Max 100 Copies","blue")
    
    if duplications_M2 is not None:
        #1.9 Method Two for calling Duplications
        htmlManager.addHtmlNewSection(m2Duplications,"m2Duplications","Method 2: 7 windows where at least 6 over +4 St.Dev, Regions greater than 10.000 bp","green")
    
    #1.10 Close Html Report
    htmlManager.closeHtmlReport()
    
    #2. SAVE HTML DOCUMENT
    htmlManager.saveDocument(html_file)
        
    #3. CREATE CASCADE STYLE SHEET
    vCSScontent = []
    cssManager = BuildHtml(vCSScontent)
    cssManager.buildStyleSheet()

    cssManager.saveDocument(os.path.dirname(os.path.abspath(html_file)) + "/style.css")
    
    #4.CREATE JSON
    jsonDataDocument = {}
    jsonDataDocument.update(mappingStatistics.getAllDataValues())
    jsonDataDocument.update(canavarStatistics.getHeaderValues())
    jsonDataDocument.update(cnControlRegions.getHeaderValues())
    jsonDataDocument.update(cnControlRegionsPlot.getHeaderValues())
    jsonDataDocument.update(normCNControlRegions.getHeaderValues())
    if duplications_M1 is not None:
        jsonDataDocument.update(m1Duplications.getHeaderValues())
    if duplications_M2 is not None:    
        jsonDataDocument.update(m2Duplications.getHeaderValues())
    
    with open(json_file, 'w') as of:
        json.dump(jsonDataDocument, of, indent=2)
  
        
def create_bam_report(metrics=None,html_file=None,json_file=None,sample_description=None):
    """ Generate HTML Report from picard tool callectAlginmentSummary Metrics
    
    Parameters
    ----------
    metrics: Metrics File from Picard Tools
    html_file: HTML document file to be generated
    json_file: JSON File to be generated
    sample_description: Sample Description
    """  
    #MrCanavar stats
    bamStats = BamMappingStats(source_file=metrics,is_json = False)
    
    #1. HTML REPORT CONSTRUCTION
    #1.1 Header HTML
    vHtmlContent = []
    htmlManager = ReportBamHtml(vHtmlContent)
    htmlManager.addHtmlReportHeader()
    
    #1.1 Sample Description
    if sample_description != None:
        htmlManager.addDescriptionSectionTable([sample_description])
        
    #1.2 Section Bam Report Stats
    htmlManager.addHtmlNewSection(bamStats,"bamAlignment","BAM Alignment File Stats","blue")
    
     #1.10 Close Html Report
    htmlManager.closeHtmlReport()
    
    #2. SAVE HTML DOCUMENT
    htmlManager.saveDocument(html_file)
        
    #3. CREATE CASCADE STYLE SHEET
    vCSScontent = []
    cssManager = BuildHtml(vCSScontent)
    cssManager.buildStyleSheet()

    cssManager.saveDocument(os.path.dirname(os.path.abspath(html_file)) + "/style.css")
    
    #4.CREATE JSON
    jsonDataDocument = {}
    jsonDataDocument.update(bamStats.getHeaderValues())
    
    with open(json_file, 'w') as of:
        json.dump(jsonDataDocument, of, indent=2)

def create_duplicates_report(metrics=None,html_file=None,json_file=None,sample_description=None):
    """ Generate HTML Report from picard tool callectAlginmentSummary Metrics
    
    Parameters
    ----------
    metrics: Metrics File from Picard Tools
    html_file: HTML document file to be generated
    json_file: JSON File to be generated
    sample_description: Sample Description
    """  
    rmStats = PcrDuplicatesStats(source_file=metrics,is_json = False)
    
     #1. HTML REPORT CONSTRUCTION
    #1.1 Header HTML
    vHtmlContent = []
    htmlManager = ReportRmDupHtml(vHtmlContent)
    htmlManager.addHtmlReportHeader()
    
    #1.1 Sample Description
    if sample_description != None:
        htmlManager.addDescriptionSectionTable([sample_description])
        
    #1.2 Section Bam Report Stats
    htmlManager.addHtmlNewSection(rmStats,"rmDup","Remove Duplicates Stats","blue")
    
    #1.10 Close Html Report
    htmlManager.closeHtmlReport()
    
    #2. SAVE HTML DOCUMENT
    htmlManager.saveDocument(html_file)
        
    #3. CREATE CASCADE STYLE SHEET
    vCSScontent = []
    cssManager = BuildHtml(vCSScontent)
    cssManager.buildStyleSheet()

    cssManager.saveDocument(os.path.dirname(os.path.abspath(html_file)) + "/style.css")
    
    #4.CREATE JSON
    jsonDataDocument = {}
    jsonDataDocument.update(rmStats.getHeaderValues())
    
    with open(json_file, 'w') as of:
        json.dump(jsonDataDocument, of, indent=2)