/*
 * main.cpp
 *
 *  Created on: 23 Mar, 2015
 *      Author: marcos
 */

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <sstream>
#include <string.h>


using namespace std;

static std::string int2String(int number){ std::stringstream ss; ss << number; return ss.str();}

/**
 * \brief Print Help File
 */
void showVersion()
{
    cout << "FASTQ Splitter C++ 0.2" << endl;
    cout << "Marcos Fernandez" << endl;
}

/**
 * \brief Show Help File
 */
void showHelp()
{
    cout << "Usage: FastqSplitter [options] <file> " << endl;
    cout << "    Options:" << endl;
    cout << "       --n-parts <N>        - Divide into <N> parts" << endl;
    cout << "       --version            - Show version." << endl;
    cout << "       --help               - Show help." << endl;
    cout << "       --prefix             - Output prefix for the files to be created. By default would be tha basename of the input file" << endl;
    cout << "       --gz                 - Is fastq Gzipped." << endl;
}


/**
 * \brief Get The number total of reads
 * \param inFileName File Name
 * \param totalNumberOfReads Total Number of reads
 * \param isGzip Is fasts gzipped
 * \return 0 in success otherwise 1
 */
int getTotalNumberReads (const string & inFileName,unsigned long int & totalNumberOfReads,bool isGzip)
{
	char header[4 * 256];
	char sequence[4 * 256];
	char plus[4 * 256];
	char qualities[4 * 256];
	char * token = "";

	unsigned long int lineCount = 0;

	totalNumberOfReads = 0;

	FILE *fp;

	if (!isGzip) fp = fopen(inFileName.c_str(), "r");
    else fp = (FILE *) gzopen(inFileName.c_str(), "r");

    if (fp == NULL)
    {
    	fprintf(stderr, "Sorry!! Unable to open file %s.", inFileName.c_str());
		return 1;
	}

    while (1)
    {
    	//1 Read entire line from the Fastq file
        if (isGzip)
    	{
        	gzgets((gzFile) fp, header, (4 * 256));
        	gzgets((gzFile) fp, sequence, (4 * 256));
        	gzgets((gzFile) fp, plus, (4 * 256));
        	gzgets((gzFile) fp, qualities, (4 * 256));
    	    if (gzeof((gzFile) fp)) break;
    	}
    	else
    	{
    		token = fgets(header, (4 * 256), fp);
    		token = fgets(sequence, (4 * 256), fp);
    		token = fgets(plus, (4 * 256), fp);
    		token = fgets(qualities, (4 * 256), fp);
    	    if (feof(fp)) break;
    	}

        //1.2 Chomp lines
        if (header[strlen(header) - 1] == '\n') header[strlen(header) - 1] = '\0';
        if (sequence[strlen(sequence) - 1] == '\n') sequence[strlen(sequence) - 1] = '\0';
        if (plus[strlen(plus) - 1] == '\n') plus[strlen(plus) - 1] = '\0';
        if (qualities[strlen(qualities) - 1] == '\n') qualities[strlen(qualities) - 1] = '\0';

        //1. READ HEADER
        lineCount ++;
        if (strcmp(header, "") == 0)
        {
        	cerr << "Error parsing line " << lineCount << " FASTQ entry does not start with @" << endl;
        	return 1;
        }
		//2.READ SEQUENCE
        lineCount ++;
        //3.READ PLUS
        lineCount ++;
        //4.READ QUALITIES
        lineCount ++;
        if(strlen(sequence) != strlen(qualities))
        {
             	cerr << "Warning: Misformatted FASTQ entry in input line " << lineCount << ": quality length (" << strlen(qualities);
               	cerr << ") differs from sequence length (" << strlen(sequence) << ")" << endl;
        }


        totalNumberOfReads ++;

        //Reset to control remaining
        strcpy(header, "");
        strcpy(sequence, "");
        strcpy(plus, "");
        strcpy(qualities, "");
    }

    //Remaining lines
    if (strcmp(header, "") != 0)
    {
    	//1. READ HEADER
    	lineCount ++;
    	//2.READ SEQUENCE
		lineCount ++;
		//3.READ PLUS
		lineCount ++;
		//4.READ QUALITIES
		lineCount ++;
		if(strlen(sequence) != strlen(qualities))
		{
				cerr << "Warning: Misformatted FASTQ entry in input line " << lineCount << ": quality length (" << strlen(qualities);
				cerr << ") differs from sequence length (" << strlen(sequence) << ")" << endl;
		}

		totalNumberOfReads ++;
    }


    //6. Close File
    if (isGzip) gzclose((gzFile)fp);
    else fclose(fp);

    return 0;
}

/**
 * \brief Print Starting Message
 * \param nParts Number of parts to divide the file
 * \param nPartsMeasured Number of parts to measure
 * \param measure Type of Measure
 */
void startingMessage(unsigned int nParts)
{
	cout << " => dividing into " << nParts << " part" ;
	if (nParts)
	{
		cout << "s";
	}
	else
	{
		cout << "";
	}
	cout << endl;
}
/**
 * \brief Split File
 * \param nParts Number of parts to divide the file
 * \param sOutputPrefix output prefix
 * \param isGzip Is input gzipped
 */
int splitFile(const string & inFileName,const string & sOutputPrefix,unsigned int nParts,bool isGzip)
{
	FILE * fout;

	char header[4 * 256];
	char sequence[4 * 256];
	char plus[4 * 256];
	char qualities[4 * 256];
    char * token = "";

	strcpy(header, "");
	strcpy(sequence, "");
	strcpy(plus, "");
	strcpy(qualities, "");

	string base = inFileName;
	if(!sOutputPrefix.empty()) base = sOutputPrefix;


	unsigned int totalMeasuredSize = 0;

	unsigned long int totalNumberOfReads = 0;

	if (getTotalNumberReads(inFileName,totalNumberOfReads,isGzip) != 0)
	{
		return 1;
	}
	cout << "Number of reads: " << totalNumberOfReads << endl;
	startingMessage(nParts);

	unsigned long int readsPerPart = round((float)((float) totalNumberOfReads / (float) nParts));

	FILE *fin;

	if (!isGzip) fin = fopen(inFileName.c_str(), "r");
	else fin = (FILE *) gzopen(inFileName.c_str(), "r");

	if (fin == NULL)
	{
		fprintf(stderr, "Sorry!! Unable to open file %s.", inFileName.c_str());
		return 1;
	}

	bool keepReading = true;

	unsigned long int remainingReads = totalNumberOfReads;
	unsigned long int currentWrote = 0;


	for (unsigned int part=1; part <= nParts; part++)
	{
		string part_file = base + ".part-" + int2String(part);

		fout=fopen(part_file.c_str(), "w");

		while (keepReading)
		{
			//1 Read entire line from the Fastq file
			if (isGzip)
			{
				gzgets((gzFile) fin, header, (4 * 256));
				gzgets((gzFile) fin, sequence, (4 * 256));
				gzgets((gzFile) fin, plus, (4 * 256));
				gzgets((gzFile) fin, qualities, (4 * 256));
				if (gzeof((gzFile) fin)) {keepReading = false; break;}
			}
			else
			{
				token = fgets(header, (4 * 256), fin);
				token = fgets(sequence, (4 * 256), fin);
				token = fgets(plus, (4 * 256), fin);
				token = fgets(qualities, (4 * 256), fin);
				if (gzeof((gzFile) fin)) {keepReading = false; break;}
			}

			//1.2 Chomp lines
	        if (header[strlen(header) - 1] == '\n') header[strlen(header) - 1] = '\0';
	        if (sequence[strlen(sequence) - 1] == '\n') sequence[strlen(sequence) - 1] = '\0';
			if (plus[strlen(plus) - 1] == '\n') plus[strlen(plus) - 1] = '\0';
			if (qualities[strlen(qualities) - 1] == '\n') qualities[strlen(qualities) - 1] = '\0';
			if (strcmp(header, "") != 0 && strcmp(qualities, "") != 0)
			{
				fwrite (header , sizeof(char), strlen(header), fout);
				fwrite ("\n" , sizeof(char), 1, fout);
				fwrite (sequence , sizeof(char), strlen(sequence), fout);
				fwrite ("\n" , sizeof(char), 1, fout);
				fwrite (plus , sizeof(char), strlen(plus), fout);
				fwrite ("\n" , sizeof(char), 1, fout);
				fwrite (qualities , sizeof(char), strlen(qualities), fout);
				fwrite ("\n" , sizeof(char), 1, fout);
			}

			currentWrote ++;
			remainingReads --;

			//Reset
			strcpy(header, "");
			strcpy(sequence, "");
			strcpy(plus, "");
			strcpy(qualities, "");

			if( (part < nParts) && (currentWrote < readsPerPart)) keepReading = true;
			else if( (part == nParts) && (remainingReads != 0))  keepReading = true;
			else keepReading = false;
		}

		currentWrote = 0;
		keepReading = true;


		fclose(fout);
	}


	//6. Close File
	if (isGzip) gzclose((gzFile)fin);
	else fclose(fin);

	return 0;
}


/**
 * \brief main function
 */
int main(int argc, char *argv[])
{
	bool bShowHelp = false;
	bool bShowVersion = false;
	bool bIsGzip = false;
	unsigned int nParts = 0;
	string sFileInputName = "";
	string sOutputPrefix = "";

    for (int i = 1; i < argc; i++)
    {
        if (string(argv[i]).compare("--help") == 0)
	    {
        	bShowHelp = true;
	    }
        else if (string(argv[i]).compare("--n-parts") == 0)
        {
        	nParts = atoi(argv[i + 1]);
	    }
        else if (string(argv[i]).compare("--prefix") == 0)
        {
        	sOutputPrefix = string(argv[i + 1]);
        }
        else if (string(argv[i]).compare("--gz") == 0)
        {
        	bIsGzip = true;
        }
        else if (string(argv[i]).compare("--version")== 0)
        {
        	bShowVersion = true;
	    }
    }

    sFileInputName = string(argv[argc-1]);

    if(bShowHelp)
    {
    	showHelp();
    	return 0;
    }

    if(bShowVersion)
    {
    	showVersion();
    	return 0;
    }

    if(nParts == 0)
    {
    	cout << "Sorry!! No number of parts to divide the file was specified!!" << endl;
    	return 1;
    }

    if(sFileInputName.empty())
    {
    	cout << "Sorry!! No input file name specified!!" << endl;
    	return 1;
    }

    return splitFile(sFileInputName,sOutputPrefix,nParts,bIsGzip);

    return 0;
}

