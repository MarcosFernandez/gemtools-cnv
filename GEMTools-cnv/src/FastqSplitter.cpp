/*
 * main.cpp
 *
 *  Created on: 3 Dec, 2014
 *      Author: marcos
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <cstdlib>
#include <ctime>
#include <cstdio>

using namespace std;

class Fastq
{
	string readHeader;
	string readSequence;
	string readPlus;
	string readQualities;
	size_t charSize;
	size_t carryReturnSize;
	bool printed;


    public:
	Fastq ()
        {
	    	readHeader = "";
	    	readSequence = "";
	    	readPlus = "";
	    	readQualities = "";
	    	charSize = sizeof(char);
	    	carryReturnSize = sizeof("\n");
	    	printed = false;
        }

	string getHeader(){ return readHeader;};
	string getSequence(){return readSequence; };
	string getPlus(){return readPlus; };
	string getQualities(){return readQualities; };

        unsigned int headerLength (){ return readHeader.length();};
        unsigned int sequenceLength (){ return readSequence.length();};
        unsigned int plusLength (){ return readPlus.length();};
        unsigned int qualitiesLength (){ return readQualities.length();};

        bool headerEmpty (){ return readHeader.empty();};
        bool sequenceEmpty (){ return readSequence.empty();};
        bool plusEmpty (){ return readPlus.empty();};
        bool qualitiesEmpty (){ return readQualities.empty();};

        bool fastqEmpty()
        {
        	return readHeader.empty();
        }

        bool checkHeader()
        {
        	 if (readHeader.at(0) != '@')
        	 {
        		 return false;
        	 }
        	 return true;
        }

        bool checkSequenceQualitiesLength()
        {
        	if(this->sequenceLength() != this->qualitiesLength())
        	{
        		return false;
        	}
        	return true;
        }

        void readFastq(ifstream &ifs)
        {
        	getline (ifs,readHeader);
        	getline (ifs,readSequence);
        	getline (ifs,readPlus);
        	getline (ifs,readQualities);
        	printed = false;
        }

        void saveFastq(FILE * fout)
        {
            if(!printed)
        	{
            	fprintf(fout, "%s\n", readHeader.c_str());
            	fprintf(fout, "%s\n", readSequence.c_str());
            	fprintf(fout, "%s\n", readPlus.c_str());
            	fprintf(fout, "%s\n", readQualities.c_str());
            	printed = true;
        	}
        }
};

static std::string int2String(int number){ std::stringstream ss; ss << number; return ss.str();}

/**
 * \brief Print Help File
 */
void showVersion()
{
    cout << "FASTQ Splitter C++ 0.1" << endl;
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
}


/**
 * \brief GetFiles size
 * \param inFileName File Name
 * \param totalNseq Total number of sequences
 */
void getFileSize (const string & inFileName,unsigned int & totalNseq)
{
	ifstream infile;
	infile.open(inFileName.c_str());

	Fastq fastq;

	unsigned int lineCount = 0;

	if (infile.is_open())
	{
			if (infile.good())
			{
				fastq.readFastq(infile);
			}

			while (!fastq.fastqEmpty())
			{
				//1. READ HEADER
				lineCount ++;
				if (!fastq.checkHeader())
                {
                	cout << "Error parsing line " << lineCount << " FASTQ entry does not start with @" << endl;
                	return;
                }
				//2.READ SEQUENCE
                lineCount ++;
                //3.READ PLUS
                lineCount ++;
				//4.READ QUALITIES
                lineCount ++;
                if(!fastq.checkSequenceQualitiesLength())
                {
                	cerr << "Warning: Misformatted FASTQ entry in input line " << lineCount << ": quality length (" << fastq.qualitiesLength();
                	cerr << ") differs from sequence length (" << fastq.sequenceLength() << ")" << endl;
                }

                totalNseq ++;

				fastq.readFastq(infile);
			}
	}

	infile.close();
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
 * \param totalNseq Total Number of Sequences
 */
void splitFile(const string & inFileName,unsigned int nParts,unsigned int & totalNseq)
{
	FILE * fout;
	Fastq fastq;
	ifstream infile;
	infile.open(inFileName.c_str());

	cout << inFileName.c_str() << endl;

	string base = inFileName;

	totalNseq = 0;

	getFileSize(inFileName,totalNseq);
	cout << ": " << totalNseq << " Reads " << endl;

	startingMessage(nParts);
	unsigned int written = 0;

    if (infile.is_open())
	{
			if (infile.good())
			{
				for (unsigned int part=1; part <= nParts; part++)
				{
					unsigned int should_write =  round((float)((float) totalNseq * (float) part / (float) nParts));
					string part_file = base + ".part-" + int2String(part);

					fout=fopen(part_file.c_str(), "w");
					while (written < should_write)
					{
						fastq.readFastq(infile);
						written += 1;
						if(!fastq.fastqEmpty()) fastq.saveFastq(fout);
					}

					fclose(fout);
				}
			}
	}

    infile.close();
}


/**
 * \brief main function
 */
int main(int argc, char *argv[])
{
	bool bShowHelp = false;
	bool bShowVersion = false;
	unsigned int nParts = 0;
	string sFileInputName = "";

	unsigned int totalNseq = 0;

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

    splitFile(sFileInputName,nParts,totalNseq);

    return 0;
}
