/*
 * main.cpp
 *
 *  Created on: 1 Abr, 2015
 *      Author: marcos
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/**
 * \brief makeBedIntervalsAssembly application help
 */
void printHelp()
{
	fprintf (stdout, "makeBedIntervalsAssembly Creates a bed of chromosome intervals \n\n");
	fprintf (stdout, "\t-chrName <chrName>     : Chromosome name \n");
	fprintf (stdout, "\t-chrLen  <chrLen>      : Chromosome length.\n");
	fprintf (stdout, "\t-outFile  <outFile>    : Path to store bed file .\n\n");
}

/**
 * \brief main function
 */
int main(int argc, char **argv)
{
	char chrName [512];
	unsigned int chrLen = 0;
	char outFile [512];

    	//1. PARSE ARGUMENTS
	for (int i=1; i<argc; i++)
	{
		if (!strcmp(argv[i], "-chrName"))
		{
			strcpy(chrName, argv[i+1]);
		}
	    else if (!strcmp(argv[i], "-chrLen"))
		{
			chrLen = atoi(argv[i+1]);
		}
	    else if (!strcmp(argv[i], "-outFile"))
		{
	    	strcpy(outFile, argv[i+1]);
		}
	    else if (!strcmp(argv[i], "-h"))
	    {
	    	printHelp();
	    	return 0;
	    }
	}

	//2. size of kmers and step
	unsigned int kmerLen = 36;
    	unsigned int step = 5;

	//3.Open Bed output file
	FILE * fout;
	fout = fopen(outFile,"w");

	unsigned int begin = 1;
	unsigned int end = begin + kmerLen -1;

	while (end <= chrLen)
	{
		char currentLine [512];
		char sBegin [512];
		char sEnd [512];
		//chr name
		strcpy(currentLine, chrName);
		strcat(currentLine, "\t");
		//begin
		sprintf(sBegin, "%d", (begin -1));
		strcat(currentLine, sBegin);
		//end
		strcat(currentLine, "\t");
		sprintf(sEnd, "%d", (end));
		strcat(currentLine, "\n");

		//Print to File
		fwrite (currentLine , 1, strlen(currentLine), fout);

		begin += step;
		end = begin + kmerLen - 1;
	}

	fclose(fout);
}


