#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <sstream>
#include <algorithm>
#include <math.h>

using namespace std;

static std::vector<std::string> split(std::string s, char delim, std::vector<std::string> elems)
{
	std::stringstream ss(s);
	std::string item;
	while(std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

static std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	return split(s, delim, elems);
}

/**
 * \brief Prints application basic arguments to stdout
 */
void printHelp()
{
	std::cout << "USAGE MergeCoord " << std::endl;
	std::cout << "-f1 input file 1 with coordinates chrom start end" << std::endl;
	std::cout << "-h1 there is a header in -f1" << std::endl;
	std::cout << "optional, by default it is not" << std::endl;
	std::cout << "-f2 optional, input file 2 with coordinates chrom start end" << std::endl;
	std::cout << "-h2 there is a header in -f2" << std::endl;
	std::cout << "optional, by default it is not" << std::endl;
	std::cout << "-d  optional, by default 0, maximum distance to merge blocks if they are separated by a length smaller than d" << std::endl;
	std::cout << "-outF output file \n" << std::endl;

}


/**
 * \brief Checks parameters quality and coherence
 * \return boolean value which indicates the arguments evaluation
 */
bool parametersControl(string & sInputFile1,string & sInputFile2,string & sOutputFile)
{
	bool test = true;

    if(sInputFile1.empty())
    {
    	std::cout << "Sorry!! -f1 Was not defined." << std::endl;
    	test = false;
    }

    if(sOutputFile.empty())
    {
    	std::cout << "Sorry!! -outF Was not defined." << std::endl;
	    test = false;
    }

    return test;
}


/**
 * \brief Process arguments
 * \param vArguments vector of arguments
 * \param sInputFile1  string Input file 1
 * \param isHeaderInput1 bool File One has header
 * \param sInputFile2  string Input file 2
 * \param isHeaderInput2 bool File two has header
 * \param distance int distance
 * \param sOutputFile string Output name file
 */
bool processArguments(const vector <string> & vArguments,string & sInputFile1,bool & isHeaderInput1,string & sInputFile2,bool & isHeaderInput2,unsigned long int & distance,string & sOutputFile)
{
	if(vArguments.size() <= 1 || vArguments.at(1) == "h")
	{
		printHelp();
		return false;
	}

	int iParameterPosition = 0;

	//ARGUMENTS MEANING ASIGNATION
	for (vector <string>::const_iterator it = vArguments.begin(); it != vArguments.end(); it++)
	{
			string sArgument = (*it);

			if (sArgument.at(0) == '-')
			{
					if (sArgument.compare("-f1") == 0)
					{
						sInputFile1 = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-h1") == 0)
					{
						isHeaderInput1 = true;
					}
					else if (sArgument.compare("-f2") == 0)
					{
						sInputFile2 = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-h2") == 0)
					{
						isHeaderInput2 = true;
					}
					else if (sArgument.compare("-d") == 0)
					{
						distance = (unsigned) atoi(vArguments.at(iParameterPosition+1).c_str());
					}
					else if (sArgument.compare("-outF") == 0)
					{
						sOutputFile = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-h") == 0 || sArgument.compare("-help") == 0 )
					{
						printHelp();
						return false;
					}
					else
					{
						std::cout << sArgument  << " is not a valid argument. Use -h to get more help!!" << std::endl;
						return false;
					}

			}
			iParameterPosition ++;
	}

	return parametersControl(sInputFile1,sInputFile2,sOutputFile);

}

/**
 * \brief Load corrdinates from file to map structure
 * \param sInputFile - string Input Name File
 * \param coordinate - map structure
 */
void loadCoordinates(const string & sInputFile,map <string,vector < pair< long int, long int> > > & coordinates, bool bHasHeader)
{
	ifstream inFile;
	inFile.open(sInputFile.c_str());
	string lineRead = "";

    if (inFile.is_open())
    {
        if (!inFile.good())
		{
		    cout << "Sorry!! Was not possible to open file " << sInputFile << endl;
        }

        if(bHasHeader)
        {
        	getline (inFile,lineRead);
        }

        while (getline (inFile,lineRead))
        {
        	if(lineRead.empty()) continue;

			vector <string> fields = split(lineRead,'\t');

			pair < long int, long int> stEnd ( atoi(fields[1].c_str()),  atoi(fields[2].c_str()));

			if (stEnd.second < stEnd.first)
			{
				long int tmp = stEnd.first;
				stEnd.first = stEnd.second;
				stEnd.second = tmp;
			}

			//if not exists chromosome in map
			if ( coordinates.find(fields[0]) == coordinates.end() )
			{
				vector < pair <long int,long int> > chrLocations;
				chrLocations.push_back(stEnd);
				coordinates [fields[0]] = chrLocations;
			}
			else
			{
			    //Chromosome already exists in map, so append coordinates
				coordinates [fields[0]].push_back(stEnd);
			}
        }
    }
    else
    {
    	cout << "Sorry!! Was not possible to open file " << sInputFile << endl;
    }

    inFile.close();
}


/**
 * \brief main function
 */
int main(int argc, char *argv[])
{
	//1. Global parameters
	string sInputFile1 = "";
	bool isHeaderInput1 = false;
	string sInputFile2 = "";
	bool isHeaderInput2 = false;
	unsigned long int distance = 0;
	string sOutputFile = "";

	//2. Process arguments
	vector <string> vArguments;

	for (int i=0; i < argc; i++)
	{
			vArguments.push_back(string(argv[i]));
	}

	if (!processArguments(vArguments,sInputFile1,isHeaderInput1,sInputFile2,isHeaderInput2,distance,sOutputFile))
	{
		return -1;
	}


	//3. Process Input Files
	map <string,vector < pair  < long int, long int> > > coordinates;
	loadCoordinates(sInputFile1,coordinates,isHeaderInput1);
	if (!sInputFile2.empty())
	{
		loadCoordinates(sInputFile2,coordinates,isHeaderInput2);
	}

	//4. Merge Coordinates a store data in output file
	ofstream outFile;
	outFile.open(sOutputFile.c_str());

    for (std::map <string,vector < pair  <long int, long int> > >::iterator it=coordinates.begin(); it!=coordinates.end(); ++it)
    {
    	vector < pair  <long int,long int> > locations = it->second;
    	std::sort(locations.begin(), locations.end());

    	//Merging chromosome locations
    	for (unsigned int i = 0; i < locations.size()-1; i++)
    	{
    		if (locations[i+1].first < (1 + locations[i].second + distance) )
    		{
    			locations[i+1].first = locations[i].first;

    			if(locations[i+1].second < locations[i].second)
    			{
    				locations[i+1].second = locations[i].second;
    			}

    			locations[i].first = -1;
    			locations[i].second = -1;
    		}
    	}

    	//Print to Out file Chromosome Locations
    	for (vector < pair  <long int,long int> > ::iterator loc = locations.begin(); loc!=locations.end(); ++loc)
    	{
    		if (loc->first != -1 && loc->second != -1)
    		{
    			outFile << it->first << '\t' << loc->first << '\t' << loc->second << endl;

    		}
    	}
    }

    outFile.close();
}
