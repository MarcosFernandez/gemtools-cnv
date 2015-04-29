#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <math.h>

using namespace std;

struct appParameters
{
	string iTable; // the table for 1st portion with chrom, chromStart, chromEnd
	string jTable; // the table for 2nd portion with chrom, chromStart, chromEnd
	bool fOption; //-f option for -i comes from file with chrom, chromStart, chromEnd as first 3 columns
	bool tOption; //-t option for -j comes from file with chrom, chromStart, chromEnd as first 3 columns
	bool left;
	bool right;
	string iCond;
	string jCond;
	string outputFile;

	appParameters()
	{
    		iTable = "";
    		jTable = "";
    		fOption = false;
    		tOption = false;
    		left = false;
    		right = false;
    		iCond = "";
    		jCond = "";
    		outputFile = "";
	}
};

struct coord
{
	long int start;
	long int end;

	//Constructor
	coord(long int ini,long int fin) : start(ini), end(fin) { }


	bool operator<( const coord& other ) const
	{
	    if (start == other.start)
	    {
	    	return end < other.end;
	    }
	    else
	    {
	    	return start < other.start;
	    }
	}
};

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
 * \brief Load File in Memory structure
 * \param nameFile - const string & file path
 * \param mapChrLoc - reference to map structure chr and locations(st end)
 * \return false if it was not possible to complete the operation
 */
bool loadFileToMap(const string & nameFile,map <string, vector <coord> > & mapChrLoc)
{
	ifstream inFile;
	inFile.open(nameFile.c_str());
	string lineRead = "";

	bool success = true;

	if (inFile.is_open())
	{
		if (!inFile.good())
		{
			cout << "Sorry!! Was not possible to open file " << nameFile << endl;
			success = false;
		}

		while (getline (inFile,lineRead))
		{
			if(lineRead.empty()) continue;

			vector <string> fields = split(lineRead,'\t');

			coord region = coord(atoi(fields[1].c_str()),atoi(fields[2].c_str()));
			//if not exists chromosome in map
			if ( mapChrLoc.find(fields[0]) == mapChrLoc.end() )
			{
				vector < coord > chrLocations;
				chrLocations.push_back(region);
				mapChrLoc [fields[0]] = chrLocations;
			}
			else
			{
				//Chromosome already exists in map, so append coordinates
				mapChrLoc [fields[0]].push_back(region);
			}
		}
	}
	else
	{
		cout << "Sorry it was not possible to open file: " <<  nameFile << endl;
		success = false;
	}

	inFile.close();

	return success;
}

/**
 * \brief Get all chromosomes from a map structure
 * \param mapChrLocation - map structure chromosome locations
 * \param setChroms - set collection of chromosome names
 */
void getAllChroms(const map <string, vector<coord> > & mapChrLocation,set <string> & setChroms)
{
	for (std::map <string, vector < coord > >::const_iterator it=mapChrLocation.begin(); it!=mapChrLocation.end(); ++it)
	{
		setChroms.insert(it->first);
	}
}

/**
 * \brief Print to file descriptor out file chromosome locations
 * \param vLocations - vector of chromosome coordinates
 * \param nameChrom - Chromosome name
 * \param outFile - file descriptor managed by ofsrteam reference
 */
void printLocations(const vector <coord> & vLocations,const string & nameChrom,ofstream & outFile)
{
	for (std::vector < coord >::const_iterator it=vLocations.begin(); it!=vLocations.end(); ++it)
	{
		outFile << nameChrom << '\t' << it->start << '\t' << it->end << std::endl;
	}


}

/**
 * \brief Prints application basic arguments to stdout
 */
void printHelp()
{
	std::cout << "twoPartOvpMgsrt " << std::endl;
    std::cout << "Outputs the overlap/exclude of two tables using merge sort algorithm" << std::endl;
    std::cout << "" << std::endl;

    std::cout << " -i the table for 1st portion with chrom, chromStart, chromEnd" << std::endl;
    std::cout << " -j the table for 2nd portion with chrom, chromStart, chromEnd" << std::endl;
    std::cout << " -f option for -i comes from file with chrom, chromStart, chromEnd as first 3 columns" << std::endl;
    std::cout << " -t option for -j comes from file with chrom, chromStart, chromEnd as first 3 columns" << std::endl;
    std::cout << " -c condition for table i. Don't set it if no." << std::endl;
    std::cout << "    Given as tableName.column=value. Concate with ' AND ' if multiple" << std::endl;
    std::cout << " -d condition for table j. Don't set it if no" << std::endl;
    std::cout << " -L output the exclude from left (-i) table (part of -i table not overlapped with -j table)" << std::endl;
    std::cout << " -R output the exclude from right(-j) table (part of -j table not overlapped with -i table)" << std::endl;
    std::cout << "    If -L, -R are set simultaneously, only -L is used" << std::endl;
    std::cout << " -o output file" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "" << std::endl;

    std::cout << "Eg:" << std::endl;
    std::cout << "       twoPartOvpMgsrt         -i 1.tab -f             -j  2.tab -t          -o /tmp/both.tab" << std::endl;
    std::cout << "       twoPartOvpMgsrt         -i 1.tab -f             -j  2.tab -t     -L   -o /tmp/onlyIn1.tab" << std::endl;
    std::cout << "       twoPartOvpMgsrt         -i 1.tab -f             -j  2.tab -t     -R   -o /tmp/onlyIn2.tab" << std::endl;

    std::cout << "" << std::endl;

    std::cout << "NOTE: There can be no overlap within tables itself, ie, the tables have to have N(on)R(edundant) bases for it to work" << std::endl;

}

/**
 * \brief Control application parameters
 * \param parameters  appParameters structure which contains app parameters
 * \returns false if any mandatory parameter is not properly configured
 */
bool parametersControl(appParameters & parameters)
{
	bool test = true;

	if (parameters.iTable.empty())
	{
		std::cout << "Sorry!! -i Was not defined." << std::endl;
		test = false;
	}

	if (parameters.jTable.empty())
	{
		std::cout << "Sorry!! -j Was not defined." << std::endl;
		test = false;
	}

	if (parameters.outputFile.empty())
	{
		std::cout << "Sorry!! -o Was not defined." << std::endl;
		test = false;
	}

	if (parameters.left && parameters.right)
	{
		parameters.right = false;
	}

	return test;
}


/**
 * \brief Process arguments
 * \param vArguments vector of arguments
 * \param parameters  appParameters structure which contains app parameters
 */
bool processArguments(const vector <string> & vArguments, appParameters & parameters)
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
					if (sArgument.compare("-i") == 0)
					{
						parameters.iTable = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-j") == 0)
					{
						parameters.jTable = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-f") == 0)
					{
						parameters.fOption = true;
					}
					else if (sArgument.compare("-t") == 0)
					{
						parameters.tOption =  true;
					}
					else if (sArgument.compare("-c") == 0)
					{
						parameters.iCond = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-d") == 0)
					{
						parameters.jCond = vArguments.at(iParameterPosition+1);
					}
					else if (sArgument.compare("-L") == 0)
					{
						parameters.left = true;
					}
					else if (sArgument.compare("-R") == 0)
					{
						parameters.right = true;
					}
					else if (sArgument.compare("-o") == 0)
					{
						parameters.outputFile = vArguments.at(iParameterPosition+1);
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

	return parametersControl(parameters);
}


/**
 * \brief main function
 */
int main(int argc, char *argv[])
{
	//1. Global parameters
	appParameters parameters;

	//2. Process arguments
	vector <string> vArguments;

	for (int i=0; i < argc; i++)
	{
		vArguments.push_back(string(argv[i]));
	}

	if (!processArguments(vArguments,parameters))
	{
		return -1;
	}

	//3. Load iTable File
	map <string, vector<coord> > firstList;
	if(!loadFileToMap(parameters.iTable,firstList))
	{
		return -1;
	}

	//4. Load jTables
	map <string, vector<coord> > secondList;
	if(!loadFileToMap(parameters.jTable,secondList))
	{
		return -1;
	}

	//5. Get All Chromosomes
	set <string> allChroms;

	getAllChroms(firstList,allChroms);
	getAllChroms(secondList,allChroms);

	//6 Open output file and star merging process
	ofstream outFile;
	outFile.open(parameters.outputFile.c_str());

	for (set <string>::iterator chrom = allChroms.begin(); chrom != allChroms.end(); ++chrom)
	{
        //6.1 Sort Chromosome locations for table i (first) and j (second)
		if( firstList.find(*chrom) != firstList.end() ) //If exists chrom in map (i table)
		{
			sort(firstList[*chrom].begin(),firstList[*chrom].end());
		}

		if( secondList.find(*chrom) != secondList.end() )//If exists chrom in map (j table)
		{
			sort(secondList[*chrom].begin(),secondList[*chrom].end());
		}

		//6.A Run option what is in table i (first) but not in table j (second)
		if(parameters.left)
		{
			if( firstList.find(*chrom) == firstList.end() ) //If not exists chrom in map (i table)
			{
				continue;
			}

			if (secondList.find(*chrom) == secondList.end())//If not exists chrom in map (j table)
			{
				printLocations(firstList[*chrom],*chrom,outFile);
				continue;
			}
		}

		//6.B Run option what is in table j (second) but not in table i (first)
		if(parameters.right)
		{
			if( secondList.find(*chrom) == secondList.end() ) //If not exists chrom in map (j table)
			{
				continue;
			}

			if (firstList.find(*chrom) == firstList.end())//If not exists chrom in map (i table)
			{
				printLocations(secondList[*chrom],*chrom,outFile);
				continue;
			}
		}

		//6.C Run option what is in table i (first) and in table j (second)
		if(!parameters.left && !parameters.right)
		{
			//If is not present chromose in some table get next chromosome
			if(firstList.find(*chrom) == firstList.end() || secondList.find(*chrom) == secondList.end())
			{
				break;
			}
		}

		//Chromosome is defined for first (i) and second (j) table
		long int i = 0;
		long int j = 0;
		pair <long int,long int> leftRightStart = make_pair(-1,-1);

		while ( (i < firstList[*chrom].size()) && (j < secondList[*chrom].size()) )
		{
			coord firstRoA = firstList[*chrom][i];
			coord secondRoA = secondList[*chrom][j];

			if(firstRoA.end <= secondRoA.start)
			{
				if(parameters.left && leftRightStart.first == -1 && firstRoA.end > firstRoA.start )
				{
					outFile << *chrom << '\t' << firstRoA.start << '\t' << firstRoA.end << std::endl;
				}

				if(parameters.left && leftRightStart.first != -1 && firstRoA.end > leftRightStart.first )
				{
					outFile << *chrom << '\t' << leftRightStart.first << '\t' << firstRoA.end << std::endl;
				}

				leftRightStart.first = -1;
				i++;
				continue;
			}

			if(secondRoA.end <= firstRoA.start)
			{
				if(parameters.right && leftRightStart.second == -1 && secondRoA.end > secondRoA.start )
				{
					outFile << *chrom << '\t' << secondRoA.start << '\t' << secondRoA.end << std::endl;
				}

				if(parameters.right && leftRightStart.second != -1 && secondRoA.end > leftRightStart.second )
				{
					outFile << *chrom << '\t' << leftRightStart.second << '\t' << secondRoA.end << std::endl;
				}

				leftRightStart.second = -1;
				j++;
				continue;
			}

            	    // Pick the greatest start and least end
		    long int start = (firstRoA.start > secondRoA.start) ? firstRoA.start : secondRoA.start;
		    long int end   = (firstRoA.end > secondRoA.end) ? secondRoA.end : firstRoA.end;

		    if( !parameters.left && !parameters.right && end > start )
		    {
		       outFile << *chrom << '\t' << start << '\t' << end << std::endl;
		    }

		    if( parameters.left && firstRoA.start < start )
		    {
		       if( leftRightStart.first == -1 && start > firstRoA.start )
		       {
		    	   outFile << *chrom << '\t' << firstRoA.start << '\t' << start << std::endl;
		       }

		       if( leftRightStart.first != -1 && start > leftRightStart.first )
		       {
		    	   outFile << *chrom << '\t' << leftRightStart.first << '\t' << start << std::endl;
		       }
	  	     }

		     if( parameters.right && secondRoA.start < start )
			 {
		       if( leftRightStart.second == -1 && start > secondRoA.start )
			   {
				   outFile << *chrom << '\t' << secondRoA.start << '\t' << start << std::endl;
			   }

			   if( leftRightStart.second != -1 && start > leftRightStart.second )
			   {
				   outFile << *chrom << '\t' << leftRightStart.second << '\t' << start << std::endl;
			   }
			 }

     		 // move the one whose end is smaller
		     if( firstRoA.end == end ) i++;
		     if( secondRoA.end == end ) j++;

		     leftRightStart.first = (firstRoA.end > end) ? end : -1;
		     leftRightStart.second = (secondRoA.end > end) ? end : -1;

		}

        // Append those that are left from -i table
	    while( parameters.left && i < firstList[*chrom].size() )
	    {
		    if( leftRightStart.first != -1 )
		    {
		      	if( firstList[*chrom][i].end > leftRightStart.first )
		       	{
		       	    outFile << *chrom << '\t' << leftRightStart.first << '\t' << firstList[*chrom][i].end << std::endl;
		       	}

		       leftRightStart.first = -1;
		    }
		    else
		    {
		       	outFile << *chrom << '\t' << firstList[*chrom][i].start << '\t' << firstList[*chrom][i].end << std::endl;
		    }
		    i++;
	    }

		// Append those that are left from -j table
	    while( parameters.right && j < secondList[*chrom].size() )
	    {
		    if( leftRightStart.second != -1 )
		    {
		        if( secondList[*chrom][j].end > leftRightStart.second )
		        {
		        	outFile << *chrom << '\t' << leftRightStart.second << '\t' << secondList[*chrom][j].end << std::endl;
		        }

		        leftRightStart.second = -1;
		    }
		    else
		    {
		       outFile << *chrom << '\t' << secondList[*chrom][j].start << '\t' << secondList[*chrom][j].end << std::endl;
		    }
		    j++;
	    }


	}

	outFile.close();
    return 0;
}
