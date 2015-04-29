#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

int main(int argc, char *argv[]) {
	// Define lists
	unsigned long long int mono[5] = {0};
	unsigned long long int di[5][5] = {{0}};
	unsigned long long int din[5][5] = {{0}};
	unsigned short int cur = 5, last = 5, thrid = 5;
	unsigned long long int sequences = 1;

	if (argc == 1) {
		cerr << "galculator 1.0" << endl;
		cerr << "Counts mono, di and xNy derivation of nucleotide sequences in a fasta-file" << endl;
		cerr << "Lower and upper case letters are treated equally" << endl << endl;
		cerr << "Usage: galculator FASTA-FILE " << endl;
		cerr << "Output: #sequences #nucleotides | NUC # ... NUCaNUCb # ... NUCanNUCb # ... "<< endl;
		cerr << endl << "Questions, comments, bug reports to lechner@staff.uni-marburg.de" << endl;
		exit(EXIT_FAILURE);
	}

	// Open file
	ifstream is;
	filebuf * fasta;
	fasta = is.rdbuf();
	fasta->open (argv[1],ios::in);

	if (fasta->is_open()) {
		char letter;
		unsigned long long int pos = 0;	// Position in current entry
		bool scan = false;
		while (fasta->sgetc()!=EOF) {
			letter = fasta->sbumpc();
			pos++;
			if (scan) {
			// decide which letter
			switch(letter) {
				case 'a' :
				case 'A' :  cur = 0; break;
				case 'c' :
				case 'C' :  cur = 1; break;
				case 'g' :
				case 'G' :  cur = 2; break;
				case 't' :
				case 'T' :
				case 'u' :
				case 'U' :  cur = 3; break;
				case 'n' :
				case 'N' :  cur = 4; break;
				default :   cur = 5;
			}

			if (cur == 5) {
				if (letter == '>') {
					sequences++;
					scan = false;
					continue;
				}
				if (letter == '\n') {continue;}
				cerr << "Warning: Char " << letter << " in sequence " << sequences << " " << argv[1] << " treated as N." << endl;
				cur = 4;
			}

			// Store
			mono[cur]++;
			if (pos >= 2) di[last][cur]++;
			if (pos >= 3) din[thrid][cur]++;
			// shift
			thrid = last;
			last = cur;
			}
			else {
				if (letter == '\n') {scan = true;}
			}

		}
		fasta->close();

		// Alpha assignment
		//				 0   1   2   3   4
		char alpha[] = {'A','C','G','T','N'};

		// Count-Loop finished
		unsigned long long int alln = mono[0]+mono[1]+mono[2]+mono[3]+mono[4];

		printf("%lld %lld | ",sequences,alln);

		// Mono
		for (unsigned short int i = 0; i < 5; i++) {
				printf("%c %lld ",alpha[i],mono[i]);
		}
		// Di
		for (unsigned short int i = 0; i < 4; i++) {
			for (unsigned short int j = 0; j < 4; j++) {
				printf("%c%c %lld ",alpha[i],alpha[j],di[i][j]);
			}
		}
		// DiN
		for (unsigned short int i = 0; i < 4; i++) {
					for (unsigned short int j = 0; j < 4; j++) {
						printf("%cn%c %lld ",alpha[i],alpha[j],din[i][j]);
					}
		}
		puts("");

		} else {
			cerr << "Error: Unable to open file " << argv[1] << endl;
			exit(EXIT_FAILURE);
		}

		return EXIT_SUCCESS;
}
