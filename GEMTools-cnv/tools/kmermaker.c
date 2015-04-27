#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kmermaker.h"


/* This program splits the reads into a given length*/
int main(int argc, char **argv){
  FILE *fasta;
  int klen=0;
  int i = 0;
  char name[MAXKMER];
  char fname[100];
  char line[100];

  

  fname[0]=0;
  /*Checking arguments*/
  if((argc < 6 ) || (!strcmp(argv[1], "-h"))) 
  {
     printf("kmermaker -i nameFile -k lengthKmers -s window -f firstPosition -p hasHeader \n");
     return 0;
  }


  /*Reading arguments */
  for (i=0; i<argc; i++)
  {
    if (!strcmp(argv[i], "-i"))
      strcpy(fname, argv[i+1]);
    else if (!strcmp(argv[i], "-k"))
      klen = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-s"))
      slide = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-f"))
      first = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-pass"))
      pass_first = 1;
  }

  if (!strcmp(fname, "stdin"))
    fasta  = stdin;
  else
    fasta = fopen(fname, "r");
  
  i=0; seq[0]=0; name[0]=0; qual[0]=0;
  
  /*inicialization string of N nucleotids according to klen (size reads)*/
  for (i=0;i<klen/3;i++)
  {
    nmer[i]='N';
  }
  nmer[klen/3]=0;

  /*Parsing the file*/
  while (fscanf(fasta, "@%s\n%s\n", name, seq) > 0)
  {
    fgets(line, 100, fasta);
    fscanf(fasta, "%s\n", qual);


    breakthis(name, seq, qual, klen, first -1);
   
  }
  
  return fclose(fasta);
}

void  breakthis(char *name, char *seq, char *qual, int klen, int start){
  int i;
  int len = strlen(seq);
  char kmer[MAXKMER];
  char qmer[MAXKMER];

  /*Copies the sequence and qualities for each window*/
  for (i=start; i<=len-klen;i+=slide)
  {
    memcpy(kmer, seq+i, klen);
    memcpy(qmer, qual+i, klen);
    kmer[klen]=0; qmer[klen]=0;

    if (i==start && pass_first) continue;

    printf("@%s\n%s\n+\n%s\n", name, kmer, qmer);
    
  }
}

