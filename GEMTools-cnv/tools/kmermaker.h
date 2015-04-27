#define MAXKMER 100
#define MAXLEN  10000


void  breakthis(char *name, char *seq, char *, int klen, int start);

int slide=1;
int first=1;
char seq[MAXLEN];
char qual[MAXLEN];
char nmer[MAXKMER];
int pass_first = 0;
