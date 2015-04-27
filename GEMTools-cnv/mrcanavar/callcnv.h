#ifndef __CALLCNV
#define __CALLCNV

#include <stdio.h>
#include "globals.h"
#include "utils.h"
#include "gcnorm.h"

void call_cnv(char *, char *);
void readDepth(char *);
void  print_copy_numbers(char *);
void dump_text_windows(char *, enum WINDOWTYPE);
void conc_depth(char **, int, char *);
#endif
