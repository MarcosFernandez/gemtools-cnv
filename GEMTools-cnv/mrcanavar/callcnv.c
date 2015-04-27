#include "callcnv.h"

void call_cnv(char *depthFile, char *out_prefix){

  int i, j;

  float *gclookup;
  float *gclookup_x;

  char *fname;
  char logname[2 * MAX_STR];
  FILE *log;

  if (GENOME_CONF == NULL)
    print_error("Select genome configuration file (input) through the -conf parameter.\n");
  if (depthFile == NULL)
    print_error("Select read depth output file through the -depth parameter.\n");
  if (out_prefix == NULL)
    print_error("Select output file prefix through the -o parameter.\n");



  loadRefConfig(GENOME_CONF);

  readDepth(depthFile);


  sprintf(logname, "%s.log", out_prefix);
  log = my_fopen(logname, "w", 0);

  gclookup   = (float *) getMem(sizeof(float) * GC_BIN);
  gclookup_x = (float *) getMem(sizeof(float) * GC_BIN);

  
  /* trivial control removal */
  
  fprintf(stdout, "Control region cleanup...");
  fflush(stdout);
    
  /* add stdev calculation here */
  
  
  for (i=0; i<num_chrom; i++){
    for (j=0; j<chromosomes[i]->lw_cnt; j++){
      if (chromosomes[i]->lw[j].depth > (float) LW_MEAN * 2.0 || chromosomes[i]->lw[j].depth < (float) LW_MEAN / 10.0)
	chromosomes[i]->lw[j].isControl = 0;
      if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
	chromosomes[i]->lw[j].isControl = 0;
    }
    for (j=0; j<chromosomes[i]->sw_cnt; j++){
      if (chromosomes[i]->sw[j].depth > (float) SW_MEAN * 2.0 || chromosomes[i]->sw[j].depth < (float) SW_MEAN / 10.0)
	chromosomes[i]->sw[j].isControl = 0;
      if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
	chromosomes[i]->sw[j].isControl = 0;
    }
    
    for (j=0; j<chromosomes[i]->cw_cnt; j++){
      if (chromosomes[i]->cw[j].depth > (float) CW_MEAN * 2.0 || chromosomes[i]->cw[j].depth < (float) CW_MEAN / 10.0)
	chromosomes[i]->cw[j].isControl = 0;
      if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
	chromosomes[i]->cw[j].isControl = 0;
      }
  }
    
  
  fprintf(stdout, "\n");
  
  
  norm_until_converges(CW, gclookup, gclookup_x);
  
  norm_until_converges(LW, gclookup, gclookup_x);
  
  norm_until_converges(SW, gclookup, gclookup_x);
  //}
  
  fprintf (stdout, "Writing normalized CW depth to: %s.cw_norm.bed.\n", out_prefix);
  
  fname = (char *) getMem(sizeof (char) * (strlen(out_prefix) + strlen(".copynumber.bed") + 2));
  
  sprintf (fname, "%s.cw_norm.bed", out_prefix);
  dump_text_windows(fname, CW);
  
  fprintf (stdout, "Writing normalized LW depth to: %s.lw_norm.bed.\n", out_prefix);
  
  sprintf (fname, "%s.lw_norm.bed", out_prefix);
  dump_text_windows(fname, LW);
  
  fprintf (stdout, "Writing normalized SW depth to: %s.sw_norm.bed.\n", out_prefix);
  
  sprintf (fname, "%s.sw_norm.bed", out_prefix);
  dump_text_windows(fname, SW);
  
  sprintf (fname, "%s.copynumber.bed", out_prefix);
  fprintf (stdout, "Writing copy numbers to: %s.copynumber.bed. \n", out_prefix);
  print_copy_numbers(fname);
  
  freeMem(fname, (strlen(out_prefix) + strlen(".copynumber.bed") + 2));

    
  fprintf (log, "\nmrCaNaVaR version %s\nLast update: %s\n\n", VERSION, LAST_UPDATE);
  fprintf (log, "Calculating library %s\n", out_prefix);
  fprintf (log, "GC correction mode: %s\n", MULTGC == 1 ? "MULTIPLICATIVE" : "ADDITIVE");

  fprintf (log, "\nAfter GC Correction:\n--------------------\n");
  fprintf (log, "Sample Gender: %s.\n", GENDER == MALE ? "Male" : "Female");
  fprintf (log, "CW Average Read Depth: %f, Standard Deviation: %f\n", CW_MEAN, CW_STD);


  if (GENDER == MALE)
    fprintf (log, "CW Average chrX Read Depth: %f, Standard Deviation: %f\n", CW_MEAN_X, CW_STD_X);

  fprintf (log, "LW Average Read Depth: %f, Standard Deviation: %f\n", LW_MEAN, LW_STD);
  if (GENDER == MALE)
    fprintf (log, "LW Average chrX Read Depth: %f, Standard Deviation: %f\n", LW_MEAN_X, LW_STD_X);

  fprintf (log, "SW Average Read Depth: %f, Standard Deviation: %f\n", SW_MEAN, SW_STD);
  if (GENDER == MALE)
    fprintf (log, "SW Average chrX Read Depth: %f, Standard Deviation: %f\n", SW_MEAN_X, SW_STD_X);

  fclose(log);

  freeMem(gclookup, sizeof(float) * GC_BIN);
  freeMem(gclookup_x, sizeof(float) * GC_BIN);

}



void readDepth(char *depthFile){

  FILE *binDepth;
  int i, j;
  int retVal;

  double lw_total;
  double sw_total;
  double cw_total;

  int lw_cnt;
  int sw_cnt;
  int cw_cnt;

  int isMagicNum;


  binDepth = my_fopen(depthFile, "r", 0);

  retVal = fread(&isMagicNum, sizeof(isMagicNum), 1, binDepth);


  if (isMagicNum == magicNumDepth)
    ISNORMALIZED = 1;
  else if (isMagicNum == magicNum)
    fprintf(stdout, "Read depth file is not yet normalized.\n");
  else
    print_error("Read depth file seems to be invalid or corrupt.\n");



  lw_total = 0.0;
  sw_total = 0.0;
  cw_total = 0.0;

  lw_cnt   = 0;
  sw_cnt   = 0;
  cw_cnt   = 0;

  float lwdepth;
  float swdepth;
  float cwdepth;

  /* read LW */

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->lw_cnt; j++){
      retVal = fread(&(lwdepth), sizeof(chromosomes[i]->lw[j].depth), 1, binDepth);
      chromosomes[i]->lw[j].depth += lwdepth;
      lw_total += lwdepth;
//      lw_total += chromosomes[i]->lw[j].depth;
      chromosomes[i]->lw[j].isControl = 1;
      lw_cnt++;
    }
  }

  /* read SW */

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->sw_cnt; j++){
      retVal = fread(&(swdepth), sizeof(chromosomes[i]->sw[j].depth), 1, binDepth);
			chromosomes[i]->sw[j].depth += swdepth;
      sw_total += swdepth;
//sw_total += chromosomes[i]->sw[j].depth;
      chromosomes[i]->sw[j].isControl = 1;
      sw_cnt++;
    }
  }
// Bu forlar
  /* read CW */

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->cw_cnt; j++){
      retVal = fread(&(cwdepth), sizeof(chromosomes[i]->cw[j].depth), 1, binDepth);
			chromosomes[i]->cw[j].depth += cwdepth;
//     cw_total += chromosomes[i]->cw[j].depth;
      chromosomes[i]->cw[j].isControl = 1;
      cw_total +=cwdepth; 
      cw_cnt++;
    }
  }

  LW_MEAN = lw_total / lw_cnt;
  SW_MEAN = sw_total / sw_cnt;
  CW_MEAN = cw_total / cw_cnt;

  fprintf(stdout, "[OK] depth file %s is loaded.\n", depthFile);
  if (VERBOSE){
    fprintf(stdout, "LW_MEAN: %f\tSW_MEAN: %f\tCW_MEAN:%f\n",  LW_MEAN, SW_MEAN, CW_MEAN);
    if (retVal == 0)
      fprintf(stderr, "There was potentially an error in fread.\n");
  }
  fclose(binDepth);
}


void dump_text_windows(char *fname, enum WINDOWTYPE wt){

  FILE *txtDepth;
  int i, j;

  txtDepth = my_fopen(fname, "w", 0);

  fprintf(txtDepth, "#%s\t%s\t%s\t%s\t%s\t%s\n\n", "CHROM", "START", "END", "GC\%", "READ_DEPTH", "IS_CONTROL");

  switch (wt){
  case LW:
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->lw_cnt; j++)
        if (chromosomes[i]->lw[j].isControl == 1)
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tY\n", chromosomes[i]->name, chromosomes[i]->lw[j].start, chromosomes[i]->lw[j].end, chromosomes[i]->lw[j].gc, chromosomes[i]->lw[j].depth);
	else
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tN\n", chromosomes[i]->name, chromosomes[i]->lw[j].start, chromosomes[i]->lw[j].end, chromosomes[i]->lw[j].gc, chromosomes[i]->lw[j].depth);
    }
    break;

  case SW:
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->sw_cnt; j++)
        if (chromosomes[i]->sw[j].isControl == 1)
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tY\n", chromosomes[i]->name, chromosomes[i]->sw[j].start, chromosomes[i]->sw[j].end, chromosomes[i]->sw[j].gc, chromosomes[i]->sw[j].depth);
	else
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tN\n", chromosomes[i]->name, chromosomes[i]->sw[j].start, chromosomes[i]->sw[j].end, chromosomes[i]->sw[j].gc, chromosomes[i]->sw[j].depth);
    }
    break;

  case CW:
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->cw_cnt; j++)
        if (chromosomes[i]->cw[j].isControl == 1)
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tY\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end, chromosomes[i]->cw[j].gc, chromosomes[i]->cw[j].depth);
	else
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tN\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end, chromosomes[i]->cw[j].gc, chromosomes[i]->cw[j].depth);
    }
    break;
  }

  fclose(txtDepth);

}



void  print_copy_numbers(char *fname){

  FILE *txtDepth;
  int i, j;
  float copy_num;

  txtDepth = my_fopen(fname, "w", 0);

  fprintf(txtDepth, "#%s\t%s\t%s\t%s\t%s\n\n", "CHROM", "START", "END", "GC\%", "COPYNUMBER");

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->cw_cnt; j++){
      if (isAutosome(chromosomes[i]->cw, j, i))
	copy_num = (chromosomes[i]->cw[j].depth / CW_MEAN) * 2;
      else
	copy_num = chromosomes[i]->cw[j].depth / CW_MEAN_X;
 
      if (GENDER == FEMALE && (strstr(chromosomes[i]->name, "chrY")))
	continue;

      fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end, chromosomes[i]->cw[j].gc, copy_num);
    }
  }


  fclose(txtDepth);


}


void conc_depth(char **c_depths, int fileCount, char *depthFile){
  int i;

  if (GENOME_CONF == NULL)
    print_error("Select genome configuration file (input) through the -conf parameter.\n");
  if (c_depths == NULL) 
    print_error("Select read depth input files through the -concdepth parameter.\n");  
  if (depthFile == NULL)
    print_error("Select read depth output file through the -depth parameter.\n");

  
  loadRefConfig(GENOME_CONF);

  for(i = 0; i < fileCount; i++)
    readDepth(c_depths[i]);
  
  saveDepth(depthFile, 1);
}

