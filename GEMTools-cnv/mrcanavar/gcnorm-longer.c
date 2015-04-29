#include "gcnorm.h"

int norm_until_converges (enum WINDOWTYPE wt, float *gclookup, float *gclookup_x){
  float max;
  float max_x;
  float min;
  float min_x;

  float p_max;
  float p_max_x;
  float p_min;
  float p_min_x;

  int i, j;
  float mean    = 0.0;
  float mean_x  = 0.0;
  float stdev   = 0.0;
  float stdev_x = 0.0;
  int iter;

  float maxcut;
  float maxcut_x;
  float mincut;
  float mincut_x;

  iter = 1;
  calc_stat(wt, gclookup, gclookup_x, 0, &max, &max_x, &min, &min_x);

  switch(wt){
  case LW:
    mean    = LW_MEAN;
    mean_x  = LW_MEAN_X;
    break;
  case SW:
    mean    = SW_MEAN;
    mean_x  = SW_MEAN_X;
    break;
  case CW:
    mean    = CW_MEAN;
    mean_x  = CW_MEAN_X;
    break;
  }


  if (GENDER == AUTODETECT){
    if (VERBOSE)
      fprintf(stdout, "MEAN: %f\tMEAN_X: %f\n", mean, mean_x);


    if (mean_x / mean < 0.75){ // magical ratio
      GENDER = MALE;

      if (VERBOSE)
	fprintf(stdout, "Autodetect Gender: Male\n");

      i = findChrom("chrX", "", -1);


      switch(wt){

      case CW:
      for (j=0; j<chromosomes[i]->cw_cnt; j++)
	if (chromosomes[i]->cw[j].depth > (float) CW_MEAN_X * 2 || chromosomes[i]->cw[j].depth < (float) CW_MEAN_X / 10.0)
	  chromosomes[i]->cw[j].isControl = 0;
      break;

      case LW:
      for (j=0; j<chromosomes[i]->lw_cnt; j++)
	if (chromosomes[i]->lw[j].depth > (float) LW_MEAN_X * 2 || chromosomes[i]->lw[j].depth < (float) LW_MEAN_X / 10.0)
	  chromosomes[i]->lw[j].isControl = 0;
      break;

      case SW:
      for (j=0; j<chromosomes[i]->sw_cnt; j++)
	if (chromosomes[i]->sw[j].depth > (float) SW_MEAN_X * 2 || chromosomes[i]->sw[j].depth < (float) SW_MEAN_X / 10.0)
	  chromosomes[i]->sw[j].isControl = 0;
      break;

      }

    }
    else{
      GENDER = FEMALE;
      if (VERBOSE)
	fprintf(stdout, "Autodetect Gender: Female\n");
    }
  }

  normalize(wt, gclookup, gclookup_x, &max, &max_x, &min, &min_x);


  do{

    if (VERBOSE){
      fprintf(stdout, "Control regions %s iteration %d ", (wt == LW ? "LW" : (wt == SW ? "SW" : "CW")), iter);
      fflush(stdout);
    }

    p_min   = min;
    p_max   = max;
    p_min_x = min_x;
    p_max_x = max_x;

    calc_stat(wt, gclookup, gclookup_x, 1, &max, &max_x, &min, &min_x);

    switch(wt){
    case LW:
      mean    = LW_MEAN;
      stdev   = LW_STD;
      mean_x  = LW_MEAN_X;
      stdev_x = LW_STD_X;
      break;
    case SW:
      mean    = SW_MEAN;
      stdev   = SW_STD;
      mean_x  = SW_MEAN_X;
      stdev_x = SW_STD_X;
      break;
    case CW:
      mean    = CW_MEAN;
      stdev   = CW_STD;
      mean_x  = CW_MEAN_X;
      stdev_x = CW_STD_X;
      break;
    }

    if (GENDER == FEMALE){
      max_x   = max;
      mean_x  = mean;
      stdev_x = stdev;
      min_x   = min;
    }

    iter++;

    maxcut   = mean + 2.5 * stdev;
    mincut   = mean - 2.5 * stdev;
    maxcut_x = mean_x + 2.5 * stdev_x;
    mincut_x = mean_x - 2.5 * stdev_x;

    if (GENDER != FEMALE){
      maxcut_x = mean_x + 2 * stdev_x;
      mincut_x = mean_x - 2 * stdev_x;
    }

    if (mincut < 0.0) mincut = mean / 10.0;
    if (mincut_x < 0.0) mincut_x = mean_x / 10.0;

    if (maxcut - mean > mean - mincut)
      maxcut = mean + (mean - mincut);

    if (GENDER != FEMALE){
      if (maxcut_x - mean_x > mean_x - mincut_x)
	maxcut_x = mean_x + (mean_x - mincut_x);
    }


    if (VERBOSE){
      fprintf(stdout, "mean: %f\tstdev: %f\tmax: %f (cut: %f)\tmin: %f (cut: %f) \n", mean, stdev, max, maxcut, min, mincut);
      if (GENDER != FEMALE)
	fprintf(stdout, "mean_x: %f\tstdev_x: %f\tmax_x: %f (cut: %f)\tmin_x: %f (cut: %f)\n", mean_x, stdev_x, max_x, maxcut_x, min_x, mincut_x);
    }

    if (p_min == min && p_max == max && p_min_x == min_x && p_max_x == max_x)
      break;

  } while (max >= maxcut || max_x >= maxcut_x || min <= mincut || min_x <= mincut_x);


  fprintf(stdout, "%s Normalization completed.\n",  (wt == LW ? "LW" : (wt == SW ? "SW" : "CW")));

  return 1;
}


void normalize (enum WINDOWTYPE wt, float *gclookup, float *gclookup_x, float *max, float *max_x, float *min, float *min_x){

  int i;

  for (i=0; i<num_chrom; i++)
    norm_wins(i, wt, gclookup, gclookup_x);

  calc_stat(wt, gclookup, gclookup_x, 0, max, max_x, min, min_x);

}


void calc_stat(enum WINDOWTYPE wt, float *gclookup, float *gclookup_x, char doClean, float *_max, float *_max_x, float *_min, float *_min_x){
  int i, j, k;

  float lw_var;
  float sw_var;
  float cw_var;

  float lw_total;
  float sw_total;
  float cw_total;

  float lw_var_x;
  float sw_var_x;
  float cw_var_x;

  float lw_total_x;
  float sw_total_x;
  float cw_total_x;

  int win_cnt;
  float max;

  int win_cnt_x;
  float max_x;

  float min;
  float min_x;

  int gc_index;

  float gc_total[GC_BIN]; // total depth by GC
  int gc_wincount[GC_BIN]; // count of windows with the given GC content

  float gc_total_x[GC_BIN]; // total depth by GC
  int gc_wincount_x[GC_BIN]; // count of windows with the given GC content

  float MEAN   = 0.0;
  float MEAN_X = 0.0;


  float maxcut, mincut;
  float maxcut_x, mincut_x;

  lw_var = 0.0;
  sw_var = 0.0;
  cw_var = 0.0;

  lw_total = 0.0;
  sw_total = 0.0;
  cw_total = 0.0;

  win_cnt = 0;

  lw_var_x = 0.0;
  sw_var_x = 0.0;
  cw_var_x = 0.0;

  lw_total_x = 0.0;
  sw_total_x = 0.0;
  cw_total_x = 0.0;

  win_cnt_x = 0;

  for (i=0; i<GC_BIN; i++){
    gc_total[i]    = 0.0;
    gclookup[i]    = 0.0;
    gc_wincount[i] = 0;

    gc_total_x[i]    = 0.0;
    gclookup_x[i]    = 0.0;
    gc_wincount_x[i] = 0;
  }

  max = 0.0;
  max_x = 0.0;

  min = FLT_MAX;
  min_x = FLT_MAX;

  switch (wt){
  case LW:

    if (doClean){

      for (i=0; i<num_chrom; i++){
	for (j=0; j<chromosomes[i]->lw_cnt; j++){

	  if (chromosomes[i]->lw[j].isControl == 1){

	    /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	    if (isAutosome(NULL, -1, i)){
	      maxcut = LW_MEAN + 2.5 * LW_STD;
	      mincut = LW_MEAN - 2.5 * LW_STD;
	      if (mincut < LW_MEAN / 10.0) mincut = LW_MEAN / 10.0;

	      if (maxcut - LW_MEAN > LW_MEAN - mincut)
		maxcut = LW_MEAN + (LW_MEAN - mincut);

	      if (chromosomes[i]->lw[j].depth > maxcut || chromosomes[i]->lw[j].depth < mincut){

		/* Remove this window and its neighbors from controls */
		chromosomes[i]->lw[j].isControl = 0;

		k = j;

		while (k > 0 && chromosomes[i]->lw[k-1].end >= chromosomes[i]->lw[j].start){
		  chromosomes[i]->lw[--k].isControl = 0;
		}

		k = j;

		while (k < chromosomes[i]->lw_cnt-1 && chromosomes[i]->lw[k+1].start <= chromosomes[i]->lw[j].end){
		  chromosomes[i]->lw[++k].isControl = 0;
		}

	      }

	      else{
		lw_total += chromosomes[i]->lw[j].depth;
		win_cnt++;
	      }
	    }  // if AUTOSOME


	    /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	    else{
	      maxcut_x = LW_MEAN_X + 2 * LW_STD_X;
	      mincut_x = LW_MEAN_X - 2 * LW_STD_X;
	      if (mincut_x < LW_MEAN_X / 10.0) mincut_x = LW_MEAN_X / 10.0;

	      if (chromosomes[i]->lw[j].depth > maxcut_x || chromosomes[i]->lw[j].depth < mincut_x){

		/* Remove this window and its neighbors from controls */
		chromosomes[i]->lw[j].isControl = 0;

		k = j;

		while (k > 0 && chromosomes[i]->lw[k-1].end >= chromosomes[i]->lw[j].start){
		  chromosomes[i]->lw[--k].isControl = 0;
		}

		k = j;

		while (k < chromosomes[i]->lw_cnt-1 && chromosomes[i]->lw[k+1].start <= chromosomes[i]->lw[j].end){
		  chromosomes[i]->lw[++k].isControl = 0;
		}
	      }

	      else{
		lw_total_x += chromosomes[i]->lw[j].depth;
		win_cnt_x++;
	      }
	    }  // if chrX

	  }

	}
      }

      LW_MEAN_X = lw_total_x / win_cnt_x;
      LW_MEAN = lw_total / win_cnt;

    }  // do clean

    lw_total   = 0.0;
    win_cnt    = 0;
    lw_total_x = 0.0;
    win_cnt_x  = 0;

    for (i=0; i<num_chrom; i++){
      for (j=0; j<chromosomes[i]->lw_cnt; j++){

	if (chromosomes[i]->lw[j].isControl == 1){

	  /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	  if (isAutosome(NULL, -1, i)){
	    lw_var += (LW_MEAN - chromosomes[i]->lw[j].depth) * (LW_MEAN - chromosomes[i]->lw[j].depth);
	    win_cnt++;

	    lw_total += chromosomes[i]->lw[j].depth;
	    gc_index = chromosomes[i]->lw[j].gc * GC_BIN;
	    gc_total[gc_index] += chromosomes[i]->lw[j].depth;
	    gc_wincount[gc_index]++;

	    if (chromosomes[i]->lw[j].depth > max)
	      max = chromosomes[i]->lw[j].depth;
	    if (chromosomes[i]->lw[j].depth < min)
	      min = chromosomes[i]->lw[j].depth;

	  } // AUTOSOMES

	  /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	  else{ // if (strstr(chromosomes[i]->name, "chrX") && GENDER != FEMALE){

	    lw_var_x += (LW_MEAN_X - chromosomes[i]->lw[j].depth) * (LW_MEAN_X - chromosomes[i]->lw[j].depth);
	    win_cnt_x++;

	    lw_total_x += chromosomes[i]->lw[j].depth;
	    gc_index = chromosomes[i]->lw[j].gc * GC_BIN;
	    gc_total_x[gc_index] += chromosomes[i]->lw[j].depth;
	    gc_wincount_x[gc_index]++;

	    if (chromosomes[i]->lw[j].depth > max_x)
	      max_x = chromosomes[i]->lw[j].depth;
	    if (chromosomes[i]->lw[j].depth < min_x)
	      min_x = chromosomes[i]->lw[j].depth;

	  } // chrX


	} // if control
      } // outer for
    }     //    if (!doClean)

    LW_MEAN_X = lw_total_x / win_cnt_x;
    LW_STD_X = sqrt(lw_var_x / win_cnt_x);

    LW_MEAN = lw_total / win_cnt;
    LW_STD = sqrt(lw_var / win_cnt);

    MEAN = LW_MEAN;
    MEAN_X = LW_MEAN_X;

    break;

  case SW:


    if (doClean){
      for (i=0; i<num_chrom; i++){
	for (j=0; j<chromosomes[i]->sw_cnt; j++){

	  if (chromosomes[i]->sw[j].isControl == 1){

	    /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	    if (isAutosome(NULL, -1, i)){
	      maxcut = SW_MEAN + 2.5 * SW_STD;
	      mincut = SW_MEAN - 2.5 * SW_STD;
	      if (mincut < SW_MEAN / 10.0) mincut = SW_MEAN / 10.0;

	      if (maxcut - SW_MEAN > SW_MEAN - mincut)
		maxcut = SW_MEAN + (SW_MEAN - mincut);

	      if (chromosomes[i]->sw[j].depth > maxcut || chromosomes[i]->sw[j].depth < mincut){

		/* Remove this window and its neighbors from controls */
		chromosomes[i]->sw[j].isControl = 0;

		k = j;

		while (k > 0 && chromosomes[i]->sw[k-1].end >= chromosomes[i]->sw[j].start){
		  chromosomes[i]->sw[--k].isControl = 0;
		}

		k = j;

		while (k < chromosomes[i]->sw_cnt-1 && chromosomes[i]->sw[k+1].start <= chromosomes[i]->sw[j].end){
		  chromosomes[i]->sw[++k].isControl = 0;
		}

	      }

	      else{
		sw_total += chromosomes[i]->sw[j].depth;
		win_cnt++;
	      }
	    }  // if AUTOSOME


	    /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	    else{

	      maxcut_x = SW_MEAN_X + 2 * SW_STD_X;
	      mincut_x = SW_MEAN_X - 2 * SW_STD_X;
	      if (mincut_x < SW_MEAN_X / 10.0) mincut_x = SW_MEAN_X / 10.0;

	      if (chromosomes[i]->sw[j].depth > maxcut_x || chromosomes[i]->sw[j].depth < mincut_x){

		/* Remove this window and its neighbors from controls */
		chromosomes[i]->sw[j].isControl = 0;

		k = j;

		while (k > 0 && chromosomes[i]->sw[k-1].end >= chromosomes[i]->sw[j].start){
		  chromosomes[i]->sw[--k].isControl = 0;
		}

		k = j;

		while (k < chromosomes[i]->sw_cnt-1 && chromosomes[i]->sw[k+1].start <= chromosomes[i]->sw[j].end){
		  chromosomes[i]->sw[++k].isControl = 0;
		}
	      }

	      else{
		sw_total_x += chromosomes[i]->sw[j].depth;
		win_cnt_x++;
	      }
	    }  // if chrX

	  }

	}
      }

      SW_MEAN_X = sw_total_x / win_cnt_x;
      SW_MEAN = sw_total / win_cnt;

    }  // do clean

    sw_total   = 0.0;
    win_cnt    = 0;
    sw_total_x = 0.0;
    win_cnt_x  = 0;

    for (i=0; i<num_chrom; i++){
      for (j=0; j<chromosomes[i]->sw_cnt; j++){

	if (chromosomes[i]->sw[j].isControl == 1){

	  /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	  if (isAutosome(NULL, -1, i)){
	    sw_var += (SW_MEAN - chromosomes[i]->sw[j].depth) * (SW_MEAN - chromosomes[i]->sw[j].depth);
	    win_cnt++;

	    sw_total += chromosomes[i]->sw[j].depth;
	    gc_index = chromosomes[i]->sw[j].gc * GC_BIN;
	    gc_total[gc_index] += chromosomes[i]->sw[j].depth;
	    gc_wincount[gc_index]++;

	    if (chromosomes[i]->sw[j].depth > max)
	      max = chromosomes[i]->sw[j].depth;
	    if (chromosomes[i]->sw[j].depth < min)
	      min = chromosomes[i]->sw[j].depth;

	  } // AUTOSOMES

	  /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	  // if (strstr(chromosomes[i]->name, "chrX") && GENDER != FEMALE){
	  else{
	    sw_var_x += (SW_MEAN_X - chromosomes[i]->sw[j].depth) * (SW_MEAN_X - chromosomes[i]->sw[j].depth);
	    win_cnt_x++;

	    sw_total_x += chromosomes[i]->sw[j].depth;
	    gc_index = chromosomes[i]->sw[j].gc * GC_BIN;
	    gc_total_x[gc_index] += chromosomes[i]->sw[j].depth;
	    gc_wincount_x[gc_index]++;

	    if (chromosomes[i]->sw[j].depth > max_x)
	      max_x = chromosomes[i]->sw[j].depth;
	    if (chromosomes[i]->sw[j].depth < min_x)
	      min_x = chromosomes[i]->sw[j].depth;

	  } // chrX


	} // if control
      } // outer for
    }     //    if (!doClean)

    SW_MEAN_X = sw_total_x / win_cnt_x;
    SW_STD_X = sqrt(sw_var_x / win_cnt_x);

    SW_MEAN = sw_total / win_cnt;
    SW_STD = sqrt(sw_var / win_cnt);

    MEAN = SW_MEAN;
    MEAN_X = SW_MEAN_X;


    break;

  case CW:

    if (doClean){
      for (i=0; i<num_chrom; i++){
	for (j=0; j<chromosomes[i]->cw_cnt; j++){


	  if (chromosomes[i]->cw[j].isControl == 1){

	    /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */

	    if (isAutosome(NULL, -1, i)){
	      maxcut = CW_MEAN + 2.5 * CW_STD;
	      mincut = CW_MEAN - 2.5 * CW_STD;
	      if (mincut < CW_MEAN / 10.0) mincut = CW_MEAN / 10.0;

	      if (maxcut - CW_MEAN > CW_MEAN - mincut)
		maxcut = CW_MEAN + (CW_MEAN - mincut);

	      if (chromosomes[i]->cw[j].depth > maxcut || chromosomes[i]->cw[j].depth < mincut){

		/* Remove this window and its neighbors from controls */
		chromosomes[i]->cw[j].isControl = 0;

		if (j > 0)
		  chromosomes[i]->cw[j-1].isControl = 0;
		if (j < chromosomes[i]->cw_cnt-1)
		  chromosomes[i]->cw[j+1].isControl = 0;
	      }
	      else{
		cw_total += chromosomes[i]->cw[j].depth;
		win_cnt++;
	      }
	    }  // if AUTOSOME


	    /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	    //if (strstr(chromosomes[i]->name, "chrX") && GENDER != FEMALE){
	    else{
	      maxcut_x = CW_MEAN_X + 2 * CW_STD_X;
	      mincut_x = CW_MEAN_X - 2 * CW_STD_X;
	      if (mincut_x < CW_MEAN_X / 10.0) mincut_x = CW_MEAN_X / 10.0;

	      if (chromosomes[i]->cw[j].depth > maxcut_x || chromosomes[i]->cw[j].depth < mincut_x){

		/* Remove this window and its neighbors from controls */
		chromosomes[i]->cw[j].isControl = 0;

		if (j > 0)
		  chromosomes[i]->cw[j-1].isControl = 0;
		if (j < chromosomes[i]->cw_cnt-1)
		  chromosomes[i]->cw[j+1].isControl = 0;
	      }
	      else{
		cw_total_x += chromosomes[i]->cw[j].depth;
		win_cnt_x++;
	      }
	    }  // if chrX

	  }

	}
      }


      CW_MEAN_X = cw_total_x / win_cnt_x;

      CW_MEAN = cw_total / win_cnt;

    }  // do clean

    cw_total   = 0.0;
    win_cnt    = 0;
    cw_total_x = 0.0;
    win_cnt_x  = 0;

    for (i=0; i<num_chrom; i++){
      for (j=0; j<chromosomes[i]->cw_cnt; j++){

	if (chromosomes[i]->cw[j].isControl == 1){

	  /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */

	  if (isAutosome(NULL, -1, i)){
	    cw_var += (CW_MEAN - chromosomes[i]->cw[j].depth) * (CW_MEAN - chromosomes[i]->cw[j].depth);
	    win_cnt++;

	    cw_total += chromosomes[i]->cw[j].depth;
	    gc_index = chromosomes[i]->cw[j].gc * GC_BIN;
	    gc_total[gc_index] += chromosomes[i]->cw[j].depth;
	    gc_wincount[gc_index]++;

	    if (chromosomes[i]->cw[j].depth > max)
	      max = chromosomes[i]->cw[j].depth;
	    if (chromosomes[i]->cw[j].depth < min)
	      min = chromosomes[i]->cw[j].depth;

	  } // AUTOSOMES

	  /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */

	  else{
	    cw_var_x += (CW_MEAN_X - chromosomes[i]->cw[j].depth) * (CW_MEAN_X - chromosomes[i]->cw[j].depth);
	    win_cnt_x++;

	    cw_total_x += chromosomes[i]->cw[j].depth;
	    gc_index = chromosomes[i]->cw[j].gc * GC_BIN;
	    gc_total_x[gc_index] += chromosomes[i]->cw[j].depth;
	    gc_wincount_x[gc_index]++;

	    if (chromosomes[i]->cw[j].depth > max_x)
	      max_x = chromosomes[i]->cw[j].depth;
	    if (chromosomes[i]->cw[j].depth < min_x)
	      min_x = chromosomes[i]->cw[j].depth;

	  } // chrX


	} // if control
      } // outer for
    }     //    if (!doClean)


    CW_MEAN_X = cw_total_x / win_cnt_x;
    CW_STD_X = sqrt(cw_var_x / win_cnt_x);

    CW_MEAN = cw_total / win_cnt;
    CW_STD = sqrt(cw_var / win_cnt);

    MEAN = CW_MEAN;
    MEAN_X = CW_MEAN_X;

    break;
  }


  /* calculate the gclookup table */

  for (i=0; i<GC_BIN; i++){

    /* AUTOSOMES */

    if (gc_wincount[i] == 0 || gc_total[i] == 0){
      j=i-1;

      while (j >= 0 && gc_total[j] == 0)
	j--;

      if (gc_total[j] == 0 || gc_wincount[j] == 0){
	j=i+1;
	while (j < GC_BIN && gc_total[j] == 0) j++;
      }

      gc_total[i]    = gc_total[j];
      gc_wincount[i] = gc_wincount[j];
    }

    if (gc_total[i] != 0.0){
      gclookup[i] = gc_total[i] / (float) gc_wincount[i];

      if (MULTGC){

	gclookup[i] = MEAN / gclookup[i];

	if (gclookup[i] > MAX_GC_CORR)
	  gclookup[i] = MAX_GC_CORR;
	else if (gclookup[i] < MIN_GC_CORR)
	  gclookup[i] = MIN_GC_CORR;
      }

    }



    /* chrX */

    if (GENDER != FEMALE){

      if (gc_wincount_x[i] == 0 || gc_total_x[i] == 0){
	j=i-1;

	while (j >= 0 && gc_total_x[j] == 0)
	  j--;

	if (gc_total_x[j] == 0 || gc_wincount_x[j] == 0){
	  j=i+1;
	  while (j < GC_BIN && gc_total_x[j] == 0) j++;
	}

	gc_total_x[i]    = gc_total_x[j];
	gc_wincount_x[i] = gc_wincount_x[j];
      }

      if (gc_total_x[i] != 0.0){
	gclookup_x[i] = gc_total_x[i] / (float) gc_wincount_x[i];

	if (MULTGC){
	  gclookup_x[i] = MEAN_X / gclookup_x[i];
	  if (gclookup_x[i] > MAX_GC_CORR)
	    gclookup_x[i] = MAX_GC_CORR;
	  else if (gclookup_x[i] < MIN_GC_CORR)
	    gclookup_x[i] = MIN_GC_CORR;
	}

      }

    }


  }


  *_max   = max;
  *_max_x = max_x;

  *_min   = min;
  *_min_x = min_x;


}


void norm_wins(int chrom_id, enum WINDOWTYPE wt, float *gclookup, float *gclookup_x){
  int j;
  int gc_index;
  float new_depth;
  struct window *win = NULL;
  int win_cnt        = 0;
  float mean         = 0.0;
  float mean_x       = 0.0;

  switch (wt){
  case LW:
    win     = chromosomes[chrom_id]->lw;
    win_cnt = chromosomes[chrom_id]->lw_cnt;
    mean    = LW_MEAN; 
    mean_x  = LW_MEAN_X;
    break;
  case SW:
    win     = chromosomes[chrom_id]->sw;
    win_cnt = chromosomes[chrom_id]->sw_cnt;
    mean    = SW_MEAN; 
    mean_x  = SW_MEAN_X;
    break;
  case CW:
    win     = chromosomes[chrom_id]->cw;
    win_cnt = chromosomes[chrom_id]->cw_cnt;
    mean    = CW_MEAN; 
    mean_x  = CW_MEAN_X;
    break;
  }

  for (j=0; j<win_cnt; j++){
    gc_index  = win[j].gc * (float) GC_BIN;
    if (isAutosome(win, j, chrom_id)){
      if (MULTGC)
	new_depth = gclookup[gc_index] * win[j].depth;
      else
	new_depth = win[j].depth - (gclookup[gc_index] - mean);
      }
    else{
      if (MULTGC)
	new_depth = gclookup_x[gc_index] * win[j].depth;
      else
	new_depth = win[j].depth - (gclookup_x[gc_index] - mean_x);
    }
    
    if (new_depth < 0) new_depth = 0;
    
    win[j].depth = new_depth;
  }
  
}

