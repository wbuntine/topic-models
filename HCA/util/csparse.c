/*
 * Various data structure read/write/report routines.
 * Copyright (C) 2010-2011 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *     Various data structure read/write/report routines.
 *     Defined in "dca.h"
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"

#define DEBUG_READC
// Read C matrix as a sparse matrix from files
void read_c_sparse(char *cfile, double threshold, int maxdim,
		   int W, sparse_vec *C) {

  FILE *fr;
  int i,t1, t2;
  double value;

  // because the words are indexed from 1. 
  double *self_count=dvec(W);
  uint16_t *idx=u16vec(W);
  double *sum_vec=dvec(W);
  uint16_t *nnz_vec=u16vec(W);
#ifdef DEBUG_READC
  char **ctmp = (char **)read_vocab(cfile,W,50);
  FILE *vr;
  char *vrname = yap_makename(cfile,".out");
  vr = fopen(vrname,"w");
  if ( !vr )
    yap_sysquit("Cannot open file '%s' for write\n", vrname);
#endif
  for (i=0;i<W;i++)
    nnz_vec[i]=0;

  for (i=0;i<W;i++)
    sum_vec[i]=0;

  for (i=0;i<W;i++)
    idx[i]=0;


  fr = fopen(cfile,"r");
  if ( !fr ) 
    yap_sysquit("cfile '%s' not read\n", cfile);
  
  // get the count of each of word
  while (fscanf(fr, "%d", &t1) != EOF ) { 
    if (fscanf(fr, "%d", &t2) != EOF)
      {
	if (fscanf(fr,"%lf",&value) !=EOF)
	  {
	    if (t1==t2)
	      self_count[t1-1]=value;
	  }
      }
  }

  // get the count of co-occurance of two words. and only retain the words-pair that has co-occurance larger than a threshold and compute the sum of each column
  fseek(fr,0,SEEK_SET);
  while (fscanf(fr, "%d", &t1) != EOF ) { 
    if (fscanf(fr, "%d", &t2) != EOF)
      {
	if (fscanf(fr,"%lf",&value) !=EOF)
	  {
	    if (value*value>=threshold*self_count[t1-1]*self_count[t2-1])
	      {
		
		nnz_vec[t1-1]++;
		sum_vec[t1-1]+=value;

		if (t1 != t2)
		  {
		    nnz_vec[t2-1]++;
		    sum_vec[t2-1]+=value;
		  }
	      }	    
	  }
      }
  }

  for (i=0;i<W;i++)
    {
      if (nnz_vec[i]>0)
	{
	  C[i].id=u16vec(nnz_vec[i]);
	  C[i].value=dvec(nnz_vec[i]);
	  C[i].nnz=nnz_vec[i];
	}
    }

  // fill in the C matrix by using the sparser co-occurance matrix with normalization
  fseek(fr,0,SEEK_SET);
  while (fscanf(fr, "%d", &t1) != EOF ) { 
    if (fscanf(fr, "%d", &t2) != EOF)
      {
	if (fscanf(fr,"%lf",&value) !=EOF)
	  {
	    if (value*value>=threshold*self_count[t1-1]*self_count[t2-1])
	      {
		C[t1-1].id[idx[t1-1]]=t2-1;
		C[t1-1].value[idx[t1-1]]=value/sum_vec[t2-1];
		idx[t1-1]++;
#ifdef DEBUG_READC
		    fprintf(vr,"%s %s %f\n", ctmp[t2-1], ctmp[t1-1], value/sum_vec[t2-1]);
#endif

		if (t1 !=t2)
		  {
		    C[t2-1].id[idx[t2-1]]=t1-1;
		    C[t2-1].value[idx[t2-1]]=value/sum_vec[t1-1];
#ifdef DEBUG_READC
		    fprintf(vr,"%s %s %f\n", ctmp[t1-1], ctmp[t2-1], value/sum_vec[t1-1]);
#endif
		    idx[t2-1]++;
		  }
	      }	    
	  }
      }
  }

  fclose(fr);

  free(self_count);
  free(sum_vec);
  free(idx);
  free(nnz_vec);
#ifdef DEBUG_READC
  free(ctmp[0]); free(ctmp);
  fclose(vr);
  free(vrname);
#endif
}


