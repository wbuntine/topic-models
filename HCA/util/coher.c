/*
 * Coherency check
 * Copyright (C) 2011 Nan Ding
 *           (C) 2013-4 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Nan Ding 
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "yap.h"
#include "util.h"

/*
 *  print out the topic topk=10 words. report the PMI score. 
 */
double coherency_check(char *cofile,   /* name of COOC file */
		       char *vocfile,  /* name of vocab file */
		       char *pmifile,  /* name of PMI file */
		       int T,          /* total topics */
		       int W,          /* total words */
		       uint32_t **Nwt, /* counts [W][T] */
		       int topk)
{
  int t,i,j;
  double coherency=0.0;
  uint16_t **Nwt_dispi=u16mat(T,topk);
  uint32_t **Nwt_dispn=u32mat(T,topk);
  uint32_t *Nwt_tmp=u32vec(W);
  uint16_t *Nwt_id=u16vec(W);
  uint32_t x;
  uint16_t id;
  char **ctmp = read_vocab(vocfile,0,W,50);

  FILE *fr;
  fr = fopen(cofile,"w");
  if ( !fr ) 
    yap_sysquit("cofile '%s' not read\n", cofile);

  FILE *fr2;
  fr2 = fopen(pmifile,"r");
  if ( !fr2 ) 
    yap_sysquit("pmifile '%s' not read\n", pmifile);
  
  // for each topic, get the top 10 words by bubble sort of Nwt
  for (t=0;t<T;t++) {
    int tot=0;
    for (i=0;i<W;i++)
      tot += Nwt_tmp[i]=Nwt[i][t];
    for (i=0;i<W;i++)
      Nwt_id[i]=i;
    fprintf(fr, "#%d count=%d\n", t, tot);
    
    for (j=0;j<topk;j++) {
      for (i=1;i<W-j;i++) {
	if (Nwt_tmp[i-1]>Nwt_tmp[i]) {
	  x=Nwt_tmp[i-1];Nwt_tmp[i-1]=Nwt_tmp[i];Nwt_tmp[i]=x;
	  id=Nwt_id[i-1];Nwt_id[i-1]=Nwt_id[i];Nwt_id[i]=id;
	}
      }
      Nwt_dispi[t][j]=Nwt_id[W-1-j];
      Nwt_dispn[t][j]=Nwt_tmp[W-1-j];
      fprintf(fr,"%d %u %u %s\n",
	      t,Nwt_dispi[t][j],Nwt_dispn[t][j],ctmp[Nwt_dispi[t][j]]);
    }
  }
  
  // get the PMI matrix (which is dense). 
  int t1,t2;
  double value;
  double **PMI_mat=dmat(W,W);
  while (fscanf(fr2, "%d", &t1) != EOF ) { 
      if (fscanf(fr2, "%d", &t2) != EOF) {
	  if (fscanf(fr2,"%lf",&value) !=EOF) {
	    PMI_mat[t1-1][t2-1]=value;
	  }
      }
  }
  
  // compute PMI score for each topic
  double coh[T];
  for (t=0;t<T;t++) {
    coh[t]=0;
    for (i=0;i<topk;i++) {
      for (j=i+1;j<topk;j++) {
	value=PMI_mat[Nwt_dispi[t][i]][Nwt_dispi[t][j]];
	coh[t]+=value;
	coherency+=value;
      }
    }
    coh[t]/=(topk-1)*topk/2;
    fprintf(fr, "#%d  coher=%lf\n",t,coh[t]);
  }
  
  coherency/=T*(topk-1)*topk/2;
  
  free(Nwt_dispi[0]);
  free(Nwt_dispi);
  free(Nwt_dispn[0]);
  free(Nwt_dispn);
  free(Nwt_id);
  free(Nwt_tmp);
  free(PMI_mat[0]);
  free(PMI_mat);

  free(ctmp[0]); free(ctmp);

  fclose(fr2);
  fclose(fr);
  
  
  return coherency;
}
