/*
 * Various data structures for statistics
 * Copyright (C) 2010-2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 *     Various data structure read/write/report routines.
 *     Defined in "hca.h"
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "hca.h"
#include "data.h"
#include "stats.h"
#include "diag.h"
#include "pctl.h"
#include "probs.h"
#include "pmi.h"
#include "ehash.h"

/*
 *    computes various quality measures,
 *    sending them both to the log and to the ".par" file
 */
void hca_report(char *resstem, char *stem, int ITER, int procs,
		enum GibbsType fix, int showlike, int nopar) {
    double scale = 1;
    char *fname = NULL;
    double logprob;
    FILE *fp = NULL;
    if ( nopar==0 ) {
      fname = yap_makename(resstem,".par");
      fp = fopen(fname, "a");
      if ( !fp )
        yap_sysquit("Cannot open output '%s' file:", fname);
    }
    if ( ddP.hold_all==0 && ddP.phi==NULL ) {
      if ( !showlike )
	scale = -M_LOG2E /ddN.NT;
      logprob = likelihood();
      if ( fp )
	fprintf(fp, "logperptrain = %lf\n", scale * logprob);
      yap_message("log_2(train perp) = %lf\n", scale * logprob);    
    }
    if ( ddN.TEST>0 ) { 
#ifdef EXPERIMENTAL
      if ( ddP.lrsiter>0 ) {
	logprob = lp_test_LRS();
	yap_message("log_2(test perpLRS) = %lf\n", -M_LOG2E * logprob);
        if ( fp )
          fprintf(fp, "logperpLRStest_%d = %lf\n", ITER, -M_LOG2E * logprob);
      }
#endif
      if ( ddP.mltiter>0 ) {
	char *teststr = fix==GibbsHold?"Hold":"ML";
	logprob = lp_test_ML(procs, fix);
	yap_message("log_2(test perp%s) = %lf\n", teststr, -M_LOG2E * logprob);
	if ( fp )
          fprintf(fp, "logperp%stest_%d = %lf\n", 
                  teststr, ITER, -M_LOG2E * logprob);
      }
      if ( ddD.c && ddP.prditer>0 ) {
	logprob = lp_test_Pred(resstem);
	yap_message("test accuracy = %lf\n", logprob);
	if ( fp )
          fprintf(fp, "test accuracy = %lf\n", logprob);
      }
    }
    if ( fp )
      fclose(fp);
    if ( fname ) 
      free(fname);
}

void diag_alloc() {
  tprob_null();
  ddG.dophi = 0;
  ddG.phi_cnt = 0;
  ddG.iscodeword = NULL;
  ddG.words = NULL;
  ddG.n_words = 0; 
  ddG.didcode = 0; 
  ddG.docode = 0;
  ddG.code = NULL;
}

/*
 *  allocation of various data statistics, matrices, vectors
 */
void hca_alloc() {
  ddS.z = u16vec(ddN.N);
  if ( !ddS.z )
    yap_quit("Cannot allocate memory for data\n");
  ddS.Tlife = ddS.TDt = NULL;
  ddS.Tdt = NULL;
  ddS.TWt = NULL;
  ddS.Twt = NULL;
  ddS.TwT = NULL;
  ddS.TWt = NULL;
  if ( ddP.PYbeta && ddP.phi==NULL ) {
    ddS.Twt = u16mat(ddN.W,ddN.T);
    ddS.TwT = u32vec(ddN.W);
    ddS.TWt = u32vec(ddN.T);
  } 
  if ( ddP.PYalpha ) {
    ddS.Tdt = u16mat(ddN.D,ddN.T);
    ddS.TDt = u32vec(ddN.T);
    ddS.Tlife = u32vec(ddN.T);
  }
  ddS.NdT = u16vec(ddN.D);
  ddS.Ndt = u16mat(ddN.D,ddN.T);
  if ( ddP.phi==NULL ) 
    ddS.Nwt = u32mat(ddN.W,ddN.T);
  ddS.NWt  = u32vec(ddN.T);
  sparsemap_null();
  tprob_null();
}

void hca_free() {
  /*
   *  free
   */
  free(ddS.NWt);
  free(ddS.z);
  if ( ddS.Nwt ) {
    free(ddS.Nwt[0]); free(ddS.Nwt);
  }
  free(ddS.NdT);
  u16mat_free(ddS.Ndt, ddN.D, ddN.T);
  if ( ddP.PYbeta && ddP.phi==NULL ) {
    free(ddS.TwT);
    free(ddS.TWt);
    free(ddS.Twt[0]);  free(ddS.Twt);
  } 
  if ( ddP.PYalpha ) {
    free(ddS.TDt);
    free(ddS.Tlife);
    u16mat_free(ddS.Tdt, ddN.D, ddN.T);

  }  
  sparsemap_free();
  tprob_free();
}

/*
 *   Build a matrix doc-frequency for words.
 *   If topic<0, do for all training words,
 *   otherwise only for words of topic
 */
uint32_t **hca_dfmtx(uint32_t *words, int n_words, int topic) {
  uint32_t **mtx = u32mat(n_words,n_words);
  int i, l;
  char *win = calloc(n_words, 1);  /* boolean, stores if word in doc */
  uint32_t *docwords = u32vec(n_words);   /* stores index of word in doc */
  ehash_t hp;

  /*
   *  make hash table
   */
  if ( !win || !docwords || !mtx || ehash_init(&hp, n_words*2) )
    yap_quit("Cannot allocate table in hca_dfmtx()\n");
  for (i=0; i<n_words; i++) 
    ehash_addw(&hp, words[i], i);
  /*
   *  run through docs 
   */
  for (i=0; i<ddN.DT; i++) {
    memset(win,0,n_words);
    int n_dw = 0;
    /*
     *    record which are in
     */
    for (l=ddD.NdTcum[i]; l<ddD.NdTcum[i+1]; l++) {
      if ( topic<0 || Z_t(ddS.z[l])==topic ) {
	int i1 = ehash_findw(&hp, ddD.w[l], words);
	if ( i1!=UINT32_MAX ) {
	  if ( !win[i1] ) {
	    win[i1] = 1;
	    docwords[n_dw++] = i1;
	  }
	}
      }
    }
    /*
     *  now update stats
     */
    for (l=0; l<n_dw; l++) {
      int l2;
      int t1 = docwords[l];
      mtx[t1][t1]++;
      for (l2=0; l2<l; l2++) {
	if ( t1<docwords[l2] )
	  mtx[docwords[l2]][t1]++;
	else
	  mtx[t1][docwords[l2]]++;
      }
    }
  }
  free(docwords);
  free(win);
  ehash_free(&hp);
  return mtx;
}

/*
 *    randomise topics for docs in range
 */
void hca_rand_z(int Tinit, int firstdoc, int lastdoc) {
  int i, t;
  int firstword = ddD.NdTcum[firstdoc];
  int lastword = ddD.NdTcum[lastdoc];
  for (i=firstword; i<lastword; i++) {
    t = (int)(Tinit*rng_unit(rngp));
    if ( t>=Tinit ) t = Tinit-1;
    //   zero's r as a side effect
    ddS.z[i] = t;
    //    initally, all doc indicators are on
    Z_setr(ddS.z[i]);
  }
#if 0
  if ( ddP.bdk!=NULL ) {
    dmi_rand(&ddM, firstdoc, lastdoc);
  }
#endif
}

/*
 *   read topic assignments (z's) from file (both t and r)
 *   for documents docstart to (docend-1)
 */
void hca_read_z(char *resstem, int docstart, int docend) {
  FILE *fr;
  int i, t, r;
  int dobd = (ddP.bdk!=NULL);
  char buf[50];
  char *restartfile = yap_makename(resstem, ".zt");
  fr = fopen(restartfile,"r");
  if ( !fr )
    yap_sysquit("restart file '%s' not read\n", restartfile);
  i = 0;
  docstart = ddD.NdTcum[docstart];
  docend = ddD.NdTcum[docend];
  if ( docstart>0 ) {
    for ( ; i<docstart; i++) {
      if ( !fgets(&buf[0], 50, fr) )
	yap_quit("Cannot read %d-th entry from '%s'\n", i, restartfile);
    }
  }
  for (; i<docend; i++) {
    if ( !fgets(&buf[0], 50, fr) )
      yap_quit("Cannot read %d-th entry from '%s'\n", i, restartfile);
    if ( ( dobd && sscanf(buf," %d,%d", &t, &r)!=2 )
	 || ( !dobd && sscanf(buf," %d", &t)!=1 ) )
      yap_quit("Cannot read %d-th entry from '%s'\n", i, restartfile);
    if ( t>=ddN.T || t<0 )
      yap_quit("Illegal t=%d in '%s'\n", t, restartfile);
    ddS.z[i] = t;
    if ( dobd && r )
      Z_setr(ddS.z[i]);
  }
  fclose(fr);
  free(restartfile);
}

void hca_write_z(char *resstem) 
{
  char *fname = yap_makename(resstem, ".zt");
  FILE *fp = fopen(fname,"w");
  int i;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for write\n", fname);
  for (i = 0; i < ddN.N; i++) {
    if ( ddP.bdk!=NULL )
      fprintf(fp, "%u,%u\n", (unsigned)Z_t(ddS.z[i]), (unsigned)Z_r(ddS.z[i]));
    else
      fprintf(fp, "%u\n", (unsigned)Z_t(ddS.z[i]));
  }
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", fname);
  fclose(fp);
  free(fname);
}

/*
 *    zero everything and rebuild entirely from z[] (both t and r)
 *    but only for training docs
 */
void hca_reset_stats(char *resstem, 
		     int restart,  // restarting, read stats from file
		     int zero,     // just zero and return
		     int firstdoc, 
		     int lastdoc  // may be greater the ddN.DT, 
		                  // so mod before use
		     ) {
  int i, t;
  /*
   *  reset data counts, N??, first,
   */
  memset((void*)ddS.NdT, 0, sizeof(ddS.NdT[0])*ddN.D);
  memset((void*)ddS.NWt, 0, sizeof(ddS.NWt[0])*ddN.T);
  for (i=0; i<ddN.D; i++)
    /*   ddS.Ndt not allocated monolithically  */
    memset((void*)ddS.Ndt[i], 0, sizeof(ddS.Ndt[0][0])*ddN.T);
  if ( ddP.phi==NULL )
    memset((void*)ddS.Nwt[0], 0, sizeof(ddS.Nwt[0][0])*ddN.W*ddN.T);
  /*
   *  now reset table count stats
   */
  if ( ddP.PYbeta && ddP.phi==NULL ) {
    memset((void*)ddS.TwT, 0, sizeof(ddS.TwT[0])*ddN.W);
    memset((void*)ddS.TWt, 0, sizeof(ddS.TWt[0])*ddN.T);
    memset((void*)ddS.Twt[0], 0, sizeof(ddS.Twt[0][0])*ddN.W*ddN.T);
  }
  if ( ddP.PYalpha ) {
    memset((void*)ddS.TDt, 0, sizeof(ddS.TDt[0])*ddN.T);
    memset((void*)ddS.Tlife, 0, sizeof(ddS.Tlife[0])*ddN.T);
    for (i=0; i<ddN.D; i++)
      /*   ddS.Tdt not allocated monolithically  */
      memset((void*)ddS.Tdt[i], 0, sizeof(ddS.Tdt[0][0])*ddN.T);
  }
  ddS.TWT = 0;
  ddS.TWTnz = 0;
  ddS.TDT = 0;
  ddS.TDTnz = 0;
  if ( zero )
    return;

  /*
   *  do data counts, N??, first,
   */  
  for (i=firstdoc; i<lastdoc; i++) {
    int l;
    int usei = (i+ddN.DT) % ddN.DT;
    for (l=ddD.NdTcum[usei]; l<ddD.NdTcum[usei+1]; l++) {
      t = Z_t(ddS.z[l]);
      if ( ( ddP.bdk==NULL ) || Z_issetr(ddS.z[l]) ) {
	if ( ddP.phi==NULL ) {
	  ddS.Nwt[ddD.w[l]][t]++;
	  ddS.NWt[t]++;
	}
      }
      ddS.Ndt[usei][t]++;
      ddS.NdT[usei]++;
    }
  }

  if ( ddP.PYbeta && ddP.phi==NULL ) {
    if ( restart ) {
      char *fname = yap_makename(resstem,".twt");
      read_u16sparse(ddN.W,ddN.T,ddS.Twt,fname);
      free(fname);
      /*  
       *  check consistency
       */
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  if ( ((ddS.Nwt[i][t]>0) ^ (ddS.Twt[i][t]>0)) ||
	       ddS.Nwt[i][t]<ddS.Twt[i][t] ) 
	    yap_quit("Inconsistency:  ddS.Nwt[%d][%d]=%d, ddS.Twt[%d][%d]=%d\n",
		     i, t, (int)ddS.Nwt[i][t],i, t, (int)ddS.Twt[i][t]);
	}
      }
     } else {
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  if ( ddS.Nwt[i][t]>0 ) 
	    ddS.Twt[i][t] = 1;
	}
      }
    } 
    for (i=0; i<ddN.W; i++) {
      for (t=0; t<ddN.T; t++) {
	int tt = ddS.Twt[i][t];
	if ( tt>0 ) {
	  ddS.TWt[t] += tt;
	  ddS.TwT[i] += tt;
	}
      }
    }
    for (i=0; i<ddN.W; i++) {
      if ( ddS.TwT[i]>0 )
        ddS.TWTnz++;
    }
    for (t=0; t<ddN.T; t++) {
      ddS.TWT += ddS.TWt[t];
    }
  }
  if ( ddP.PYalpha  ) {
     if ( restart ) {
      char *fname = yap_makename(resstem,".tdt");
      read_u16sparse(ddN.DT,ddN.T,ddS.Tdt,fname);
      free(fname);
      /*  
       *  check consistency
       */
      for (i=firstdoc; i<lastdoc; i++) {
	int usei = (i+ddN.DT) % ddN.DT;
	for (t=0; t<ddN.T; t++) {
	  if ( ((ddS.Ndt[usei][t]>0) ^ (ddS.Tdt[usei][t]>0)) ||
	       (ddS.Tdt[usei][t]>ddS.Ndt[usei][t])) 
	    yap_quit("Inconsistency:  ddS.Ndt[%d][%d]=%d, ddS.Tdt[%d][%d]=%d\n",
		     usei, t, (int)ddS.Ndt[usei][t], usei, t, (int)ddS.Tdt[usei][t]);
	}
      }
    } else {
      for (i=firstdoc; i<lastdoc; i++) {
	int usei = (i+ddN.DT) % ddN.DT;
	for (t=0; t<ddN.T; t++) {
	  if ( ddS.Ndt[usei][t]>0 ) {
	    ddS.Tdt[usei][t] = 1;
	  }
	}
      }
    }
    for (i=firstdoc; i<lastdoc; i++) {
      int usei = (i+ddN.DT) % ddN.DT;
      for (t=0; t<ddN.T; t++) {
        if ( ddS.Ndt[usei][t]>0 ) {
          ddS.Tdt[usei][t] = 1;
          ddS.TDt[t]++;
        }
      }
    }
    for (t=0; t<ddN.T; t++) {
      ddS.TDT += ddS.TDt[t];
      if ( ddS.TDt[t]>0 )
        ddS.TDTnz++;
    }
  }
}

/*
 *  ensure Twt satisfies constraints
 */
void hca_correct_twt()  {
  int w, t;
  
  if ( ddP.PYbeta==0 )
    return;

  ddS.TWT = ddS.TWTnz = 0;
  memset((void*)ddS.TwT, 0, sizeof(ddS.TwT[0])*ddN.W);
  memset((void*)ddS.TWt, 0, sizeof(ddS.TWt[0])*ddN.T);
  for(t = 0; t < ddN.T; t++) {
    for(w = 0; w < ddN.W; w++) {
      if ( ddS.Twt[w][t]>ddS.Nwt[w][t] )
	ddS.Twt[w][t] = ddS.Nwt[w][t];
      if ( ddS.Twt[w][t]==0 && ddS.Nwt[w][t]>0 )
	ddS.Twt[w][t] = 1;
      ddS.TWt[t] += ddS.Twt[w][t];
      ddS.TwT[w] += ddS.Twt[w][t];
    }
    ddS.TWT += ddS.TWt[t];
  }
  for(w = 0; w < ddN.W; w++) {
    if ( ddS.TwT[w]>0 )
      ddS.TWTnz ++;
  }
}
