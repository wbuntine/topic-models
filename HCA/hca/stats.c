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
#include "check.h"
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
        double perp;
	logprob = lp_test_ML(procs, fix);
	perp = exp(-logprob);
	yap_message("log_2(test perp%s) = %lf (%lf)\n", teststr, -M_LOG2E * logprob, perp);
	if ( fp )
          fprintf(fp, "logperp%stest_%d = %lf (%lf)\n", 
                  teststr, ITER, -M_LOG2E * logprob, perp);
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
  ddS.UN = NULL;
#ifdef NG_SCALESTATS
  ddS.NGscalestats = NULL;
#endif
  ddS.sparse = NULL;
  ddS.sparseD = NULL;
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
  if ( ddP.PYalpha==H_NG ) {
    if ( ddP.ngs0 || ddP.ngs1 ) {
      int i;
      unsigned incr;
      ddS.sparseD = malloc(ddN.T*sizeof(*ddS.sparse));
      ddS.sparse = malloc(ddN.D*sizeof(*ddS.sparse));
      incr = M_bitveclen();
      ddS.sparse[0] = malloc(ddN.D*sizeof(**ddS.sparse)*incr);
      for (i=1; i<ddN.D; i++)
        ddS.sparse[i] = ddS.sparse[i-1] + incr;
    }
    ddS.UN = malloc(ddN.D*sizeof(*ddS.UN));
#ifdef NG_SCALESTATS
    ddS.NGscalestats = malloc(ddN.T*sizeof(*ddS.NGscalestats));
#endif
  }
  ddS.NdT = u16vec(ddN.D);
  ddS.Ndt = u16mat(ddN.D,ddN.T);
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
  if ( ddS.UN )  free(ddS.UN);
#ifdef NG_SCALESTATS
  if ( ddS.NGscalestats )  free(ddS.NGscalestats);
#endif
  if ( ddS.sparse ) {
    free(ddS.sparse[0]);
    free(ddS.sparse);
    free(ddS.sparseD);
  }
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

void hca_merge_stats(int k1, int k2,  uint16_t *Tdt,  uint16_t *Twt) {
  int d, i, w;
  int32_t Tdiff1=0;

  for (i=0; i<ddN.N; i++)
    if ( Z_t(ddS.z[i])==k2 ) {
      Z_sett(ddS.z[i],k1);
    }

  /*
   *  reset Alpha side stats
   */
  for (d=0; d<ddN.D; d++) {
    ddS.Ndt[d][k1] += ddS.Ndt[d][k2];
    ddS.Ndt[d][k2] = 0;
  }
  if ( ddP.PYalpha ) {
    Tdiff1 = 0;
    ddS.TDT -= ddS.TDt[k1] + ddS.TDt[k2];
    for (d=0; d<ddN.DT; d++) {
      Tdiff1 += Tdt[d];
      ddS.Tdt[d][k1] = Tdt[d];
      ddS.Tdt[d][k2] = 0;
    }
    ddS.TDt[k1] = Tdiff1;
    ddS.TDT += ddS.TDt[k1];
    ddS.TDt[k2] = 0;
    ddS.Tlife[k2] = 0;
    ddS.TDTnz--;
    for (d=ddN.DT; d<ddN.D; d++) {
      ddS.Tdt[d][k1] += ddS.Tdt[d][k2];
      ddS.Tdt[d][k2] = 0;
    }
  }

  /*
   *  reset Beta side stats
   */
  for (w=0; w<ddN.W; w++) {
    ddS.Nwt[w][k1] += ddS.Nwt[w][k2];
    ddS.Nwt[w][k2] = 0;
  }
  ddS.NWt[k1] += ddS.NWt[k2];
  ddS.NWt[k2] = 0;
  if ( ddP.PYbeta ) {
    Tdiff1 = 0;
    ddS.TWT -= ddS.TWt[k1] + ddS.TWt[k2];
    for (w=0; w<ddN.W; w++) {
      ddS.TwT[w] += Twt[w] - ddS.Twt[w][k1] - ddS.Twt[w][k2];
      Tdiff1 += Twt[w];
      ddS.Twt[w][k1] = Twt[w];
      ddS.Twt[w][k2] = 0;
    }
    ddS.TWt[k1] = Tdiff1;
    ddS.TWt[k2] = 0;
    ddS.TWT += ddS.TWt[k1];
  }
}

static void hca_reset_sparse(char *resstem, int restart, int firstdoc, int lastdoc) {
  int readOK = 0;
  char *restartfile;
  FILE *fpin;
  int i, t;
    /*
     *  now, deal with .sparse[][] and .sparseD[]
     */
    if ( restart ) {
      /*
       *  check if file exists
       */
      restartfile = yap_makename(resstem, ".ngs");
      if ( firstdoc!=0 )
	yap_quit("Cannot read '%s' from %d-th doc\n", restartfile, firstdoc);
      fpin = fopen(restartfile ,"rb"); 
      if ( !fpin ) {
	readOK = 0;
	yap_message("Cannot open file '%s', setting .sparse to default\n", restartfile);
	free(restartfile);
      } else {
	readOK = 1;
      }
    } else {
      /*
       * default on no restart is no sparsity
       */
      if ( ddS.sparse ) {
	memset((void*)ddS.sparse[0], 0U, sizeof(ddS.sparse[0][0])*
	       ddN.D*M_bitveclen());
	for (t=0; t<ddN.T; t++) 
	  ddS.sparseD[t] = 0;
      }
    }
    if ( restart && readOK ) {
      /*
       *  now read sparse bitvectors
         */
      int size;
      size = ddN.D*M_bitveclen();
      if ( fread(ddS.sparse[0], sizeof(ddS.sparse[0][0]), size, fpin) 
	   !=size )
	yap_sysquit("Cannot read bitvector from '%s'\n", restartfile);
      fclose(fpin);
      free(restartfile);
    }
    if ( restart && !readOK ) {
      /*
       *  nothing to read so make all zero topics sparse
       */
      memset((void*)ddS.sparse[0], 0, sizeof(ddS.sparse[0][0])*
	     ddN.D*M_bitveclen());
    }
    /*
     *  ensure non-zero topics are not sparse
     */
    for (i=0; i<ddN.DT && i<lastdoc; i++) {
      if ( ddD.NdT[i]<ddP.mindocsize )
	continue;
      for (t=0; t<ddN.T; t++) 
	if ( ddS.Ndt[i][t]>0 )
	  M_docsp_set(i,t);
    }
    /*
     *  fix up total vectors
     */
    for (t=0; t<ddN.T; t++) {
      ddS.sparseD[t] = 0;
    }
    for (i=0; i<ddN.DT && i<lastdoc; i++) {
      if ( ddD.NdT[i]<ddP.mindocsize )
	continue;
      for (t=0; t<ddN.T; t++) {
	if ( M_docsparse(i,t) )
	  ddS.sparseD[t]++;
      }
    }
}

static void hca_reset_UN(char *resstem, int restart, int firstdoc, int lastdoc) {
  int readOK = 0;
  char *restartfile;
  FILE *fpin;
  int i, t;
  /*
   *  deal with .UN
   */
  if ( restart ) {
    /*
     *  check if file exists
     */
    restartfile = yap_makename(resstem, ".UN");
    if ( firstdoc!=0 )
      yap_quit("Cannot read '%s' from %d-th doc\n", restartfile, firstdoc);
    fpin = fopen(restartfile ,"r"); 
    if ( !fpin ) {
      readOK = 0;
      yap_message("Cannot open file '%s', setting UN to default\n",
		  restartfile);
      free(restartfile);
    } else {
      fclose(fpin);
      readOK = 1;
    }
  }
  if ( restart && readOK ) {
    read_dvec(restartfile, lastdoc, ddS.UN);
    free(restartfile);
  } else {
    /*
     *   no restart OR no ".UN" file to restart with
     */
    double aveB = 0;
    double totA = 0;
    assert(ddP.NGbeta);
    for (t=0; t<ddN.T; t++) {
      aveB += ddP.NGbeta[t];
      totA += ddP.NGalpha[t];
    }
    aveB /= ddN.T;
    for (i=firstdoc; i<lastdoc; i++) {
      if ( ddS.NdT[i]==0 ) continue;
      ddS.UN[i] = ddS.NdT[i]/(totA+1)*aveB;
    }
  }
#ifdef NG_SCALESTATS
  {
    for (i=firstdoc; i<lastdoc; i++) {
      if ( ddS.NdT[i]==0 ) continue;
      for (t=0; t<ddN.T; t++) {
	if ( !ddS.sparse || M_docsparse(i,t) )
	  ddS.NGscalestats[t] += log(1.0+ddS.UN[i]/ddP.NGbeta[t]);
      }
      ddS.NGscalestats_recomp = 0;
    }
  }
#endif
}

static void hca_zero_N() {
  int i;
  /*
   *  reset data counts
   */
  memset((void*)ddS.NdT, 0, sizeof(ddS.NdT[0])*ddN.D);
  memset((void*)ddS.NWt, 0, sizeof(ddS.NWt[0])*ddN.T);
  for (i=0; i<ddN.D; i++)
    /*   ddS.Ndt not allocated monolithically  */
    memset((void*)ddS.Ndt[i], 0, sizeof(ddS.Ndt[0][0])*ddN.T);
  memset((void*)ddS.Nwt[0], 0, sizeof(ddS.Nwt[0][0])*ddN.W*ddN.T);
}

static void hca_zero_T() {
  int i;
  /*
   *  reset table count stats
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
}

static void hca_init_N(int firstdoc, int lastdoc) { 
  int i;
  for (i=firstdoc; i<lastdoc; i++) {
    int t, l;
    int usei = (i+ddN.DT) % ddN.DT;
    for (l=ddD.NdTcum[usei]; l<ddD.NdTcum[usei+1]; l++) {
      t = Z_t(ddS.z[l]);
      if ( ( ddP.bdk==NULL ) || Z_issetr(ddS.z[l]) ) {
	ddS.Nwt[ddD.w[l]][t]++;
	ddS.NWt[t]++;
      }
      ddS.Ndt[usei][t]++;
      ddS.NdT[usei]++;
    }
  }
}

static void hca_init_betaT(char *resstem, int restart) {
    int readOK = 0;
    int i, t;
    if ( restart ) {
      char *fname = yap_makename(resstem,".twt");
      /*  check if file is readable */
      if ( access(fname,R_OK)==0 ) {
	read_u16sparse(ddN.W,ddN.T,ddS.Twt,fname);
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
	readOK = 1;
      }
      free(fname);
    }
    if ( !readOK ) {
      /* no read, so force initialise */
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

static void hca_init_alphaT(char *resstem, int restart, int firstdoc, int lastdoc) {
   int readOK = 0;
   int i, t;
   if ( restart ) {
     char *fname = yap_makename(resstem,".tdt");
      /*  check if file is readable */
      if ( access(fname,R_OK)==0 ) {
	read_u16sparse(ddN.DT,ddN.T,ddS.Tdt,fname);
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
	readOK = 1;
      }
      free(fname);
    }
    if ( !readOK ) {
      /* no read, so force initialise */
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
  if ( ddS.UN )
    memset((void*)ddS.UN, 0, sizeof(ddS.UN[0])*ddN.D);
#ifdef NG_SCALESTATS
  if ( ddS.NGscalestats )
    memset((void*)ddS.NGscalestats, 0, sizeof(ddS.NGscalestats[0])*ddN.T);
#endif
  hca_zero_N();
  hca_zero_T();

  if ( zero )
    return;

  hca_init_N(firstdoc, lastdoc);

  if ( ddS.UN ) {
    if ( ddS.sparse ) 
      hca_reset_sparse(resstem, restart, firstdoc, lastdoc);
    hca_reset_UN(resstem, restart, firstdoc, lastdoc);
  }

  if ( ddP.PYbeta && ddP.phi==NULL ) 
    hca_init_betaT(resstem, restart);

  if ( ddP.PYalpha  ) 
    hca_init_alphaT(resstem, restart, firstdoc, lastdoc);
}

