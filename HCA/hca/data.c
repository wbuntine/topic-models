/*
 * Various data structure read/write/report routines.
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
 *     Defined in "data.h"
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "hca.h"
#include "data.h"
#include "stats.h"
#include "probs.h"

/*
 *  allocation of various data statistics, matrices, vectors
 */
void data_alloc() {
  int i;
  ddD.c = NULL;
  if ( ddP.phi==NULL && ddN.DT<5 )
    yap_quit("Only %d training examples\n", ddN.DT);

  ddD.NdT  = u16vec(ddN.D);
  ddD.df = NULL;
  ddD.n_df = 0;
  ddD.NdTcum=u32vec(ddN.D+1);
  if ( !ddD.NdT )
    yap_quit("Cannot allocate memory for vecs/matrices\n");
  for(i=0;i<ddN.N;i++)
    ddD.NdT[ddD.d[i]]++;
  ddD.NdTcum[0] = 0;
  for(i=1;i<=ddN.D;i++)
    ddD.NdTcum[i]=ddD.NdTcum[i-1]+ddD.NdT[i-1];
  ddD.NdTmax = ddD.NdT[0];
  for(i=1;i<ddN.D;i++)
    if ( ddD.NdTmax<ddD.NdT[i] )
      ddD.NdTmax = ddD.NdT[i];
  assert(ddD.NdTcum[ddN.D]==ddN.N);
  ddN.NT = ddD.NdTcum[ddN.DT];
  ddN.tokens = NULL;
}

void data_vocab(char *stem) {
    char *wname = yap_makename(stem, ".tokens");
    ddN.tokens = read_vocab(wname,0,ddN.W,50);
    free(wname);
}

/*
 *   returns 0 if no file read;
 *   quits on file format error;
 *   otherwise fills dfvec[] and returns no.. docs used
 */
int data_df(char *stem) {
  char buf[1000];
  char *wname;
  FILE *fp;
  int n_df;
  int i;
  
  if ( ddD.df )
    return ddD.n_df;
  ddD.df = u32vec(ddN.W);

  /*
   *  check .srcpar file exists
   */
  wname = yap_makename(stem, ".srcpar");
  fp = fopen(wname,"r");
  if ( fp ) { 
    fclose(fp);
    /*  it does so read dfdocs */
    {
      char *p = readsrcpar(stem,"dfdocs",50);
      if ( p ) {
	n_df = atoi(p);
	free(p);
      } else
	n_df = ddN.DT;
    }
  } else {
    /*  it does not, so set default */
    n_df = ddN.DT;
  }
  free(wname);
  
  /*  
   *  first try read dfs from "stem.words"
   */
  wname = yap_makename(stem, ".df");
  fp = fopen(wname ,"r"); 
  if ( fp ) {
    for (i = 0; i < ddN.W; i++) {
      int sl;
      unsigned df;
      if ( fgets(&buf[0],sizeof(buf)-1,fp)==NULL )
        yap_sysquit("Cannot read line %d from '%s'\n", i+1, wname);
      sl = strlen(&buf[0]);
      if ( ! iscntrl(buf[sl-1]) )
        /*   line too long  */
        yap_quit("Cannot parse line %d from '%s', too long\n", i, wname);
      if ( sscanf(&buf[0],"%u", &df) != 1 )
        yap_quit("Cannot parse line %d from '%s', no df\n", i, wname);
      ddD.df[i] = df;
    }
    fclose(fp);
  } else {
    /*  
     *  second try read dfs from "stem.df" 
     */
    free(wname);
    wname = yap_makename(stem, ".words");
    fp = fopen(wname ,"r"); 
    if ( fp ) { 
      for (i = 0; i < ddN.W; i++) {
	int sl;
	unsigned df;
	if ( fgets(&buf[0],sizeof(buf)-1,fp)==NULL )
	  yap_sysquit("Cannot read line %d from '%s'\n", i+1, wname);
	sl = strlen(&buf[0]);
	if ( ! iscntrl(buf[sl-1]) )
	  /*   line too long  */
	  yap_quit("Cannot parse line %d from '%s', too long\n", i, wname);
	if ( sscanf(&buf[0],"%*u %*s %*x %*u %u ", &df) != 1 )
	  yap_quit("Cannot parse line %d from '%s', no df\n", i, wname);
	ddD.df[i] = df;
      }
      fclose(fp);
    } else {
      uint32_t *flag = u32vec(ddN.W);
      /* 
       *   else compute from training data
       */
      for (i = 0; i < ddN.DT; i++) {
	int l;
	for (l=ddD.NdTcum[i]; l<ddD.NdTcum[i+1]; l++) {
	  int w = ddD.w[l];
	  if ( flag[w]==0 ) {
	    ddD.df[w] ++;
	    flag[w]++;
	  }
	}
	// second pass to zero flag[]
	for (l=ddD.NdTcum[i]; l<ddD.NdTcum[i+1]; l++) {
	  flag[ddD.w[l]]=0;
	}
      }
      free(flag);
    }
  }
  free(wname);
  return n_df;
}

/*
 *  read a file of ddN.D entries with 0,...,C-1
 */
void data_class(char *stem) {
    char *cfile = yap_makename(stem, ".class");
    int i, c;
    int maxc = 0;
    ddD.c = u16vec(ddN.D);
    FILE *fp = fopen(cfile ,"r"); 
    if ( !fp ) 
      yap_sysquit( "Cannot open file '%s' for read\n", cfile);
    for (i = 0; i < ddN.D; i++) {
      if ( fscanf(fp,"%d",&c) ) {
	ddD.c[i] = c;
	if ( c>maxc )
	  maxc = c;
      } else
	break;
    }
    if ( ferror(fp) )
      yap_sysquit("Error on reading file '%s' ", cfile);
    fclose(fp);
    ddN.C = maxc+1;
    free(cfile);
}

void data_free() {
  /*
   *  free
   */
  free(ddD.w);
  free(ddD.d);
  free(ddD.NdT);
  free(ddD.NdTcum);
  if ( ddN.tokens )
	free_vocab(ddN.tokens);
  if ( ddD.c )
    free(ddD.c);
  if ( ddD.df )
    free(ddD.df);
}

int data_docsize() {
  int i;
  int maxd = 100;
  for (i=0; i<ddN.D; i++)
    if ( maxd<ddD.NdT[i] )
      maxd = ddD.NdT[i];
  return maxd;
}


void data_report(int ITER, int seed) {
  if ( memallocd>5e8 )
    yap_message("mem   = %.3f (GByte)\n", 1e-9*memallocd);
  else
    yap_message("mem   = %.1f (MByte)\n", 1e-6*memallocd);
  yap_message("seed  = %d\n", seed);
  yap_message("N     = %d\n", ddN.N);
  yap_message("W     = %d\n", ddN.W);
  yap_message("D     = %d\n", ddN.D);
  yap_message("TRAIN   = %d\n", ddN.DT);
  yap_message("TEST    = %d\n", ddN.TEST);
  yap_message("T     = %d\n", ddN.T);
  yap_message("ITER  = %d\n", ITER);
}

extern void hca_write_z(char *resstem);
/*
 *  save various files suitable for a restart or for
 *  report handling, and also save the parameters
 *        STEM.ndt = Ndt matrix in standard sparse format
 *        STEM.tdt = Tdt matrix in standard sparse format
 *        STEM.nwt = Nwt matrix in standard sparse format
 *        STEM.twt = Twt matrix in standard sparse format
 *        STEM.zt  = z vector for all training docs sequentially,
 *                  one value per line
 *        STEM.UN  = U latent var for H_NG
 *        STEM.ngs = sparsity bit vectors for H_NG
 *        STEM.par = various parameters in readable form
 */
void data_checkpoint(char *resstem, char *stem, int ITER) {
    char *fname;

    fname = yap_makename(resstem,".ndt");
    write_u16sparse(ddN.DT,ddN.T,ddS.Ndt,fname);
    free(fname);
    if ( ddP.PYalpha ) {
      fname = yap_makename(resstem,".tdt");
      write_u16sparse(ddN.DT,ddN.T,ddS.Tdt,fname);
      free(fname);
    }
    if ( ddP.phi==NULL ) {
      fname = yap_makename(resstem,".nwt");
      write_u32sparse(ddN.W,ddN.T,ddS.Nwt,fname);
      free(fname);
      if ( ddP.PYbeta ) {
	fname = yap_makename(resstem,".twt");
	write_u16sparse(ddN.W,ddN.T,ddS.Twt,fname);
	free(fname);
      }
    }
    hca_write_z(resstem);
    if ( ddS.UN ) {
      fname = yap_makename(resstem,".UN");
      write_dvec(fname, ddN.D, ddS.UN);
      free(fname);
      if ( ddS.sparse ) {
	/*
	 *    write as binary
	 */
	FILE *fpout=NULL;
	int size;
	fname = yap_makename(resstem,".ngs");
	fpout = fopen(fname,"wb");
	if ( !fpout ) 
	  yap_sysquit("Cannot open file '%s' for write in data_checkpoint()\n", 
		      fname);
	size = ddN.D*M_bitveclen();
	if ( fwrite(ddS.sparse[0], sizeof(ddS.sparse[0][0]), size, fpout) 
	     !=size )
	  yap_sysquit("Cannot write bitvector to '%s' in data_checkpoint()\n", fname);
	fclose(fpout);
	free(fname);
      }
    }

    fname = yap_makename(resstem,".par");
    FILE *fp = fopen(fname, "w");
    if ( !fp )
      yap_sysquit("Cannot open output '%s' file:", fname);
    fprintf(fp, "stem = %s\n", stem);
    fprintf(fp, "N = %d\n", ddN.N);
    fprintf(fp, "NT = %d\n", ddN.NT);
    fprintf(fp, "W = %d\n", ddN.W);
    fprintf(fp, "D = %d\n", ddN.D);
    fprintf(fp, "TRAIN = %d\n", ddN.DT);
    fprintf(fp, "TEST = %d\n", ddN.TEST);
    fprintf(fp, "T = %d\n", ddN.T);
    fprintf(fp, "ITER = %d\n", ITER);

    pctl_print(fp);
    print_probs(fp);
    fclose(fp);
    free(fname);
}
