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

int data_df(char *stem, uint32_t *dfvec) {
  char buf[1000];
  char *wname = yap_makename(stem, ".srcpar");
  FILE *fp;
  int n_df;
  int i;
  /*
   *  check .srcpar file exists
   */
  fp = fopen(wname,"r");
  if ( !fp ) 
    yap_quit("Parameter file '%s' doesn't exist\n", wname);
  fclose(fp);
  free(wname);
  /*  read dfdocs */
  {
    char *p = readsrcpar(stem,"dfdocs",50);
    if ( p ) {
      n_df = atoi(p);
      free(p);
    } else
      n_df = ddN.D;
  }
  /*  read dfs */
  wname = yap_makename(stem, ".words");
  fp = fopen(wname ,"r"); 
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for read\n", wname);
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
    dfvec[i] = df;
  }
  fclose(fp);
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
      write_dvec(fname, ddN.D*ddN.T, ddS.UN);
      free(fname);
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
