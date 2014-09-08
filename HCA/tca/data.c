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
#include "tca.h"
#include "data.h"
#include "stats.h"
#include "probs.h"

char **read_vocab(char *infile, int W, int len);
void free_vocab(char **ctmp);


/*
 *  allocation of various data statistics, matrices, vectors
 */
void data_alloc() {
  int i;
  if ( ddN.DT<5 )
    yap_quit("Only %d training examples\n", ddN.DT);

  ddD.N_dT  = u16vec(ddN.D);
  ddD.N_dTcum=u32vec(ddN.D+1);
  if ( !ddD.N_dT )
    yap_quit("Cannot allocate memory for vecs/matrices\n");
  for(i=0;i<ddN.N;i++)
    ddD.N_dT[ddD.d[i]]++;
  ddD.N_dTcum[0] = 0;
  for(i=1;i<=ddN.D;i++)
    ddD.N_dTcum[i]=ddD.N_dTcum[i-1]+ddD.N_dT[i-1];
  ddD.N_dTmax = ddD.N_dT[0];
  for(i=1;i<ddN.D;i++)
    if ( ddD.N_dTmax<ddD.N_dT[i] )
      ddD.N_dTmax = ddD.N_dT[i];
  assert(ddD.N_dTcum[ddN.D]==ddN.N);
  ddN.NT = ddD.N_dTcum[ddN.DT];
  ddN.tokens = NULL;
}

void data_vocab(char *stem) {
    char *wname = yap_makename(stem, ".tokens");
    ddN.tokens = read_vocab(wname,ddN.W,50);
    free(wname);
}

void data_free() {
  /*
   *  free
   */
  free(ddD.w);
  free(ddD.d);
  free(ddD.e);
  free(ddD.j);
  free(ddD.esize);
  free(ddD.N_dT);
  free(ddD.N_dTcum);
  if ( ddN.tokens )
    free_vocab(ddN.tokens);
}

int data_docsize() {
  int i;
  int maxd = 100;
  for (i=0; i<ddN.D; i++)
    if ( maxd<ddD.N_dT[i] )
      maxd = ddD.N_dT[i];
  return maxd;
}
  
void data_read_epoch(char *stem) {
  char *wname;
  FILE  *fp;
  unsigned ein, d;
  int i, range;
  wname = yap_makename(stem, ".epoch");
  fp = fopen(wname,"r");
  if ( !fp )
    yap_sysquit("Cannot open epoch file '%s'\n", wname);
  if ( fscanf(fp,"%u ", &ein) != 1 )
    yap_sysquit("Cannot read dimension from '%s'\n", wname);
  ddN.E = ein;
  ddD.e = u16vec(ddN.D);
  ddD.j = u16vec(ddN.D);
  ddD.esize = u32vec(ddN.E*2);
  for (d=0, i=0; i<ein; i++) {
    int j;
    if ( fscanf(fp," %u", &range) != 1 )
      yap_sysquit("Cannot read %n-th epoch range from '%s'\n", i+1, wname);
    for (j=0; j<range; j++) {
      ddD.e[d] = i;
      ddD.j[d] = j;
      d++;
    }
    ddD.esize[i] = range;
  }
  if (d!=ddN.D && d!=ddN.D-ddN.TEST) 
	yap_quit("Total count of docs in epochs, %d, not equal to docs %d/%d\n",
		 d, ddN.D, ddN.D-ddN.TEST);
  if ( d==ddN.D-ddN.TEST && d<ddN.D) {
    /*    more to read */
    if ( fscanf(fp,"%u ", &ein) != 1 )
      yap_sysquit("Cannot read dimension from '%s'\n", wname);
    if ( ddN.E != ein )
      yap_quit("Training epochs not equal to test epochs\n");
    for (i=0; i<ein; i++) {
      int j;
      if ( fscanf(fp," %u", &range) != 1 )
	yap_sysquit("Cannot read %n-th epoch range from '%s'\n", i+1, wname);
      for (j=0; j<range; j++) {
	ddD.e[d] = i;
	ddD.j[d] = j;
	d++;
      }
      ddD.esize[ein+i] = range;
    }
  }
  if (d!=ddN.D) 
	yap_quit("Total count of docs in epochs, %d, not equal to docs %d\n",
		 d, ddN.D);
  fclose(fp);
  free(wname); 
  if ( verbose>1 )
	yap_message("Read %d epochs\n", ddN.E);
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
  yap_message("E     = %d\n", ddN.E);
  yap_message("TRAIN   = %d\n", ddN.DT);
  yap_message("TEST    = %d\n", ddN.TEST);
  yap_message("T     = %d\n", ddN.T);
  yap_message("ITER  = %d\n", ITER);
}

extern void tca_write_z(char *resstem);
/*
 *  save various files suitable for a restart or for
 *  report handling, and also save the parameters
 *        STEM.n_dt = n_dt matrix in standard sparse format
 *        STEM.c_dt = c_dt matrix in standard sparse format
 *        STEM.cp_et = cp_et matrix in standard sparse format
 *        STEM.m_evt = m_evt matrix in standard sparse format
 *        STEM.s_evt = s_evt matrix in standard sparse format
 *        STEM.zt  = z vector for all training docs sequentially,
 *                  one value per line
 *        STEM.par = various parameters in readable form
 */ 
void data_checkpoint(char *resstem, char *stem, int ITER) {
    char *fname;
#if 0 
    fname = yap_makename(resstem,".ndt");
    write_u16sparse(ddN.DT,ddN.T,ddS.n_dt,fname);
    free(fname);
#endif
    fname = yap_makename(resstem,".cdt");
    write_u16sparseco(ddN.DT,ddN.T,ddS.c_dt,fname,1);
    free(fname);
    fname = yap_makename(resstem,".cpet");
    write_u32sparseco(ddN.E,ddN.T,ddS.cp_et,fname,1);
    free(fname);
#if 1
    fname = yap_makename(resstem,".mevt");
    write_u32sparsetri(ddN.E,ddN.W,ddN.T,ddS.m_evt,fname,0);
    free(fname);
#endif
    fname = yap_makename(resstem,".sevt");
    write_u32sparsetri(ddN.E,ddN.W,ddN.T,ddS.s_evt,fname,1);
    free(fname);
    tca_write_z(resstem);

    fname = yap_makename(resstem,".par");
    FILE *fp = fopen(fname, "w");
    if ( !fp )
      yap_sysquit("Cannot open output '%s' file:", fname);
    fprintf(fp, "stem = %s\n", stem);
    fprintf(fp, "N = %d\n", ddN.N);
    fprintf(fp, "NT = %d\n", ddN.NT);
    fprintf(fp, "W = %d\n", ddN.W);
    fprintf(fp, "D = %d\n", ddN.D);
    fprintf(fp, "E = %d\n", ddN.E);
    fprintf(fp, "TRAIN = %d\n", ddN.DT);
    fprintf(fp, "TEST = %d\n", ddN.TEST);
    fprintf(fp, "T = %d\n", ddN.T);
    fprintf(fp, "ITER = %d\n", ITER);
    pctl_print(fp);
    print_probs(fp);
    fclose(fp);
    free(fname);
}
