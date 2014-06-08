/*
 * Saving document topic probs
 * Copyright (C) 2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include "hca.h"
#include "stats.h"
#include "util.h" 
#include "yap.h" 
#include "diag.h" 

extern int verbose;

/*
 *    report sparsity for words in map
 */
void tprob_init() {
  if ( ddP.probiter>0 )
    ddG.prob = fmat(ddN.DT,ddN.T);
  if ( ddP.tprobiter>0 )
    ddG.tprob = fmat(ddN.D-ddN.DT,ddN.T);
  ddG.doprob = 0;
  ddG.didprob = 0;
  ddG.didtprob = 0;
}

void tprob_null() {
  ddG.prob = NULL;
  ddG.doprob = 0;
  ddG.didprob = 0;
  ddG.tprob = NULL;
  ddG.didtprob = 0;
}

void tprob_free() {
  if ( ddG.tprob ) {
    free(ddG.tprob[0]);free(ddG.tprob);
  }
  if ( ddG.prob ) {
    free(ddG.prob[0]);free(ddG.prob);
  }
  tprob_null();
}

void tprob_report(char *resstem, double epsilon) {
  FILE *fp;
  int n, k;
  char *fname;
  if ( ddP.teststem )
    fname = yap_makename(ddP.teststem,".testprob");
  else
    fname = yap_makename(resstem,".testprob");
  fp = fopen(fname,"w");
  if ( !fp )
    yap_sysquit("Cannot open doc-topic output file '%s'\n", fname);
  for (n=ddN.DT; n<ddN.D; n++) {
    float tot = 0;
    fprintf(fp,"%d:", n);
    if ( ddD.c ) {
      if ( ddN.C==2 ) {
	if ( ddD.c[n]==0 )
	  fprintf(fp," -1");
	else
	  fprintf(fp," +1");
      } else
	fprintf(fp," %d", ddD.c[n]);
    }
    for (k=0; k<ddN.T; k++) {
      double val = ddG.tprob[n-ddN.DT][k]/ddG.didtprob/ddD.NdT[n];
      if ( val>epsilon ) 
	fprintf(fp," %d:%f", k, val);
      tot += val;
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  if ( verbose )
    yap_message("Written test theta estimate as tagged sparse matrix to '%s'\n",
                fname);
  free(fname);
}

void tprob_load(char *resstem) {
  FILE *fp;
  int n;
  char *buf;
  int bufsize = 20*ddN.T+10;
  char *fname;
  buf = malloc(bufsize+1);
  ddP.theta = fmat(ddN.D-ddN.DT,ddN.T);
  if ( !ddP.theta || !buf)
    yap_quit("Out of memory in tprob_load()\n");
  if ( ddP.teststem )
    fname = yap_makename(ddP.teststem,".testprob");
  else
    fname = yap_makename(resstem,".testprob");
  fp = fopen(fname,"r");
  if ( !fp )
    yap_sysquit("Cannot open doc-topic input file '%s'\n", fname);
  for (n=0; n<ddN.D-ddN.DT; n++) {
    float f;
    int nin;
    char *bufptr;
    if ( fgets(&buf[0],bufsize,fp)==NULL )
      yap_sysquit("Cannot read line %d from '%s'\n", n+1, fname);
    if ( ! iscntrl(buf[strlen(&buf[0])-1]) )
      /*   line too long  */
      yap_quit("Cannot parse line %d from '%s', too long\n", n+1, fname);
    bufptr = &buf[0];
    if ( sscanf(bufptr,"%d", &nin)!=1 ) 
      yap_quit("Cannot read line %d from input file '%s'\n", n+1, fname);
    if ( n!=nin ) 
      yap_quit("Bad doc %d!=%d from input file '%s'\n", n, nin, fname);
    bufptr = strchr(bufptr,':');
    if ( !bufptr ) 
      yap_quit("Bad line %d from input file '%s'\n", n+1, fname);
    bufptr = strrchr(bufptr,' ');
    if ( !bufptr ) 
      yap_quit("Bad line %d from input file '%s'\n", n+1, fname);
    bufptr ++;
    while ( sscanf(bufptr,"%d:%f", &nin, &f)==2 ) {
      ddP.theta[n][nin] = f;
      bufptr = strchr(bufptr,' ');
      if ( !bufptr ) 
        break;
      bufptr = strrchr(bufptr,' ');
      if ( !bufptr ) 
        yap_quit("Bad line %d from input file '%s'\n", n+1, fname);
      bufptr ++;
    }
  }
  fclose(fp);
  free(fname);
  free(buf);
}



void prob_report(char *resstem, double epsilon) {
  FILE *fp;
  int n, k;
  char *fname = yap_makename(resstem,".theta");
  fp = fopen(fname,"w");
  if ( !fp )
    yap_sysquit("Cannot open doc-topic output file '%s'\n", fname);
  for (n=0; n<ddN.DT; n++) {
    float tot = 0;
    fprintf(fp,"%d:", n);
    if ( ddD.c ) {
      if ( ddN.C==2 ) {
	if ( ddD.c[n]==0 )
	  fprintf(fp," -1");
	else
	  fprintf(fp," +1");
      } else
	fprintf(fp," %d", ddD.c[n]);
    }
    for (k=0; k<ddN.T; k++) {
      double val = ddG.prob[n][k]/ddG.didprob/ddD.NdT[n];
      if ( val>epsilon ) 
	fprintf(fp," %d:%f", k, val);
      tot += val;
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  if ( verbose )
    yap_message("Written training theta estimate as tagged sparse matrix to '%s'\n",
                fname);
  free(fname);
}

