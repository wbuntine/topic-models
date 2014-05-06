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
  free(fname);
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
  free(fname);
}

