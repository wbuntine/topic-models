/*
 * Probability utilities
 * Copyright (C) 2010-2011 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine@nicta.com.au
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
#include "data.h"
#include "hca.h"
#include "stats.h"


/*
 *   no fancy statistics here, just the empirical
 *   topic counts smoothed with alpha;
 *   in fact its meaningless statistically, just a
 *   useful diagnostic
 */
static void get_probs_alpha(double *vp) {
  int t, d, NWT=0;
  if ( ddS.Ndt ) {
    for (d=0; d<ddN.DT; d++)
      NWT += ddS.NdT[d];
    for (t=0; t<ddN.T; t++) {
      int tot=0;
      for (d=0; d<ddN.DT; d++)
        tot += ddS.Ndt[d][t];
      vp[t] = (ddP.alphapr[t]+tot)/(ddP.alphatot+NWT);
    }
  } else if ( ddP.alphapr ) {
    for (t=0; t<ddN.T; t++) {
      vp[t] = ddP.alphapr[t];
    }
  }
}

void get_probs(double *vp) {
  int zerod = 1;
  int t;
  double tot = 0;
  if ( ddP.PYalpha==H_None ) {
    get_probs_alpha(vp);
    return;
  } 
  for (t=0; t<ddN.T; t++) {
    /*
     *   all this trickery so only the first NULL topic
     *   gets the remainder probability
     */
    if ( ddP.PYalpha!=H_HPDD || ddS.TDt[t]>0 || zerod ) {
      tot += vp[t] = alphabasetopicprob(t);
      if (zerod) zerod = 0;
    } else
      vp[t] = 0;
  }
#ifndef NDEBUG
  if ( fabs(tot-1.0)>1e-4 ) {
    yap_message("get_probs() probs doesn't normalise, get %lf\n", tot);
  }
#endif
  for (t=0; t<ddN.T; t++) 
    vp[t] /= tot;
}

void yap_probs() {
  int t;
  int empty = 0;
  double ent = 0;
  double *vp = malloc(sizeof(*vp)*ddN.T);
  get_probs(vp);
  yap_message("probs = ");
  for (t=0; t<ddN.T; t++) 
    if ( vp[t]>0 ) {
      yap_message(" %lf", vp[t]);
      ent -= vp[t]*log(vp[t]);
    } else {
      empty++;
      yap_message(" -");
    }
  yap_message("\nfactor = %lf, empty = %d, ent = %lf\n", 
	      (ddP.PYalpha)?ddP.bpar:ddP.alphatot,
              empty, exp(ent));
  free(vp);
}

void print_probs(FILE *fp) {
  int t, num = 0;
  double *vp = malloc(sizeof(*vp)*ddN.T);
  get_probs(vp);
  fprintf(fp, "factor = %lf\nprobs = ", 
          (ddP.PYalpha)?ddP.bpar:ddP.alphatot);
  for (t=0; t<ddN.T; t++) 
    if ( vp[t]>0 ) {
      fprintf(fp, " %lf", vp[t]);
      num++;
    } else
      fprintf(fp, " -");
  fprintf(fp, "\n");
  fprintf(fp, "# topics: %d\n", num);
  free(vp);
}

